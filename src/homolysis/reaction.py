from kimmdy.recipe import (
    Break,
    Relax,
    Recipe,
    RecipeCollection,
)
from kimmdy.plugins import ReactionPlugin
from kimmdy.tasks import TaskFiles
from kimmdy.utils import (
    morse_transition_rate,
    get_atomnrs_from_plumedid,
    get_atominfo_from_atomnrs,
    get_bondprm_from_atomtypes,
    get_edissoc_from_atomnames,
)
from kimmdy.parsing import (
    read_top,
    read_plumed,
    read_distances_dat,
    read_edissoc,
)


class Homolysis(ReactionPlugin):
    """Homolytic bond breaking leading to 2 radicals.
    Implementation for time-varying rates
    """

    def get_recipe_collection(self, files: TaskFiles):
        logger = files.logger
        logger.debug("Getting recipe for reaction: homolysis")

        # Initialization of filepaths
        files.input["itp"] = self.config.itp
        files.input["edis"] = self.config.edis
        frequency_factor = self.config.arrhenius_equation.frequency_factor
        temperature = self.config.arrhenius_equation.temperature

        # Initialization of objects from files
        distances = read_distances_dat(files.input["plumed_out"])
        plumed = read_plumed(files.input["plumed"])
        top = self.runmng.top
        ffbonded = read_top(files.input["itp"])
        edissoc = read_edissoc(files.input["edis"])

        # handle already averaged distances
        t0, t1 = distances["time"][0], distances["time"][-1]
        if len(distances["time"]) == 1:
            logger.debug(
                "Plumed output contains single line, assuming averaged distance."
            )
            t0 = 0.0
            if t1 < 0.0001:
                logger.warning(
                    "First and last time in plumed output = 0.0"
                    "Check distance input.\n"
                    "If input is averaged already, set time to end time in ps."
                )
        if abs(t1 - t0) < 0.0001:
            logger.warning(
                "First and last time in plumed output equal!" "Check distance input."
            )

        recipes = []
        e_dis_d = {}
        for plumedid, dists in distances.items():
            if plumedid == "time":
                continue
            atomnrs = get_atomnrs_from_plumedid(plumedid, plumed)
            if self.config.check_bound:
                if not atomnrs[1] in top.atoms[atomnrs[0]].bound_to_nrs:
                    continue
            # get from plumedid to b0 and kb of the bond via atomtypes
            atomtypes, atomnames = get_atominfo_from_atomnrs(atomnrs, top)
            b0, kb = get_bondprm_from_atomtypes(atomtypes, ffbonded)

            residue = top.atoms[atomnrs[0]].residue

            if not (e_dis := e_dis_d.get("".join(atomnames + [residue]))):
                e_dis = get_edissoc_from_atomnames(atomnames, edissoc, residue)
                e_dis_d["".join(atomnames + [residue])] = e_dis

            k_avg, _ = morse_transition_rate(
                [sum(dists) / len(dists)], b0, e_dis, kb, frequency_factor, temperature
            )
            # averaging distances works here because we typically have
            # one conformational state per calculation

            recipes.append(
                Recipe(
                    recipe_steps=[
                        Break(atom_id_1=atomnrs[0], atom_id_2=atomnrs[1]),
                        Relax(),
                    ],
                    rates=[*k_avg],
                    timespans=[(t0, t1)],
                )
            )

        return RecipeCollection(recipes)
