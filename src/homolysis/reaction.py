from pathlib import Path

import numpy as np
from kimmdy.parsing import read_edissoc
from kimmdy.plugin_utils import (
    bondstats_from_csv,
    bondstats_to_csv,
    calculate_beta,
    calculate_bondstats,
    calculate_harmonic_forces,
    get_edissoc_from_atomnames,
    harmonic_rates_from_forces,
    harmonic_transition_rate,
    morse_rates_from_forces,
    morse_transition_rate,
)
from kimmdy.plugins import ReactionPlugin
from kimmdy.recipe import Break, Recipe, RecipeCollection, Relax
from kimmdy.tasks import TaskFiles


class Homolysis(ReactionPlugin):
    """Homolytic bond breaking leading to 2 radicals.
    Implementation for time-varying rates
    """

    def get_recipe_collection(self, files: TaskFiles):
        logger = files.logger
        self.logger = logger
        logger.debug("Getting recipe for reaction: homolysis")

        xtc = files.input["xtc"]
        trr = files.input["trr"]
        if xtc is not None:
            trj = xtc
        elif trr is not None:
            trj = trr
        else:
            m = "No xtc trr file found. The homolysis plugin requires a trajectory file."
            logger.error(m)
            raise ValueError(m)

        self.bondstatsfile = trj.with_name(".kimmdy.bondstats")
        md_instance = trj.stem
        timings = self.runmng.timeinfos.get(md_instance)
        if timings is None:
            m = f"No timings from mdp file found for {md_instance}"
            logger.error(m)
            raise ValueError(m)
        self.xtc_trr_ratio = timings.trr_nst / timings.xtc_nst
        self.dt = timings.dt
        top = self.runmng.top

        did_read_bondstats = self.use_cached_bondstats()
        if not did_read_bondstats:
            plumed_out = files.input["plumed_out"]
            plumed_in = files.input["plumed"]
            if plumed_out is None or plumed_in is None:
                m = f"External force not specified but no plumed file found"
                logger.error(m)
                raise ValueError(m)
            self.bondstats = calculate_bondstats(
                top=top,
                plumed_in=plumed_in,
                plumed_out=plumed_out,
                dt=self.config.dt_distances,
                edissoc_dat=files.input["edissoc.dat"],
            )
            self.cache_bondstats()

        path_edissoc = files.input["edissoc.dat"]
        if path_edissoc is None:
            m = "edissoc.dat file not found. The homolysis plugin requires a file with dissociation energies defined in the kimmdy config as `edissoc`."
            logger.error(m)
            raise ValueError(m)
        edissoc = read_edissoc(path_edissoc)
        t0 = 0.0
        t1 = timings.t_max

        recipes = []
        for bondkey, stats in self.bondstats.items():
            ai = top.atoms[bondkey[0]]
            aj = top.atoms[bondkey[1]]
            typekey = (ai.type, aj.type)
            bondtype = top.ff.bondtypes.get(typekey)
            if self.config.check_bound:
                if not aj.nr in ai.bound_to_nrs:
                    continue
            if bondtype is None:
                typekey = (aj.type, ai.type)
                bondtype = top.ff.bondtypes.get(typekey)
            if bondtype is None or bondtype.c0 is None or bondtype.c1 is None:
                m = f"Bond type {typekey} for atoms {ai.nr} and {aj.nr} not found in force field."
                logger.error(m)
                raise ValueError(m)
            b0 = float(bondtype.c0)
            kb = float(bondtype.c1)

            if self.config.b0_overwrite != 0.0:
                b0 = self.config.b0_overwrite

            if edissoc is not None:
                edis = get_edissoc_from_atomnames(
                    atomnames=[ai.atom, aj.atom], edissoc=edissoc, residue=ai.residue
                )
            else:
                m = f"Could not find dissociation energy for atoms with ids {ai.nr} and {aj.nr}. Using default."
                logger.debug(m)
                edis = 500

            # averaging forces works here because we typically have
            # one conformational state per calculation
            dist = stats["mean_d"]
            if self.config.use_morse:
                force = stats["mean_f"]
            else:
                # for harmonic bonds, we use the mean harmonic force
                # which is the force at the mean distance
                # since the harmonic force is linear
                maybe_force = stats.get("harmonic_f")
                if maybe_force is None:
                    force: float = calculate_harmonic_forces(
                        ds=np.array([dist]), b0=b0, k=kb
                    )[0]
                else:
                    force = maybe_force

            beta = calculate_beta(kb=kb, edis=edis)

            if self.config.b0_overwrite != 0.0:
                # if b0 is overwritten, we need to recalculate the force
                # and don't use it from the bondstats
                if self.config.use_morse:
                    ks, fs = morse_transition_rate(
                        r_curr=[dist],
                        r_0=b0,
                        dissociation_energy=edis,
                        k_f=kb,
                        frequency_factor=self.config.arrhenius_equation.frequency_factor,
                        temperature=self.config.arrhenius_equation.temperature,
                    )
                else:
                    ks, fs = harmonic_transition_rate(
                        r_curr=[dist],
                        r_0=b0,
                        k_f=kb,
                        dissociation_energy=edis,
                        frequency_factor=self.config.arrhenius_equation.frequency_factor,
                    )
            else:
                if self.config.f0_overwrite != 0.0:
                    force = force - self.config.f0_overwrite

                # set negative average forces to 0
                force = max(force, 0)

                fs = np.asarray([force])
                if self.config.use_morse:
                    k = morse_rates_from_forces(
                        fs=fs,
                        b0=b0,
                        edis=edis,
                        beta=beta,
                        frequency_factor=self.config.arrhenius_equation.frequency_factor,
                        temperature=self.config.arrhenius_equation.temperature,
                    )
                else:
                    k = harmonic_rates_from_forces(
                        fs=fs,
                        edis=edis,
                        kb=kb,
                        frequency_factor=self.config.arrhenius_equation.frequency_factor,
                        temperature=self.config.arrhenius_equation.temperature,
                    )
                ks = list(k)

            recipes.append(
                Recipe(
                    recipe_steps=[
                        Break(atom_id_1=bondkey[0], atom_id_2=bondkey[1]),
                        Relax(),
                    ],
                    rates=[*ks],
                    timespans=[(t0, t1)],
                )
            )

        return RecipeCollection(recipes)

    def use_cached_bondstats(self) -> bool:
        if self.config.recompute_bondstats:
            return False
        if not Path(self.bondstatsfile).exists():
            m = f"bondstatsfile {self.bondstatsfile} does not exist. Not using cached bondstats."
            self.logger.info(m)
            return False

        m = f"bondstatsfile {self.bondstatsfile} found. Using cached bondstats."
        self.logger.info(m)
        self.bondstats = bondstats_from_csv(self.bondstatsfile)
        return True

    def cache_bondstats(self) -> None:
        bondstats_to_csv(self.bondstats, self.bondstatsfile)
