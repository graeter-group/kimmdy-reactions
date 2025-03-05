from kimmdy.topology.atomic import Atom
from kimmdy.constants import ATOMTYPE_BONDORDER_FLAT
from kimmdy.recipe import (
    Break,
    Bind,
    Relax,
    Recipe,
    RecipeCollection,
    RecipeStep,
)
from kimmdy.plugins import ReactionPlugin
import MDAnalysis as mda
import random as rng


def find_radical(atoms: list[Atom]):
    for atom in atoms:
        if atom.is_radical:
            return atom
        bo = ATOMTYPE_BONDORDER_FLAT.get(atom.type)
        if bo and bo > len(atom.bound_to_nrs):
            return atom

    return None


class NaiveHAT(ReactionPlugin):
    """Naive HAT reaction, selects all neighboring hydrogens and assigns random rates."""

    def get_recipe_collection(self, files) -> RecipeCollection:
        logger = files.logger
        logger.info("Starting naive HAT reaction")
        top = self.runmng.top

        gro = files.input["gro"]
        trr = files.input["trr"]
        u = mda.Universe(str(gro), str(trr))

        if not top.radicals:
            radical = find_radical(list(top.atoms.values()))
            if radical:
                top.radicals[radical.nr] = radical

        t1 = u.trajectory[0].time
        t2 = u.trajectory[-1].time
        full_timespan = (t1, t2)
        recipes = []
        if len(top.radicals) > 0:
            for rad in top.radicals.values():
                hs = []
                froms = []
                for nr in rad.bound_to_nrs:
                    atom = top.atoms[nr]
                    for nr2 in atom.bound_to_nrs:
                        atom2 = top.atoms[nr2]
                        if atom2.type.startswith("H"):
                            froms.append(atom.nr)
                            hs.append(atom2.nr)
                if len(hs) == 0:
                    continue
                r = rad.nr
                for h, f in zip(hs, froms):
                    rate = rng.random()
                    logger.debug(f"radical: {rad}")
                    logger.debug(f"h: {top.atoms[h]}")
                    logger.debug(f"from: {top.atoms[f]}")
                    # int(x) - 1 to be zero based because h,f,r are from topology
                    recipe = Recipe(
                        recipe_steps=[
                            Break(atom_id_1=f, atom_id_2=h),
                            Bind(atom_id_1=h, atom_id_2=r),
                            Relax(),
                        ],
                        rates=[rate],
                        timespans=[full_timespan],
                    )
                    recipes.append(recipe)

            if len(recipes) == 0:
                # empty recipe that does nothing as a fallback
                recipe = Recipe(
                    recipe_steps=[RecipeStep()],
                    rates=[1],
                    timespans=[full_timespan],
                )
                recipes.append(recipe)

            return RecipeCollection(recipes)

        return RecipeCollection([])
