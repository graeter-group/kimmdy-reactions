from collections import Counter
from pathlib import Path

import numpy as np
from kimmdy.parsing import read_edissoc
from kimmdy.plugin_utils import (
    atom_nrs_from_selection,
    bondstats_from_csv,
    bondstats_to_csv,
    calculate_beta,
    calculate_bondstats_from_trajectory,
    calculate_harmonic_forces,
    eyring_prefactor,
    get_edissoc_from_atomnames,
    harmonic_rates_from_forces,
    harmonic_transition_rate,
    load_trajectory_universe,
    morse_rates_from_forces,
    morse_transition_rate,
    select_bonds,
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

        trajectory_format = self.config.trajectory
        trj = files.input[trajectory_format]
        if trj is None:
            mdp_option = (
                "nstxout-compressed" if trajectory_format == "xtc" else "nstxout"
            )
            m = (
                f"No {trajectory_format} file found. The homolysis plugin computes "
                f"bond distances from the trajectory selected with `trajectory` "
                f"({trajectory_format}). Make sure the preceding MD writes it "
                f"(`{mdp_option}` > 0 in the mdp file)."
            )
            logger.error(m)
            raise ValueError(m)
        logger.info(f"Using trajectory {trj.parent.name}/{trj.name}")

        self.bondstatsfile = trj.with_name(".kimmdy.bondstats")
        md_instance = trj.stem
        timings = self.runmng.timeinfos.get(md_instance)
        if timings is None:
            m = f"No timings from mdp file found for {md_instance}"
            logger.error(m)
            raise ValueError(m)
        top = self.runmng.top

        did_read_bondstats = self.use_cached_bondstats()
        if not did_read_bondstats:
            structure = files.input["gro"]
            if structure is None:
                m = "No gro file found. The homolysis plugin requires a structure file matching the trajectory."
                logger.error(m)
                raise ValueError(m)
            # xtc files may only contain the compressed-x-grps group.
            # NOTE: this lookup duplicates kimmdy-hat (kimmdy_hat/reaction.py)
            # and could be merged into a shared util in the future.
            compressed_group = None
            mdp = self.runmng.mdps.get(md_instance)
            if mdp is not None:
                compressed_group = mdp.get("compressed-x-grps")
            u = load_trajectory_universe(
                structure=structure,
                trajectory=trj,
                compressed_group=compressed_group,
            )
            # bonds are taken from the current topology, so bonds broken
            # in earlier reactions are no longer monitored
            atom_nrs = atom_nrs_from_selection(u, self.config.bond_selection)
            n_selected = len(select_bonds(top=top, atom_nrs=atom_nrs))
            bonds = select_bonds(
                top=top,
                atom_nrs=atom_nrs,
                exclude=self.config.bond_exclude.split(),
            )
            counts = Counter(
                f"{top.atoms[i].atom}-{top.atoms[j].atom}" for i, j in bonds
            )
            logger.info(
                f"bond_selection '{self.config.bond_selection}' matches "
                f"{len(atom_nrs)} atoms with {n_selected} bonds, bond_exclude "
                f"'{self.config.bond_exclude}' excludes {n_selected - len(bonds)}. "
                f"Monitoring {len(bonds)} bonds: "
                + ", ".join(f"{k}: {n}" for k, n in counts.most_common())
            )
            if not bonds:
                logger.warning(
                    "No bonds selected for homolysis, returning no recipes. "
                    "Check bond_selection and bond_exclude."
                )
                self.bondstats = {}
                self.cache_bondstats()
                return RecipeCollection([])
            self.bondstats = calculate_bondstats_from_trajectory(
                top=top,
                u=u,
                bonds=bonds,
                stride=self.config.stride,
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

        temperature, prefactor = self.rate_parameters(logger)
        recipes = []
        for bondkey, stats in self.bondstats.items():
            ai = top.atoms[bondkey[0]]
            aj = top.atoms[bondkey[1]]
            if self.config.check_bound and aj.nr not in ai.bound_to_nrs:
                continue
            typekey = (ai.type, aj.type)
            bondtype = top.ff.bondtypes.get(typekey)
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
                        frequency_factor=prefactor,
                        temperature=temperature,
                    )
                else:
                    ks, fs = harmonic_transition_rate(
                        r_curr=[dist],
                        r_0=b0,
                        k_f=kb,
                        dissociation_energy=edis,
                        frequency_factor=prefactor,
                        temperature=temperature,
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
                        frequency_factor=prefactor,
                        temperature=temperature,
                    )
                else:
                    k = harmonic_rates_from_forces(
                        fs=fs,
                        edis=edis,
                        kb=kb,
                        frequency_factor=prefactor,
                        temperature=temperature,
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

    def rate_parameters(self, logger) -> tuple[float, float]:
        """Temperature in K and prefactor in 1/ps of the Eyring equation.

        A deprecated `arrhenius_equation` is converted to the Eyring equation
        with the transmission coefficient reproducing its frequency factor.
        """
        arrhenius = getattr(self.config, "arrhenius_equation", None)
        if arrhenius is not None:
            temperature = arrhenius.temperature
            prefactor = arrhenius.frequency_factor
            kappa = prefactor / eyring_prefactor(temperature)
            logger.warning(
                "arrhenius_equation is deprecated, use eyring_equation instead. "
                f"Converting it to temperature={temperature:g} K and "
                f"transmission_coefficient={kappa:.4g}."
            )
        else:
            eyring = self.config.eyring_equation
            temperature = eyring.temperature
            prefactor = eyring_prefactor(temperature, eyring.transmission_coefficient)
        logger.info(
            f"Using Eyring equation with temperature={temperature:g} K and "
            f"prefactor={prefactor:.4g} 1/ps."
        )
        return temperature, prefactor

    def use_cached_bondstats(self) -> bool:
        if self.config.recompute_bondstats:
            return False
        if not Path(self.bondstatsfile).exists():
            m = f"bondstatsfile {self.bondstatsfile} does not exist. Not using cached bondstats."
            self.logger.info(m)
            return False

        self.bondstats = bondstats_from_csv(self.bondstatsfile)
        m = (
            f"bondstatsfile {self.bondstatsfile} found. Using cached bondstats of "
            f"{len(self.bondstats)} bonds; bond_selection and bond_exclude are not applied."
        )
        self.logger.info(m)
        return True

    def cache_bondstats(self) -> None:
        bondstats_to_csv(self.bondstats, self.bondstatsfile)
