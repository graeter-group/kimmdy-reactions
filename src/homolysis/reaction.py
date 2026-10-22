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
    calculate_forces,
    calculate_harmonic_forces,
    eyring_prefactor,
    get_edissoc_from_atomnames,
    harmonic_rates_from_forces,
    load_trajectory_universe,
    morse_rates_from_forces,
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
        overwrite = self.config.b0_overwrite != 0.0 or self.config.f0_overwrite != 0.0
        if overwrite:
            logger.warning(
                "Overwriting bond forces with "
                f"b0_overwrite={self.config.b0_overwrite} nm and "
                f"f0_overwrite={self.config.f0_overwrite} kJ/mol/nm "
                "(mainly for testing)."
            )
        recipes = []
        if self.config.check_bound:
            unbound = []
            for bondkey in self.bondstats:
                if bondkey[1] not in top.atoms[bondkey[0]].bound_to_nrs:
                    unbound.append(bondkey)
            for bondkey in unbound:
                del self.bondstats[bondkey]
            logger.info(
                f"check_bound: skipping {len(unbound)} bonds that are no longer bound."
            )
        for bondkey, stats in self.bondstats.items():
            ai = top.atoms[bondkey[0]]
            aj = top.atoms[bondkey[1]]
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
            edis = get_edissoc_from_atomnames(
                atomnames=[ai.atom, aj.atom], edissoc=edissoc, residue=ai.residue
            )
            beta = calculate_beta(kb=kb, edis=edis)

            # averaging forces works here because we typically have
            # one conformational state per calculation
            if self.config.use_morse:
                force = stats["mean_f"]
            elif stats.get("harmonic_f") is not None:
                force = stats["harmonic_f"]
            else:
                # harmonic forces are linear, so the force at the
                # mean distance is the mean force
                force = calculate_harmonic_forces(
                    ds=np.array([stats["mean_d"]]), b0=b0, k=kb
                )[0]
            if overwrite:
                force = self.overwrite_force(force, stats["mean_d"], b0, kb, edis, beta)
            # set negative average forces to 0
            force = max(force, 0)

            fs = np.array([force])
            if self.config.use_morse:
                ks = morse_rates_from_forces(
                    fs=fs,
                    b0=b0,
                    edis=edis,
                    beta=beta,
                    frequency_factor=prefactor,
                    temperature=temperature,
                )
            else:
                ks = harmonic_rates_from_forces(
                    fs=fs,
                    edis=edis,
                    kb=kb,
                    frequency_factor=prefactor,
                    temperature=temperature,
                )

            recipes.append(
                Recipe(
                    recipe_steps=[
                        Break(atom_id_1=bondkey[0], atom_id_2=bondkey[1]),
                        Relax(),
                    ],
                    rates=list(ks),
                    timespans=[(t0, t1)],
                )
            )

        return RecipeCollection(recipes)

    def overwrite_force(
        self,
        force: float,
        mean_d: float,
        b0: float,
        kb: float,
        edis: float,
        beta: float,
    ) -> float:
        """Apply b0_overwrite and f0_overwrite (mainly for testing) to a bond force."""
        if self.config.b0_overwrite != 0.0:
            # bondstats forces belong to the force field b0, so with an
            # offset the force is evaluated at the mean distance
            ds = np.array([mean_d])
            b0 = b0 + self.config.b0_overwrite
            if self.config.use_morse:
                force = calculate_forces(ds=ds, b0=b0, edis=edis, beta=beta)[0]
            else:
                force = calculate_harmonic_forces(ds=ds, b0=b0, k=kb)[0]
        return force - self.config.f0_overwrite

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
