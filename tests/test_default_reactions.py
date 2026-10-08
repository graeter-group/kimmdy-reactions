"""
Tests for default_reactions reaction plugin.

Is skipped if KIMMDY was installed without the plugin.
"""

import pytest
from kimmdy.runmanager import RunManager

pytest.importorskip("homolysis")
pytest.importorskip("hat_naive")
pytest.importorskip("dummyreaction")

from dataclasses import dataclass
from typing import Callable
import os
import shutil
from pathlib import Path
import numpy as np
from homolysis.reaction import Homolysis
from kimmdy.plugins import discover_plugins
from kimmdy.config import Config
from kimmdy.recipe import Break
from kimmdy.topology.topology import Topology
from kimmdy.parsing import (
    read_plumed,
    read_top,
    read_distances_dat,
    read_edissoc,
)
from kimmdy.plugin_utils import (
    get_atomnrs_from_plumedid,
    get_atominfo_from_atomnrs,
    get_bondprm_from_atomtypes,
    get_edissoc_from_atomnames,
    morse_transition_rate,
)
from kimmdy.tasks import TaskFiles

discover_plugins()


@dataclass
class DummyRunmanager(RunManager):
    top: Topology
    config: Config


@dataclass
class DummyFiles(TaskFiles):
    get_latest: Callable = lambda p: f"DummyCallable"


@pytest.fixture
def homolysis_files(tmp_path: Path):
    file_dir = Path(__file__).parent / "test_default_reactions"
    shutil.copytree(file_dir, tmp_path, dirs_exist_ok=True)
    os.chdir(tmp_path)

    assetsdir = Path(__file__).parent / "assets"
    Path(tmp_path / "amber99sb-star-ildnp.ff").symlink_to(
        assetsdir / "amber99sb-star-ildnp.ff",
        target_is_directory=True,
    )

    top = Topology(read_top(Path("topol.top")))
    plumed = read_plumed(Path("plumed.dat"))
    distances = read_distances_dat(Path("distances.dat"))
    distances_avg = read_distances_dat(Path("distances_avg.dat"))
    ffbonded = read_top(Path("ffbonded.itp"))
    edissoc = read_edissoc(Path("edissoc.dat"))

    files = {
        "top": top,
        "plumed": plumed,
        "distances": distances,
        "distances_avg": distances_avg,
        "ffbonded": ffbonded,
        "edissoc": edissoc,
    }

    return files


## test homolysis
def test_get_atomnrs(homolysis_files):
    atomnrs = get_atomnrs_from_plumedid("d1", homolysis_files["plumed"])
    assert atomnrs == ["7", "9"]


def test_get_atomtypes(homolysis_files):
    atomnrs = ["7", "9"]
    atomtypes, atomnames = get_atominfo_from_atomnrs(atomnrs, homolysis_files["top"])
    assert atomtypes == ["N", "CT"]
    assert atomnames == ["N", "CA"]


def test_lookup_bondprm(homolysis_files):
    b0, kb = get_bondprm_from_atomtypes(["CT", "C"], homolysis_files["ffbonded"])
    assert abs(b0 - 0.15220) < 1e-9
    assert abs(kb - 265265.6) < 1e-9


def test_lookup_edissoc(homolysis_files):
    e_dis = get_edissoc_from_atomnames(
        ["CA", "C"], homolysis_files["edissoc"], residue="GLY"
    )
    assert abs(e_dis - 341) < 1e-9


def test_fail_lookup_bondprm(homolysis_files):
    with pytest.raises(KeyError, match="Did not find bond parameters for atomtypes"):
        b0, kb = get_bondprm_from_atomtypes(
            ["X", "Z"],
            homolysis_files["ffbonded"],
        )


def test_morse_transition_rate(homolysis_files):
    b0, kb = get_bondprm_from_atomtypes(["CT", "C"], homolysis_files["ffbonded"])
    e_dis = get_edissoc_from_atomnames(
        ["CA", "C"], homolysis_files["edissoc"], "general"
    )

    rs_ref = list(np.linspace(0.9, 1.3, 8) * b0)
    ks, fs = morse_transition_rate(rs_ref, b0, e_dis, kb)

    ks_ref = np.asarray(
        [
            0.00000000e00,
            0.00000000e00,
            1.01587680e-40,
            2.89436858e-11,
            3.69638684e-03,
            2.17458162e-01,
            2.88000000e-01,
            2.88000000e-01,
        ]
    )
    fs_ref = np.asarray(
        [
            -6357.20305659,
            -2100.01133653,
            540.87430205,
            2094.72356994,
            2927.66683909,
            3291.55831528,
            3362.580289,
            3362.580289,
        ]
    )
    assert all(np.isclose(ks, ks_ref))
    assert all(np.isclose(fs, fs_ref))


def _homolysis_from_trajectory(top: Topology, **options):
    """Run the homolysis plugin on a protein-only xtc of npt.gro."""
    import logging
    from types import SimpleNamespace

    import MDAnalysis as mda
    from kimmdy.plugin_utils import SOL_RESNAMES
    from kimmdy.runmanager import TimeInfo

    structure = mda.Universe("npt.gro")
    protein = structure.select_atoms(f"not resname {' '.join(SOL_RESNAMES)}")
    with mda.Writer("prod.xtc", protein.n_atoms) as w:
        for _ in range(2):
            w.write(protein)

    config = SimpleNamespace(
        trajectory="xtc",
        bond_selection="backbone",
        bond_exclude="C-N H* O*",
        dt_distances=0.0,
        recompute_bondstats=True,
        check_bound=True,
        use_morse=True,
        b0_overwrite=0.0,
        f0_overwrite=0.0,
        arrhenius_equation=SimpleNamespace(frequency_factor=0.288, temperature=300),
    )
    for k, v in options.items():
        setattr(config, k, v)
    runmng = SimpleNamespace(
        config=SimpleNamespace(reactions=SimpleNamespace(homolysis=config)),
        top=top,
        timeinfos={
            "prod": TimeInfo(nsteps=1000, dt=0.002, trr_nst=100, xtc_nst=100, t_max=2.0)
        },
        mdps={"prod": {"compressed-x-grps": "Protein"}},
    )
    files = SimpleNamespace(
        logger=logging.getLogger("test_homolysis"),
        input={
            "xtc": Path("prod.xtc"),
            "trr": None,
            "gro": Path("npt.gro"),
            "edissoc.dat": Path("edissoc.dat"),
        },
    )
    rc = Homolysis("homolysis", runmng).get_recipe_collection(files)
    bonds = set()
    for recipe in rc.recipes:
        step = recipe.recipe_steps[0]
        assert isinstance(step, Break)
        bonds.add(tuple(sorted((step.atom_id_1, step.atom_id_2), key=int)))
    return bonds


def test_homolysis_bonds_from_trajectory_match_plumed(homolysis_files):
    from kimmdy.plugin_utils import read_plumed_input

    plumed_bonds = set(read_plumed_input(Path("plumed.dat")).keys())
    bonds = _homolysis_from_trajectory(homolysis_files["top"], bond_exclude="H* O*")
    assert bonds == plumed_bonds


def test_homolysis_does_not_monitor_broken_bonds(homolysis_files):
    top = homolysis_files["top"]
    bonds = _homolysis_from_trajectory(top)
    broken = sorted(bonds, key=lambda b: int(b[0]))[0]
    top.break_bond(broken)
    bonds_after = _homolysis_from_trajectory(top)
    assert broken not in bonds_after
    assert bonds_after == bonds - {broken}


def test_homolysis_requires_selected_trajectory(homolysis_files):
    with pytest.raises(ValueError, match="trajectory: xtc"):
        _homolysis_from_trajectory(homolysis_files["top"], trajectory="trr")
