"""
Tests for default_reactions reaction plugin.

Is skipped if KIMMDY was installed without the plugin.
"""

import logging

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
    Plumed_dict,
    read_top,
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

# backbone bonds of npt.gro as in its former plumed.dat
# (backbone minus "H* O*")
BACKBONE_BONDS = {
    tuple(pair.split("-"))
    for pair in """
    5-7 7-9 9-20 20-22 22-24 24-44 44-46 46-48 48-51 51-53 53-63 63-65 65-67 67-78
    78-80 80-82 82-84 84-87 87-89 89-99 99-101 101-103 103-105 105-118 118-120 120-122
    122-125 125-127 127-137 137-139 139-141 141-152 152-154 154-156 156-158 158-161
    161-163 163-165 165-180 180-182 182-184 184-190 190-192 192-194 194-197 197-199
    199-209 209-211 211-213 213-224 224-226 226-228 228-230 230-233 233-235 235-237
    237-248 248-250 250-252 252-259 259-261 261-263 263-266 266-268 268-270 270-290
    290-292 292-294 294-305 305-307 307-309 309-312 312-314 314-316 316-323 323-325
    325-336 336-338 338-340 340-342 342-345 345-347 347-349 349-355 355-357 357-359
    359-370 370-372 372-374 374-377 377-379 379-381 381-388 388-390 390-401 401-403
    403-405 405-407 407-410 410-412 412-414 414-434 434-436 436-438 438-446 446-448
    448-450 450-453 453-455 455-457 457-463 463-465 465-476 476-478 478-480 480-482
    482-485 485-487 487-489 489-495 495-497 497-499 499-517 517-519 519-521 521-524
    524-526 526-528 528-536 536-538 538-540 540-560 560-562 562-564 564-567 567-569
    569-571 571-582 582-584 594-596 596-598 598-604 604-606 606-608 608-628 628-630
    630-632 632-635 635-637 637-647 647-649 649-651 651-662 662-664 664-666 666-668
    668-671 671-673 673-675 675-681 681-683 683-685 685-697 697-699 699-701 701-704
    704-706 706-708 708-715 715-717 717-728 728-730 730-732 732-734 734-737 737-739
    739-741 741-753 753-755 755-757 757-767 767-769 769-771 771-774 774-776 776-778
    778-784 784-786 786-797 797-799 799-801 801-803 803-806 806-808 808-810 810-821
    821-823 823-825 825-831 831-833 833-835 835-838 838-840 840-842 842-862 862-864
    864-866 866-874 874-876 876-878 878-881 881-883 883-885 885-895 895-897 897-908
    908-910 910-912 912-914 914-917 917-919 919-921 921-928 928-930 930-932 932-940
    940-942 942-944 944-947 947-949 949-959 959-961 961-963 963-974 974-976 976-978
    978-980 980-983 983-985 985-987 987-1007 1007-1009 1009-1011 1011-1019 1019-1021
    1021-1023 1023-1026 1026-1028 1028-1030 1030-1043 1043-1045 1045-1056 1056-1058
    1058-1060 1060-1062 1062-1065 1065-1067 1067-1069 1069-1082 1082-1084 1084-1086
    1086-1104 1104-1106 1106-1108 1108-1111 1111-1113 1113-1115 1115-1126 1126-1128
    1128-1130 1130-1150 1150-1152 1162-1164 1164-1166 1166-1173 1173-1175 1175-1177
    1177-1180 1180-1182 1182-1184 1184-1195 1195-1197 1197-1199 1199-1219 1219-1221
    1221-1223 1223-1226 1226-1228 1228-1238 1238-1240 1240-1242 1242-1253 1253-1255
    1255-1257 1257-1259 1259-1262 1262-1264 1264-1274 1274-1276 1276-1278 1278-1280
    1280-1293 1293-1295 1295-1297 1297-1300 1300-1302 1302-1312 1312-1314 1314-1316
    1316-1327 1327-1329 1329-1331 1331-1333 1333-1336 1336-1338 1338-1340 1340-1355
    1355-1357 1357-1359 1359-1365 1365-1367 1367-1369 1369-1372 1372-1374 1374-1384
    1384-1386 1386-1388 1388-1399 1399-1401 1401-1403 1403-1405 1405-1408 1408-1410
    1410-1412 1412-1423 1423-1425 1425-1427 1427-1434 1434-1436 1436-1438 1438-1441
    1441-1443 1443-1445 1445-1465 1465-1467 1467-1469 1469-1480 1480-1482 1482-1484
    1484-1487 1487-1489 1489-1491 1491-1498 1498-1500 1500-1511 1511-1513 1513-1515
    1515-1517 1517-1520 1520-1522 1522-1524 1524-1530 1530-1532 1532-1534 1534-1545
    1545-1547 1547-1549 1549-1552 1552-1554 1554-1556 1556-1563 1563-1565 1565-1576
    1576-1578 1578-1580 1580-1582 1582-1585 1585-1587 1587-1589 1589-1609 1609-1611
    1611-1613 1613-1621 1621-1623 1623-1625 1625-1628 1628-1630 1630-1632 1632-1638
    1638-1640 1640-1651 1651-1653 1653-1655 1655-1657 1657-1660 1660-1662 1662-1664
    1664-1670 1670-1672 1672-1674 1674-1692 1692-1694 1694-1696 1696-1699 1699-1701
    1701-1703 1703-1711 1711-1713
    """.split()
}


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
    ffbonded = read_top(Path("ffbonded.itp"))
    edissoc = read_edissoc(Path("edissoc.dat"))

    files = {
        "top": top,
        "ffbonded": ffbonded,
        "edissoc": edissoc,
    }

    return files


@pytest.fixture
def deprecation_log(caplog, monkeypatch):
    """Capture deprecation warnings logged by kimmdy.plugin_utils."""
    # the kimmdy logging config disables propagation, which hides records from caplog
    monkeypatch.setattr(logging.getLogger("kimmdy"), "propagate", True)
    caplog.set_level(logging.WARNING, logger="kimmdy.plugin_utils")
    return caplog


## test homolysis
def test_get_atomnrs(deprecation_log):
    plumed = Plumed_dict(
        other=[],
        labeled_action={"d1": {"keyword": "DISTANCE", "atoms": ["9", "7"]}},
        prints=[],
    )
    atomnrs = get_atomnrs_from_plumedid("d1", plumed)
    assert atomnrs == ["7", "9"]
    assert "Deprecated function get_atomnrs_from_plumedid" in deprecation_log.text


def test_get_atomtypes(homolysis_files, deprecation_log):
    atomnrs = ["7", "9"]
    atomtypes, atomnames = get_atominfo_from_atomnrs(atomnrs, homolysis_files["top"])
    assert atomtypes == ["N", "CT"]
    assert atomnames == ["N", "CA"]
    assert "Deprecated function get_atominfo_from_atomnrs" in deprecation_log.text


def test_lookup_bondprm(homolysis_files, deprecation_log):
    b0, kb = get_bondprm_from_atomtypes(["CT", "C"], homolysis_files["ffbonded"])
    assert abs(b0 - 0.15220) < 1e-9
    assert abs(kb - 265265.6) < 1e-9
    assert "Deprecated function get_bondprm_from_atomtypes" in deprecation_log.text


def test_lookup_edissoc(homolysis_files):
    e_dis = get_edissoc_from_atomnames(
        ["CA", "C"], homolysis_files["edissoc"], residue="GLY"
    )
    assert abs(e_dis - 341) < 1e-9


def test_fail_lookup_bondprm(homolysis_files, deprecation_log):
    with pytest.raises(KeyError, match="Did not find bond parameters for atomtypes"):
        b0, kb = get_bondprm_from_atomtypes(
            ["X", "Z"],
            homolysis_files["ffbonded"],
        )
    assert "Deprecated function get_bondprm_from_atomtypes" in deprecation_log.text


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
    bonds = _homolysis_from_trajectory(homolysis_files["top"], bond_exclude="H* O*")
    assert bonds == BACKBONE_BONDS


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
