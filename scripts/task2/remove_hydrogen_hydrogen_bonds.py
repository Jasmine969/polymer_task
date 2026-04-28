"""Create the corrected single-chain ITP by removing direct H-H bonds.

The task-2 workflow deliberately separates extraction from correction:

1. ``extract_first_chain_itp.py`` writes ``CLS_single_raw.itp`` without changing
   the original bonded topology.
2. Visualization and valence checks identify the chemically invalid direct
   H-H bond in that raw file.
3. This script removes only direct H-H entries from ``[ bonds ]`` and writes
   ``CLS_single_corrected.itp`` for the final validation.

The ``[ pairs ]`` section is not modified. In this system the atom pair 1-9 is
still a valid 1-4 interaction after the direct H-H bond is removed.
"""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import scripts.utils.gmx_io as gmxio

SINGLE_CHAIN_DIR = PROJECT_ROOT / "single_chain"
RAW_ITP = SINGLE_CHAIN_DIR / "CLS_single_raw.itp"
CORRECTED_ITP = SINGLE_CHAIN_DIR / "CLS_single_corrected.itp"
MOLECULE_NAME = "CLS"


def main() -> None:
    """Read the raw single-chain ITP and write the corrected ITP files."""
    topology = gmxio.load_itp(RAW_ITP, parametrize=False)
    removed_hh_bonds = remove_hydrogen_hydrogen_bonds(topology)

    gmxio.save_itp(CORRECTED_ITP, topology, molecule_name=MOLECULE_NAME)

    print(f"Read raw topology: {RAW_ITP}")
    print(f"Wrote corrected topology: {CORRECTED_ITP}")
    print(f"Atoms: {len(topology.atoms)}")
    print(f"Bonds: {len(topology.bonds)}")
    print(f"Removed H-H bonds: {removed_hh_bonds}")
    print(f"Pairs: {len(topology.adjusts)}")
    print(f"Angles: {len(topology.angles)}")
    print(f"Dihedrals: {len(topology.dihedrals)}")


def remove_hydrogen_hydrogen_bonds(topology) -> int:
    """Remove direct bonds whose two endpoint atom names both indicate hydrogen.

    Parameters
    ----------
    topology:
        ParmEd ``GromacsTopologyFile`` loaded from ``CLS_single_raw.itp``.

    Returns
    -------
    int
        Number of removed direct H-H bonds. The caller prints this count so the
        task report can document exactly what changed.
    """
    bad_bonds = [
        bond
        for bond in topology.bonds
        if is_hydrogen_atom(bond.atom1.name) and is_hydrogen_atom(bond.atom2.name)
    ]

    for bond in bad_bonds:
        topology.bonds.remove(bond)

    return len(bad_bonds)


def is_hydrogen_atom(atom_name: str) -> bool:
    """Return whether a GROMACS atom name looks like a hydrogen atom.

    The provided ITP uses names such as ``H1`` and ``H102``. Checking the first
    letter is sufficient for this project and avoids relying on atom types,
    which are force-field labels rather than element symbols.
    """
    return atom_name.upper().startswith("H")


if __name__ == "__main__":
    main()
