"""Extract the first-chain raw topology from the provided 10-chain CLS ITP.

This script intentionally keeps the extracted bonded interactions unchanged.
Task 2 first needs a raw single-chain ITP that preserves the original error, so
the chemical correction is performed by a separate script after visualization
and valence-rule checks have located the bad bond.
"""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import scripts.utils.gmx_io as gmxio

INPUT_ITP = PROJECT_ROOT / "topology" / "CLS_10chains.itp"
INPUT_GRO = PROJECT_ROOT / "single_chain" / "CLS_single.gro"
SINGLE_CHAIN_DIR = PROJECT_ROOT / "single_chain"
RAW_ITP = SINGLE_CHAIN_DIR / "CLS_single_raw.itp"
ATOMS_PER_CHAIN = 212
MOLECULE_NAME = "CLS"


def main() -> None:
    """Create the raw single-chain ITP file from the first 212 source atoms."""
    raw_chain = extract_first_chain()
    gmxio.save_itp(RAW_ITP, raw_chain, molecule_name=MOLECULE_NAME)

    print(f"Wrote raw topology: {RAW_ITP}")
    print(f"Raw atoms: {len(raw_chain.atoms)}")
    print(f"Raw bonds: {len(raw_chain.bonds)}")
    print(f"Pairs: {len(raw_chain.adjusts)}")
    print(f"Angles: {len(raw_chain.angles)}")
    print(f"Dihedrals: {len(raw_chain.dihedrals)}")

    check_gro_matches_itp_first_atoms(INPUT_GRO, INPUT_ITP, ATOMS_PER_CHAIN)


def extract_first_chain():
    """Load the source ITP and return the first-chain ParmEd topology."""
    topology = gmxio.load_itp(INPUT_ITP, parametrize=False)
    return gmxio.first_atoms(topology, ATOMS_PER_CHAIN)


def check_gro_matches_itp_first_atoms(gro_path: Path, itp_path: Path, atom_count: int) -> None:
    """Validate that the single-chain GRO atoms match the source ITP atom rows."""
    gro = gmxio.load_gro(gro_path)
    itp = gmxio.load_itp(itp_path, parametrize=False)
    gro_atoms = [atom_record(atom) for atom in gro.atoms]
    itp_atoms = [atom_record(atom) for atom in itp.atoms[:atom_count]]

    if len(gro_atoms) != atom_count:
        raise AssertionError(
            f"{gro_path} contains {len(gro_atoms)} atoms, expected {atom_count}."
        )
    if len(itp_atoms) != atom_count:
        raise AssertionError(
            f"{itp_path} has only {len(itp_atoms)} atom rows in the first-chain range."
        )

    mismatches: list[str] = []
    for index, (gro_atom, itp_atom) in enumerate(zip(gro_atoms, itp_atoms), start=1):
        compared_fields = ("atomnr", "resid", "resname", "atomname")
        for field in compared_fields:
            if gro_atom[field] != itp_atom[field]:
                mismatches.append(
                    f"atom {index}: {field} differs "
                    f"(gro={gro_atom[field]!r}, itp={itp_atom[field]!r})"
                )

    if mismatches:
        details = "\n".join(mismatches[:10])
        raise AssertionError(
            "single_chain/CLS_single.gro does not match the first "
            f"{atom_count} [ atoms ] rows of topology/CLS_10chains.itp:\n{details}"
        )

    print("GRO/ITP atom correspondence check: PASS")
    print(f"Checked atoms: {atom_count}")
    print(
        "First atom: "
        f"{gro_atoms[0]['atomnr']} {gro_atoms[0]['resid']} "
        f"{gro_atoms[0]['resname']} {gro_atoms[0]['atomname']}"
    )
    print(
        "Last atom: "
        f"{gro_atoms[-1]['atomnr']} {gro_atoms[-1]['resid']} "
        f"{gro_atoms[-1]['resname']} {gro_atoms[-1]['atomname']}"
        )


def atom_record(atom) -> dict[str, int | str]:
    """Return the atom identity fields used to compare GRO and ITP order."""
    return {
        "atomnr": atom.idx + 1,
        "resid": atom.residue.idx + 1,
        "resname": atom.residue.name,
        "atomname": atom.name,
    }


if __name__ == "__main__":
    main()
