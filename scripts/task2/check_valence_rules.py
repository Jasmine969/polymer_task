"""Check simple chemistry valence rules in a GROMACS ITP bond graph.

This is a topology sanity check, not a full bond-order perception tool.
GROMACS ITP files define bonded interactions but do not store explicit single,
double, or aromatic bond orders. The checks here therefore use coordination
counts only: H should have exactly 1 bond, O should have 1 to 2 bonds, and C
should have 2 to 4 bonds.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
# ITP_PATH = ROOT / "single_chain" / "CLS_single_raw.itp"
ITP_PATH = ROOT / "single_chain" / "CLS_single_corrected.itp"

MAX_BONDS = {
    "C": 4,
    "O": 2,
}

MIN_BONDS = {
    "C": 2,
    "O": 1,
}


def element_from_atom_name(atom_name: str) -> str:
    """Infer the element symbol from a GROMACS atom name such as H102 or C1."""
    letters = "".join(char for char in atom_name if char.isalpha())
    if not letters:
        return ""
    return letters[0].upper()


def read_atoms_and_bonds(path: Path) -> tuple[dict[int, dict[str, str]], list[tuple[int, int]]]:
    """Read atom names and bond pairs from an ITP file.

    Parameters
    ----------
    path:
        GROMACS ITP file whose ``[ atoms ]`` and ``[ bonds ]`` sections should
        be checked.

    Returns
    -------
    tuple[dict[int, dict[str, str]], list[tuple[int, int]]]
        A mapping from atom number to parsed atom metadata, plus a list of
        bonded atom-number pairs. Only the fields needed for simple valence
        checks are extracted.
    """
    section: str | None = None
    atoms: dict[int, dict[str, str]] = {}
    bonds: list[tuple[int, int]] = []

    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("["):
            section = stripped.strip("[]").strip().lower()
            continue

        data = line.split(";", 1)[0].strip()
        if not data:
            continue

        fields = data.split()
        if section == "atoms" and len(fields) >= 5 and fields[0].isdigit():
            atom_id = int(fields[0])
            atom_name = fields[4]
            atoms[atom_id] = {
                "type": fields[1],
                "resnr": fields[2],
                "resname": fields[3],
                "name": atom_name,
                "element": element_from_atom_name(atom_name),
            }
        elif section == "bonds" and len(fields) >= 2 and fields[0].isdigit() and fields[1].isdigit():
            bonds.append((int(fields[0]), int(fields[1])))

    return atoms, bonds


def main() -> None:
    check_itp(ITP_PATH)



def check_itp(itp_path: Path) -> None:
    """Report atoms whose coordination violates simple element rules.

    Parameters
    ----------
    itp_path:
        ITP file to inspect. The output includes the number of atoms, the number
        of bonds, and any atoms that break the H/O/C coordination rules.
        Carbon and oxygen are checked against both lower and upper coordination
        limits because dangling C/O atoms are chemically suspicious even when
        they do not exceed their maximum valence.
    """
    atoms, bonds = read_atoms_and_bonds(itp_path)
    adjacency = {atom_id: [] for atom_id in atoms}

    for atom1, atom2 in bonds:
        adjacency[atom1].append(atom2)
        adjacency[atom2].append(atom1)

    violations: list[str] = []
    for atom_id, atom in atoms.items():
        element = atom["element"]
        bond_count = len(adjacency[atom_id])
        neighbors = ", ".join(
            f"{neighbor}:{atoms[neighbor]['name']}"
            for neighbor in sorted(adjacency[atom_id])
        )

        if element == "H" and bond_count != 1:
            violations.append(
                f"{atom_id:4d} {atom['name']:>4s} H has {bond_count} bonds -> {neighbors}"
            )
        elif element in MIN_BONDS and bond_count < MIN_BONDS[element]:
            violations.append(
                f"{atom_id:4d} {atom['name']:>4s} {element} has {bond_count} bonds "
                f"(min {MIN_BONDS[element]}) -> {neighbors}"
            )
        elif element in MAX_BONDS and bond_count > MAX_BONDS[element]:
            violations.append(
                f"{atom_id:4d} {atom['name']:>4s} {element} has {bond_count} bonds "
                f"(max {MAX_BONDS[element]}) -> {neighbors}"
            )

    print(f"Input: {itp_path}")
    print(f"Atoms: {len(atoms)}")
    print(f"Bonds: {len(bonds)}")
    print(f"Violations: {len(violations)}")

    for violation in violations:
        print(violation)


if __name__ == "__main__":
    main()
