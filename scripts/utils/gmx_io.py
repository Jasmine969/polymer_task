"""Small ParmEd-based helpers for GROMACS coordinate and topology files."""

from pathlib import Path

import parmed as pmd
from parmed import Structure
from parmed.gromacs import GromacsTopologyFile


# ==================== general tools =============================
def first_atoms(structure: Structure, atom_count: int):
    """Return a ParmEd structure containing the first ``atom_count`` atoms."""
    if len(structure.atoms) < atom_count:
        raise ValueError(
            f"Structure contains {len(structure.atoms)} atoms, "
            f"fewer than requested {atom_count} atoms."
        )
    return structure[:atom_count]


# ========================= gro file ================================
def load_gro(path: Path) -> Structure:
    """Load a GROMACS .gro file with ParmEd."""
    return pmd.load_file(str(path))


def save_gro(path: Path, structure: Structure) -> None:
    """Save a ParmEd structure as a GROMACS .gro file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    structure.save(str(path), overwrite=True)


# ============================== itp file ========================
def load_itp(path: Path, parametrize: bool = False) -> GromacsTopologyFile:
    """Load a GROMACS .itp/.top file with ParmEd.

    ``parametrize=False`` is the safer default for this project because the
    provided polymer ITP comments out the atomtypes section.
    """
    return GromacsTopologyFile(str(path), parametrize=parametrize)


def save_itp(path: Path, topology: GromacsTopologyFile, molecule_name: str | None = None) -> None:
    """Save a ParmEd GromacsTopologyFile as a standalone ITP file.

    ParmEd names multi-residue molecules ``system1`` when writing an ITP. When a
    GROMACS top file needs a specific molecule name, ``molecule_name`` rewrites
    the first moleculetype data line after ParmEd writes the file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    topology.write(str(path), itp=True)

    if molecule_name is not None:
        _replace_first_moleculetype_name(path, molecule_name)


def read_defaults(path: Path) -> tuple[str, str]:
    """Return the comment and data line from a GROMACS ``[ defaults ]`` block.

    The returned tuple is ``(comment_line, data_line)``.  The comment line may
    be empty if the source topology does not document the columns.  A missing
    data line is treated as an error because a combined topology header cannot
    be generated without non-bonded defaults.
    """

    block = read_directive_block(path, "defaults")
    comments = [line for line in block if line.strip().startswith(";")]
    data_lines = [
        line.strip()
        for line in block
        if line.strip() and not line.strip().startswith(";")
    ]
    if not data_lines:
        raise ValueError(f"{path} does not contain a [ defaults ] data line.")
    return (comments[-1].strip() if comments else "", data_lines[0])


def read_atomtype_rows(path: Path) -> list[str]:
    """Return data rows from a GROMACS ``[ atomtypes ]`` block.

    Inline comments are removed and blank/comment-only lines are skipped.  The
    raw row text is returned so callers can preserve source values without
    converting floating-point formats.
    """

    rows: list[str] = []
    for line in read_directive_block(path, "atomtypes"):
        row = strip_inline_comment(line).strip()
        if row and not row.startswith(";"):
            rows.append(row)
    return rows


def merge_atomtypes(paths: list[Path]) -> list[str]:
    """Merge ``[ atomtypes ]`` rows from topology files by atom type name.

    Later duplicate atom type names are skipped after checking that their
    parameters match the first occurrence.  This intentionally uses an ordered
    dictionary-style pass rather than a table merge so source ordering is stable,
    six-column water atomtypes can be normalized to the seven-column project
    style, and parameter conflicts raise clear errors.
    """

    rows_by_name: dict[str, str] = {}
    keys_by_name: dict[str, tuple[str, ...]] = {}

    for path in paths:
        for row in read_atomtype_rows(path):
            fields = row.split()
            atomtype_name = fields[0]
            comparison_key = atomtype_comparison_key(fields)
            if atomtype_name in rows_by_name:
                if keys_by_name[atomtype_name] != comparison_key:
                    raise ValueError(
                        f"Conflicting atomtype {atomtype_name!r} in {path}: "
                        f"{row!r} differs from {rows_by_name[atomtype_name]!r}"
                    )
                continue
            rows_by_name[atomtype_name] = format_atomtype_row(fields)
            keys_by_name[atomtype_name] = comparison_key

    return list(rows_by_name.values())


def write_forcefield_header(
        path: Path,
        defaults_comment: str,
        defaults_line: str,
        atomtype_rows: list[str],
        header_comments: list[str] | None = None,
) -> None:
    """Write a small GROMACS force-field header with defaults and atomtypes.

    This is useful when molecule ITP files must be included after all atom type
    definitions.  ``atomtype_rows`` should already be deduplicated and formatted,
    for example by ``merge_atomtypes``.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    lines = list(header_comments or [])
    if lines:
        lines.append("")
    lines.extend(
        [
            "[ defaults ]",
            defaults_comment or "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ",
            defaults_line,
            "",
            "[ atomtypes ]",
            "; name  bond_type     mass     charge   ptype   sigma         epsilon",
        ]
    )
    lines.extend(atomtype_rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_directive_block(path: Path, directive: str) -> list[str]:
    """Extract lines inside one top-level GROMACS directive block.

    The block ends when another ``[ directive ]`` begins or when a preprocessor
    include/define line starts.  The latter matters for compact topology files
    where ``[ atomtypes ]`` is followed directly by ``#include``.
    """

    lines = path.read_text(encoding="utf-8").splitlines()
    block: list[str] = []
    in_block = False

    for line in lines:
        stripped = line.strip()
        if in_block and stripped.startswith("#"):
            break
        if stripped.startswith("[") and stripped.endswith("]"):
            current = stripped.strip("[]").strip().lower()
            if in_block and current != directive.lower():
                break
            in_block = current == directive.lower()
            continue
        if in_block:
            block.append(line)

    return block


def atomtype_comparison_key(fields: list[str]) -> tuple[str, ...]:
    """Normalize atomtype fields for duplicate-parameter comparisons."""

    if len(fields) == 6:
        # Some files omit the explicit bond_type column for water atomtypes.
        name, mass, charge, ptype, sigma, epsilon = fields
        return (name, name, mass, charge, ptype, sigma, epsilon)
    if len(fields) >= 7:
        return tuple(fields[:7])
    raise ValueError(f"Unexpected atomtype row fields: {fields}")


def format_atomtype_row(fields: list[str]) -> str:
    """Format one atomtype row in the seven-column style used by this project."""

    if len(fields) == 6:
        name, mass, charge, ptype, sigma, epsilon = fields
        bond_type = name
    elif len(fields) >= 7:
        name, bond_type, mass, charge, ptype, sigma, epsilon = fields[:7]
    else:
        raise ValueError(f"Unexpected atomtype row fields: {fields}")

    return (
        f"{name:<7} {bond_type:<7} {mass:>10} {charge:>8}   "
        f"{ptype:<1}      {sigma:>11}   {epsilon:>11}"
    )


def strip_inline_comment(line: str) -> str:
    """Remove the inline comment portion of one topology line."""

    return line.split(";", maxsplit=1)[0]


def _replace_first_moleculetype_name(path: Path, molecule_name: str) -> None:
    """Rewrite the first [ moleculetype ] data row in a ParmEd-written ITP."""
    lines = path.read_text(encoding="utf-8").splitlines()
    in_moleculetype = False

    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.lower() == "[ moleculetype ]":
            in_moleculetype = True
            continue

        if not in_moleculetype or not stripped or stripped.startswith(";"):
            continue

        fields = stripped.split()
        nrexcl = fields[1] if len(fields) > 1 else "3"
        lines[index] = f"{molecule_name:<16}{nrexcl}"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return

    raise ValueError(f"{path} does not contain a [ moleculetype ] data row.")


if __name__ == '__main__':
    struc = pmd.load_file('../../systems/water/initial.pdb')
    print(struc)
    save_gro(Path('F:\\polymer_task\\scripts\\utils\\test.gro'), struc)
