"""Small ParmEd-based helpers for GROMACS coordinate and topology files."""

from pathlib import Path

import parmed as pmd
from parmed import Structure
from parmed.gromacs import GromacsTopologyFile


def load_gro(path: Path) -> Structure:
    """Load a GROMACS .gro file with ParmEd."""
    return pmd.load_file(str(path))


def save_gro(path: Path, structure: Structure) -> None:
    """Save a ParmEd structure as a GROMACS .gro file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    structure.save(str(path), overwrite=True)


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


def first_atoms(structure: Structure, atom_count: int):
    """Return a ParmEd structure containing the first ``atom_count`` atoms."""
    if len(structure.atoms) < atom_count:
        raise ValueError(
            f"Structure contains {len(structure.atoms)} atoms, "
            f"fewer than requested {atom_count} atoms."
        )
    return structure[:atom_count]


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
