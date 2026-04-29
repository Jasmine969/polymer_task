"""Helpers for exporting PDB files from project structures.

The functions here cover two project use cases:

* PackMol templates, where coordinates are enough and no ``CONECT`` records are
  needed.
* Visual-inspection PDB files, where a GROMACS topology is available and
  topology-derived ``CONECT`` records should be written explicitly.
"""

from pathlib import Path

from parmed import Structure
from parmed.gromacs import GromacsTopologyFile

import scripts.utils.gmx_io as gmxio


def load_pdb_source(gro_path: Path, top_path: Path | None = None) -> Structure:
    """Load the structure that should be written to PDB.

    Parameters
    ----------
    gro_path:
        Coordinate file used for atom positions.
    top_path:
        Optional GROMACS topology.  When provided, ParmEd combines topology and
        coordinates so the resulting structure contains bonded connectivity.
        When omitted, only the coordinate file is loaded, which is suitable for
        PackMol molecule templates where bonds are not needed.
    """

    if top_path is None:
        return gmxio.load_gro(gro_path)
    return GromacsTopologyFile(str(top_path), xyz=str(gro_path))


def save_pdb(path: Path, structure: Structure, include_conect: bool = False) -> None:
    """Save a ParmEd structure as PDB, optionally adding CONECT records.

    ParmEd writes the coordinate and residue records.  If ``include_conect`` is
    true, this helper rebuilds ``CONECT`` records from ``structure.bonds`` after
    writing the file, which is useful for visualization tools that need explicit
    bond records in the PDB file.
    """

    path.parent.mkdir(parents=True, exist_ok=True)
    # Some GRO-only ParmEd structures do not define ``symmetry`` until a writer
    # asks for it.  Defining it explicitly keeps PDB export stable after in-memory
    # coordinate edits such as centering a PackMol template.
    if not hasattr(structure, "symmetry"):
        structure.symmetry = None
    structure.save(str(path), overwrite=True)

    if include_conect:
        append_conect_records(path, structure)


def export_pdb(gro_path: Path, output_path: Path, top_path: Path | None = None) -> Structure:
    """Export a PDB file from GROMACS inputs and return the written structure.

    If ``top_path`` is provided, the topology is used to define bonds and the
    output PDB receives ``CONECT`` records.  If ``top_path`` is omitted, the PDB
    is written from coordinates only and no connectivity records are generated.
    """

    structure = load_pdb_source(gro_path=gro_path, top_path=top_path)
    save_pdb(output_path, structure, include_conect=top_path is not None)
    return structure


def append_conect_records(path: Path, structure: Structure) -> None:
    """Append PDB CONECT records from ParmEd bond connectivity.

    ParmEd may write atom coordinates without explicit ``CONECT`` records for
    every bond.  OVITO's "load bonds from file" option needs these records, so
    this function rebuilds them from the topology bond graph and appends them
    before the final ``END`` line.
    """

    adjacency = {atom.idx + 1: set() for atom in structure.atoms}

    for bond in structure.bonds:
        atom1 = bond.atom1.idx + 1
        atom2 = bond.atom2.idx + 1
        adjacency[atom1].add(atom2)
        adjacency[atom2].add(atom1)

    lines = [
        line
        for line in path.read_text(encoding="ascii").splitlines()
        if not line.startswith("CONECT") and line.strip() != "END"
    ]

    for atom_number in sorted(adjacency):
        neighbors = sorted(adjacency[atom_number])
        for start in range(0, len(neighbors), 4):
            chunk = neighbors[start:start + 4]
            # PDB CONECT records store at most four bonded neighbors per line.
            lines.append(
                f"CONECT{atom_number:5d}"
                + "".join(f"{neighbor:5d}" for neighbor in chunk)
            )

    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
