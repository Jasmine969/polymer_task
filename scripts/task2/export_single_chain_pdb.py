"""Export task-2 single-chain CLS systems to PDB with ParmEd.

The PDB is intended for visual inspection. It is built from the same topology
and coordinates used by ``grompp`` so atom names, residue names, and bonded
connectivity come from the GROMACS inputs rather than distance guessing.
"""

from pathlib import Path

from parmed.gromacs import GromacsTopologyFile

ROOT = Path(__file__).resolve().parents[2]
GRO_PATH = ROOT / "single_chain" / "CLS_single.gro"
# # ===== RAW =======
# TOP_PATH = ROOT / "single_chain" / "topol_raw.top"
# OUTPUT_PATH = ROOT / "single_chain" / "CLS_single_raw.pdb"
# ====== CORRECTED =======
TOP_PATH=ROOT / "single_chain" / "topol_corrected.top"
OUTPUT_PATH=ROOT / "single_chain" / "CLS_single_corrected.pdb"


def main() -> None:
    export_pdb(TOP_PATH, GRO_PATH, OUTPUT_PATH)


def export_pdb(top_path: Path, gro_path: Path, output_path: Path) -> None:
    """Write one PDB file from a GROMACS top file and the single-chain GRO.

    Parameters
    ----------
    top_path:
        Main GROMACS topology file. In this project it is either the raw
        topology for reproducing the error or the corrected topology for final
        validation.
    gro_path:
        Single-chain coordinate file shared by both topologies.
    output_path:
        PDB file to write. Existing files are overwritten so repeated runs keep
        the visual-inspection files synchronized with the ITP files.
    """
    structure = GromacsTopologyFile(str(top_path), xyz=str(gro_path))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    structure.save(str(output_path), overwrite=True)
    append_conect_records(output_path, structure)

    print(f"Wrote {output_path}")
    print(f"Topology: {top_path}")
    print(f"Atoms: {len(structure.atoms)}")
    print(f"Residues: {len(structure.residues)}")
    print(f"Bonds: {len(structure.bonds)}")


def append_conect_records(path: Path, structure) -> None:
    """Append PDB CONECT records from ParmEd bond connectivity.

    ParmEd may write atom coordinates without explicit ``CONECT`` records for
    every bond. OVITO's "load bonds from file" option needs these records, so
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


if __name__ == "__main__":
    main()
