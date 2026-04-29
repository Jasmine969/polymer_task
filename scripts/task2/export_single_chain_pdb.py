"""Export task-2 single-chain CLS systems to PDB with topology connectivity.

The PDB is intended for visual inspection.  It is built from the same topology
and coordinates used by ``grompp`` so atom names, residue names, and bonded
connectivity come from the GROMACS inputs rather than distance guessing.
"""

from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

import scripts.utils.pdb_io as pdbio


GRO_PATH = PROJECT_ROOT / "single_chain" / "CLS_single.gro"
# # ===== RAW =======
# TOP_PATH = PROJECT_ROOT / "single_chain" / "topol_raw.top"
# OUTPUT_PATH = PROJECT_ROOT / "single_chain" / "CLS_single_raw.pdb"
# ====== CORRECTED =======
TOP_PATH = PROJECT_ROOT / "single_chain" / "topol_corrected.top"
OUTPUT_PATH = PROJECT_ROOT / "single_chain" / "CLS_single_corrected.pdb"


def main() -> None:
    """Write the corrected single-chain PDB with topology-derived CONECT."""

    structure = pdbio.export_pdb(
        gro_path=GRO_PATH,
        output_path=OUTPUT_PATH,
        top_path=TOP_PATH,
    )

    print(f"Wrote {OUTPUT_PATH}")
    print(f"Topology: {TOP_PATH}")
    print(f"Atoms: {len(structure.atoms)}")
    print(f"Residues: {len(structure.residues)}")
    print(f"Bonds: {len(structure.bonds)}")


if __name__ == "__main__":
    main()
