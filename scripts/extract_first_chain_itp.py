"""Extract the first-chain topology from the provided 10-chain CLS ITP."""

from pathlib import Path

import utils.gmx_io as gmxio


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_ITP = PROJECT_ROOT / "topology" / "CLS_10chains.itp"
OUTPUT_ITP = PROJECT_ROOT / "single_chain" / "CLS_single.itp"
ATOMS_PER_CHAIN = 212
MOLECULE_NAME = "CLS"


def main() -> None:
    """Create ``single_chain/CLS_single.itp`` from the first 212 atoms."""
    topology = gmxio.load_itp(INPUT_ITP, parametrize=False)
    single_chain = gmxio.first_atoms(topology, ATOMS_PER_CHAIN)

    gmxio.save_itp(OUTPUT_ITP, single_chain, molecule_name=MOLECULE_NAME)

    print(f"Wrote {OUTPUT_ITP}")
    print(f"Atoms: {len(single_chain.atoms)}")
    print(f"Bonds: {len(single_chain.bonds)}")
    print(f"Pairs: {len(single_chain.adjusts)}")
    print(f"Angles: {len(single_chain.angles)}")
    print(f"Dihedrals: {len(single_chain.dihedrals)}")


if __name__ == "__main__":
    main()
