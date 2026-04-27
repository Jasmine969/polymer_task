"""Check whether the 10-chain GRO file is organized as repeated atom-name cycles.

The single-chain extraction assumes that the original coordinate file is laid
out as 212 atoms per polymer chain. This script verifies that assumption by
checking that each 212-atom block has the same atom-name sequence and that atom
names are unique within each block.
"""

from pathlib import Path

import utils.gmx_io as gmxio


ROOT = Path(__file__).resolve().parents[1]
INPUT_GRO_PATH = ROOT / "polymer" / "CLS_10chains.gro"
ATOMS_PER_CHAIN = 212


def require(condition: bool, message: str) -> None:
    """Raise an AssertionError with ``message`` if ``condition`` is false."""
    if not condition:
        raise AssertionError(message)


def main() -> None:
    """Validate that atom names repeat cleanly in 212-atom chain blocks."""
    gro = gmxio.load_gro(INPUT_GRO_PATH)
    total_atoms = len(gro.atoms)

    # A clean split is required before treating each 212-atom block as one
    # chain. Otherwise atom ranges like 1-212 would cut through a molecule.
    require(
        total_atoms % ATOMS_PER_CHAIN == 0,
        f"{total_atoms} atoms cannot be evenly split into {ATOMS_PER_CHAIN}-atom blocks.",
    )

    chain_count = total_atoms // ATOMS_PER_CHAIN
    atom_names = [atom.name for atom in gro.atoms]

    # Split the flattened atom-name list into one block per expected chain.
    blocks = [
        atom_names[start: start + ATOMS_PER_CHAIN]
        for start in range(0, total_atoms, ATOMS_PER_CHAIN)
    ]

    # Chain 1 is the template. Every later block should match it exactly if the
    # file is organized as repeated copies of the same polymer chain.
    reference = blocks[0]
    for chain_index, block in enumerate(blocks, start=1):
        duplicate_names = sorted(
            name for name in set(block) if block.count(name) > 1
        )
        require(
            not duplicate_names,
            f"Chain {chain_index} has duplicated atom names: {duplicate_names}",
        )

        require(
            block == reference,
            f"Atom-name sequence in chain {chain_index} differs from chain 1.",
        )

    print("Atom-name cycle check: PASS")
    print(f"Input file: {INPUT_GRO_PATH}")
    print(f"Total atoms: {total_atoms}")
    print(f"Atoms per cycle: {ATOMS_PER_CHAIN}")
    print(f"Cycle count: {chain_count}")
    print("Atom-name sequence is identical in every 212-atom cycle.")
    print("Atom names are unique within each 212-atom cycle.")
    print(f"First atom name in each cycle: {reference[0]}")
    print(f"Last atom name in each cycle: {reference[-1]}")


if __name__ == "__main__":
    main()
