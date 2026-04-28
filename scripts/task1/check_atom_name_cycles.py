"""Check whether the 10-chain GRO file is organized as repeated chain cycles.

The single-chain extraction assumes that the original coordinate file is laid
out as 212 atoms per polymer chain. This script verifies that assumption by
checking that each 212-atom block has the same atom-name sequence and that atom
names are unique within each block. It also checks the residue organization
inside each chain:

* 10 residues per chain
* residue names alternate as CLSA-CLSB-CLSA-CLSB-CLSA-CLSB-CLSA-CLSB-CLSA-CLSB
* the two terminal residues contain 22 atoms each
* each middle residue contains 21 atoms
"""

from pathlib import Path
import scripts.utils.gmx_io as gmxio


ROOT = Path(__file__).resolve().parents[2]
INPUT_GRO_PATH = ROOT / "polymer" / "CLS_10chains.gro"
ATOMS_PER_CHAIN = 212
RESIDUES_PER_CHAIN = 10
EXPECTED_RESIDUE_NAMES = ["CLSA", "CLSB"] * 5
EXPECTED_RESIDUE_ATOM_COUNTS = [22] + [21] * 8 + [22]


def require(condition: bool, message: str) -> None:
    """Raise an AssertionError with ``message`` if ``condition`` is false."""
    if not condition:
        raise AssertionError(message)


def residue_groups(chain_atoms) -> list[tuple[str, int]]:
    """Return consecutive residue groups as ``(residue_name, atom_count)``."""
    groups: list[tuple[str, int]] = []

    current_residue = None
    current_name = ""
    current_count = 0

    for atom in chain_atoms:
        residue = atom.residue
        if residue is not current_residue:
            if current_residue is not None:
                groups.append((current_name, current_count))
            current_residue = residue
            current_name = residue.name
            current_count = 1
        else:
            current_count += 1

    if current_residue is not None:
        groups.append((current_name, current_count))

    return groups


def main() -> None:
    """Validate atom-name and residue cycles in 212-atom chain blocks."""
    gro = gmxio.load_gro(INPUT_GRO_PATH)
    total_atoms = len(gro.atoms)

    # A clean split is required before treating each 212-atom block as one
    # chain. Otherwise atom ranges like 1-212 would cut through a molecule.
    require(
        total_atoms % ATOMS_PER_CHAIN == 0,
        f"{total_atoms} atoms cannot be evenly split into {ATOMS_PER_CHAIN}-atom blocks.",
    )

    chain_count = total_atoms // ATOMS_PER_CHAIN
    atom_blocks = [
        gro.atoms[start: start + ATOMS_PER_CHAIN]
        for start in range(0, total_atoms, ATOMS_PER_CHAIN)
    ]
    atom_name_blocks = [[atom.name for atom in block] for block in atom_blocks]

    # Chain 1 is the template. Every later block should match it exactly if the
    # file is organized as repeated copies of the same polymer chain.
    reference = atom_name_blocks[0]
    for chain_index, block in enumerate(atom_name_blocks, start=1):
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

    for chain_index, atom_block in enumerate(atom_blocks, start=1):
        groups = residue_groups(atom_block)
        residue_names = [name for name, _ in groups]
        residue_atom_counts = [count for _, count in groups]

        require(
            len(groups) == RESIDUES_PER_CHAIN,
            f"Chain {chain_index} has {len(groups)} residues, expected {RESIDUES_PER_CHAIN}.",
        )
        require(
            residue_names == EXPECTED_RESIDUE_NAMES,
            f"Chain {chain_index} residue-name pattern is {residue_names}, expected {EXPECTED_RESIDUE_NAMES}.",
        )
        require(
            residue_atom_counts == EXPECTED_RESIDUE_ATOM_COUNTS,
            f"Chain {chain_index} residue atom counts are {residue_atom_counts}, expected {EXPECTED_RESIDUE_ATOM_COUNTS}.",
        )

    print("Atom-name cycle check: PASS")
    print(f"Input file: {INPUT_GRO_PATH}")
    print(f"Total atoms: {total_atoms}")
    print(f"Atoms per cycle: {ATOMS_PER_CHAIN}")
    print(f"Cycle count: {chain_count}")
    print("Atom-name sequence is identical in every 212-atom cycle.")
    print("Atom names are unique within each 212-atom cycle.")
    print(f"Residues per cycle: {RESIDUES_PER_CHAIN}")
    print(f"Residue-name pattern per cycle: {'-'.join(EXPECTED_RESIDUE_NAMES)}")
    print(f"Residue atom counts per cycle: {EXPECTED_RESIDUE_ATOM_COUNTS}")
    print(f"First atom name in each cycle: {reference[0]}")
    print(f"Last atom name in each cycle: {reference[-1]}")


if __name__ == "__main__":
    main()
