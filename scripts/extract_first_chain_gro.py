"""Extract the first polymer chain from the provided 10-chain CLS GRO file."""

from pathlib import Path

import utils.gmx_io as gmxio


ROOT = Path(__file__).resolve().parents[1]
INPUT_GRO_PATH = ROOT / "polymer" / "CLS_10chains.gro"
OUTPUT_GRO_PATH = ROOT / "single_chain" / "CLS_single.gro"

# From the original gro/itp pattern: 10 chains contain 2120 atoms in total,
# so one chain contains 212 atoms. The first chain is residue 1-10, atom 1-212.
ATOMS_PER_CHAIN = 212
EXPECTED_RESIDUE_START = 1
EXPECTED_RESIDUE_END = 10
EXPECTED_NEXT_RESIDUE = 11
EXPECTED_NEXT_ATOM = 213


def require(condition: bool, message: str) -> None:
    """Raise an AssertionError with ``message`` if ``condition`` is false."""
    if not condition:
        raise AssertionError(message)


def check_single_chain(original, single) -> None:
    """Validate that the extracted coordinate file is a first-chain structure."""
    require(
        len(single.atoms) == ATOMS_PER_CHAIN,
        f"Single-chain atom count is {len(single.atoms)}, expected {ATOMS_PER_CHAIN}.",
    )
    require(
        single.box is not None and (single.box == original.box).all(),
        "Single-chain box differs from the original box.",
    )

    atom_numbers = [atom.number for atom in single.atoms]
    residue_numbers = [atom.residue.number for atom in single.atoms]
    atom_names = [atom.name for atom in single.atoms]
    original_names = [atom.name for atom in original.atoms[:ATOMS_PER_CHAIN]]

    require(atom_names == original_names, "Single-chain atom names do not match original atoms 1-212.")
    require(atom_numbers == list(range(1, ATOMS_PER_CHAIN + 1)), "Atom numbers are not continuous from 1 to 212.")
    require(min(residue_numbers) == EXPECTED_RESIDUE_START, "Unexpected first residue number.")
    require(max(residue_numbers) == EXPECTED_RESIDUE_END, "Unexpected last residue number.")

    # The next atom in the original file should be the start of the second chain.
    next_atom = original.atoms[ATOMS_PER_CHAIN]
    require(next_atom.number == EXPECTED_NEXT_ATOM, "The next original atom is not atom 213.")
    require(next_atom.residue.number == EXPECTED_NEXT_RESIDUE, "The next original residue is not residue 11.")
    require(next_atom.name == "H1", "The next original atom is not the expected new-chain H1 atom.")

    print("Single-chain .gro self-consistency check: PASS")
    print(f"Original atoms: {len(original.atoms)}")
    print(f"Single-chain atoms: {len(single.atoms)}")
    print(f"Single-chain atom range: {atom_numbers[0]}-{atom_numbers[-1]}")
    print(f"Single-chain residue range: {min(residue_numbers)}-{max(residue_numbers)}")
    print(
        "Next chain starts at "
        f"residue {next_atom.residue.number} {next_atom.residue.name}, "
        f"atom {next_atom.number} {next_atom.name}"
    )
    print("Atom names match original atoms 1-212.")
    print("Box is preserved from the original .gro file.")


def main() -> None:
    """Create ``single_chain/CLS_single.gro`` from the first 212 atoms."""
    original = gmxio.load_gro(INPUT_GRO_PATH)
    single_chain = gmxio.first_atoms(original, ATOMS_PER_CHAIN)

    gmxio.save_gro(OUTPUT_GRO_PATH, single_chain)
    print(f"Wrote {OUTPUT_GRO_PATH} with {ATOMS_PER_CHAIN} atoms.")

    # Re-read the generated file so the check validates the actual written output.
    written = gmxio.load_gro(OUTPUT_GRO_PATH)
    check_single_chain(original=original, single=written)


if __name__ == "__main__":
    main()
