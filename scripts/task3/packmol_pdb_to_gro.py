"""Convert Task 3 PackMol PDB outputs into GROMACS GRO coordinate files.

PackMol writes coordinates as PDB files.  GROMACS can often read PDB directly,
but using GRO files keeps this workflow consistent with the existing
single-chain coordinate files.  This script expects the PackMol outputs to be:

* ``systems/water/initial.pdb``
* ``systems/water_ethanol/initial.pdb``

It then writes matching ``initial.gro`` files with molecule metadata taken from
the known Task 3 molecule order:

* water system: ``CLS 1 + SOL 1500``
* mixed system: ``CLS 1 + SOL 1050 + ETH 450``

Run it from the repository root after PackMol has created the PDB files:

    python scripts/task3/packmol_pdb_to_gro.py
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import parmed as pmd


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import scripts.utils.gmx_io as gmxio


@dataclass(frozen=True)
class TemplateAtom:
    """Atom metadata copied into the final GRO file."""

    resnr: int
    resname: str
    atomname: str


SYSTEMS = {
    "water": {
        "box_nm": 5.0,
        "blocks": [("CLS", 1), ("SOL", 1500)],
    },
    "water_ethanol": {
        "box_nm": 5.8,
        "blocks": [("CLS", 1), ("SOL", 1050), ("ETH", 450)],
    },
}


def load_gro_template(path: Path) -> list[TemplateAtom]:
    """Read atom names and residue labels from a GRO template with ParmEd.

    The coordinates in the template are irrelevant for this conversion.  Only
    residue numbers, residue names, and atom names are copied so the final GRO
    file matches the topology molecule order expected by GROMACS.
    """

    structure = gmxio.load_gro(path)
    return [
        TemplateAtom(
            resnr=int(atom.residue.number),
            resname=atom.residue.name,
            atomname=atom.name,
        )
        for atom in structure.atoms
    ]


def load_packmol_coordinates(path: Path) -> list[tuple[float, float, float]]:
    """Read PackMol PDB coordinates with ParmEd and convert Angstrom to nm."""

    structure = pmd.load_file(str(path))
    if structure.coordinates is None:
        raise ValueError(f"{path} does not contain readable coordinates.")

    return [
        (float(x_angstrom) / 10.0, float(y_angstrom) / 10.0, float(z_angstrom) / 10.0)
        for x_angstrom, y_angstrom, z_angstrom in structure.coordinates
    ]


def build_templates() -> dict[str, list[TemplateAtom]]:
    """Load molecule templates in the same atom order used by PackMol."""

    water = [
        TemplateAtom(1, "SOL", "OW"),
        TemplateAtom(1, "SOL", "HW1"),
        TemplateAtom(1, "SOL", "HW2"),
    ]
    return {
        "CLS": load_gro_template(ROOT / "single_chain" / "CLS_single.gro"),
        "SOL": water,
        "ETH": load_gro_template(ROOT / "solvent" / "ethanol.gro"),
    }


def expanded_metadata(blocks: list[tuple[str, int]], templates: dict[str, list[TemplateAtom]]) -> list[TemplateAtom]:
    """Expand molecule templates to the full atom metadata list for one system."""

    metadata: list[TemplateAtom] = []
    next_resnr = 1

    for molecule_name, count in blocks:
        template = templates[molecule_name]
        residue_offset = max(atom.resnr for atom in template)

        for _ in range(count):
            # Keep the polymer internal residue pattern, while assigning each
            # solvent molecule a distinct residue number for easier inspection.
            for atom in template:
                metadata.append(
                    TemplateAtom(
                        resnr=next_resnr + atom.resnr - 1,
                        resname=atom.resname,
                        atomname=atom.atomname,
                    )
                )
            next_resnr += residue_offset

    return metadata


def write_gro(path: Path, metadata: list[TemplateAtom], coords: list[tuple[float, float, float]], box_nm: float) -> None:
    """Write a complete GRO file after validating atom counts.

    ParmEd is used above to read the source coordinates and template GRO files.
    The final write is kept local because this file deliberately rewrites the
    metadata from the known Task 3 molecule order instead of preserving PackMol's
    PDB residue labels verbatim.
    """

    if len(metadata) != len(coords):
        raise ValueError(
            f"{path}: metadata atom count {len(metadata)} does not match "
            f"PDB coordinate count {len(coords)}"
        )

    lines = ["Task 3 initial system converted from PackMol PDB", f"{len(coords):5d}"]
    for atomnr, (atom, coord) in enumerate(zip(metadata, coords), start=1):
        x_nm, y_nm, z_nm = coord
        lines.append(
            f"{atom.resnr % 100000:5d}{atom.resname[:5]:<5}{atom.atomname[:5]:>5}"
            f"{atomnr % 100000:5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}"
        )
    lines.append(f"{box_nm:10.5f}{box_nm:10.5f}{box_nm:10.5f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def validate_coordinates_inside_box(
    system_name: str,
    coords: list[tuple[float, float, float]],
    box_nm: float,
) -> None:
    """Check that PackMol coordinates are compatible with the target box.

    This catches interrupted or malformed PackMol outputs before they overwrite
    an existing GRO file.  A tiny tolerance is allowed for text rounding at the
    box boundary.
    """

    tolerance_nm = 0.001
    for atom_index, coord in enumerate(coords, start=1):
        if any(value < -tolerance_nm or value > box_nm + tolerance_nm for value in coord):
            coord_text = ", ".join(f"{value:.4f}" for value in coord)
            raise ValueError(
                f"{system_name} atom {atom_index} coordinate ({coord_text}) nm "
                f"is outside the expected 0-{box_nm:.3f} nm box."
            )


def main() -> None:
    """Convert any available Task 3 PackMol PDB files."""

    templates = build_templates()
    converted = 0

    for system_name, spec in SYSTEMS.items():
        system_dir = ROOT / "systems" / system_name
        pdb_path = system_dir / "initial.pdb"
        if not pdb_path.exists():
            print(f"Skipping {system_name}: {pdb_path} does not exist yet.")
            continue

        coords = load_packmol_coordinates(pdb_path)
        validate_coordinates_inside_box(system_name, coords, spec["box_nm"])
        metadata = expanded_metadata(spec["blocks"], templates)
        write_gro(system_dir / "initial.gro", metadata, coords, spec["box_nm"])
        print(f"Converted {pdb_path} -> {system_dir / 'initial.gro'}")
        converted += 1

    print(f"Converted {converted} system(s).")


if __name__ == "__main__":
    main()
