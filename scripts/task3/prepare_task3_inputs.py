"""Prepare clean Task 3 topology and packing inputs.

This script converts the successful ACPYPE ethanol output into files that are
easy to include in the two solvent systems required by Task 3. It deliberately
keeps the original ``solvent/smiles_molecule.acpype`` directory untouched and
creates reproducible topology, coordinate, PackMol, and MDP input files.
"""

from pathlib import Path
import math
import sys


ROOT = Path(__file__).resolve().parents[2]
ACPYPE_DIR = ROOT / "solvent" / "smiles_molecule.acpype"

sys.path.insert(0, str(ROOT))

import scripts.utils.gmx_io as gmxio
import scripts.utils.pdb_io as pdbio


def write_clean_ethanol_itp() -> None:
    """Create ``solvent/ethanol.itp`` from ACPYPE's ethanol ITP.

    Input:
        ``solvent/smiles_molecule.acpype/smiles_molecule_GMX.itp``.
    Output:
        ``solvent/ethanol.itp``.

    The output starts at ``[ moleculetype ]`` so atom types can be collected in
    ``solvent/forcefield_gaff_spce.itp`` before any molecule definitions.
    """

    source = ACPYPE_DIR / "smiles_molecule_GMX.itp"
    target = ROOT / "solvent" / "ethanol.itp"
    lines = source.read_text(encoding="utf-8").splitlines()

    cleaned: list[str] = [
        "; Ethanol molecule topology prepared from ACPYPE smiles_molecule_GMX.itp.",
        "; Atom types are defined in solvent/forcefield_gaff_spce.itp.",
        "",
    ]

    # GROMACS requires all [ atomtypes ] sections to appear before any
    # [ moleculetype ], so the ethanol molecule file starts at [ moleculetype ].
    keep = False
    for line in lines:
        stripped = line.strip()
        if stripped.lower() == "[ moleculetype ]":
            keep = True
        # 略过[ moleculetype ]之前的行
        if not keep:
            continue
        # molecule name由smiles_molecule（actype默认）替换为ETH
        if stripped.startswith("smiles_molecule") and "3" in stripped:
            cleaned.append("ETH              3")
        # 残基名由UNL替换为ETH
        elif " UNL " in line:
            cleaned.append(line.replace(" UNL ", " ETH "))
        # 原样拷贝
        else:
            cleaned.append(line)

    target.write_text("\n".join(cleaned).rstrip() + "\n", encoding="utf-8")


def write_clean_ethanol_coordinates() -> None:
    """Create ETH-named ethanol coordinate files from ACPYPE's GRO.

    Input:
        ``solvent/smiles_molecule.acpype/smiles_molecule_GMX.gro``.
    Outputs:
        ``solvent/ethanol.gro`` and ``solvent/ethanol.pdb``.

    The GRO file is a reusable molecule coordinate file. The PDB file is the
    PackMol template used in the water/ethanol system.
    """

    source = ACPYPE_DIR / "smiles_molecule_GMX.gro"
    ethanol = gmxio.load_gro(source)
    relabel_structure(ethanol, resname="ETH")
    ethanol_gro = ROOT / "solvent" / "ethanol.gro"
    gmxio.save_gro(ethanol_gro, ethanol)

    pdbio.export_pdb(
        gro_path=ethanol_gro,
        output_path=ROOT / "solvent" / "ethanol.pdb",
    )


def relabel_structure(structure, resname: str) -> None:
    """Set all residue names in a ParmEd structure to ``resname`` in place."""

    for residue in structure.residues:
        residue.name = resname


def write_origin_centered_cls_pdb() -> None:
    """
    Create a centered single-chain PDB for PackMol placement.
    单链是从多条链中提取出来的，y坐标很大，位置很偏，把它挪到原点附近（以原点为中心），以便给PackMol使用
    """

    # Input: single_chain/CLS_single.gro.
    # Output: single_chain/CLS_single_packmol.pdb.
    cls_chain = gmxio.load_gro(ROOT / "single_chain" / "CLS_single.gro")
    coordinates = cls_chain.coordinates
    center = coordinates.mean(axis=0)

    # PackMol works more reliably when the molecule starts near the origin; the
    # final position is controlled by the packmol.inp inside-box constraints.
    cls_chain.coordinates = coordinates - center
    relabel_structure(cls_chain, resname="CLS")
    pdbio.save_pdb(
        ROOT / "single_chain" / "CLS_single_packmol.pdb",
        cls_chain,
    )


def write_water_pdb_template() -> None:
    """Write ``solvent/water.pdb`` from SPC/E SETTLE geometry.

    Reference input:
        ``solvent/water.itp`` supplies the SETTLE geometry:
        ``doh = 0.1000 nm`` and ``dhh = 0.1630 nm``.
    Output:
        ``solvent/water.pdb``.

    The simple ``solvent/water_corrected.itp`` topology is maintained manually.
    This function only computes the PackMol water coordinate template. PDB uses
    Angstrom units, so the SETTLE distances are converted from nm to Angstrom.
    """

    oxygen_hydrogen_nm = 0.1000
    hydrogen_hydrogen_nm = 0.1630
    oxygen_hydrogen_angstrom = oxygen_hydrogen_nm * 10.0
    hydrogen_hydrogen_angstrom = hydrogen_hydrogen_nm * 10.0

    # Place OW at the origin and HW1 on the +x axis. HW2 is computed from the
    # two SETTLE distances, so both O-H distances and the H-H distance match the
    # topology geometry.
    hw2_x = (
        2.0 * oxygen_hydrogen_angstrom**2 - hydrogen_hydrogen_angstrom**2
    ) / (2.0 * oxygen_hydrogen_angstrom)
    hw2_y = math.sqrt(oxygen_hydrogen_angstrom**2 - hw2_x**2)

    (ROOT / "solvent" / "water.pdb").write_text(
        "\n".join(
            [
                "REMARK SPC/E water molecule for PackMol",
                f"HETATM    1   OW SOL A   1    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           O",
                f"HETATM    2  HW1 SOL A   1    {oxygen_hydrogen_angstrom:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           H",
                f"HETATM    3  HW2 SOL A   1    {hw2_x:8.3f}{hw2_y:8.3f}{0.0:8.3f}  1.00  0.00           H",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def write_forcefield_file() -> None:
    """Merge force-field defaults and atom types used by Task 3 topologies.

    Reference inputs:
        ``single_chain/topol_corrected.top`` supplies the GAFF defaults and the
        polymer atom types. ACPYPE's original ethanol ITP supplies ethanol atom
        types. ``solvent/water.itp`` supplies the OW/HW atom types.
    Output:
        ``solvent/forcefield_gaff_spce.itp``.

    Keeping this header separate avoids GROMACS directive-order errors from
    including molecule files that also contain ``[ atomtypes ]``.
    """

    default_source = ROOT / "single_chain" / "topol_corrected.top"
    atomtype_sources = [
        default_source,
        ACPYPE_DIR / "smiles_molecule_GMX.itp",
        ROOT / "solvent" / "water.itp",
    ]

    defaults_comment, defaults_line = gmxio.read_defaults(default_source)
    atomtype_rows = gmxio.merge_atomtypes(atomtype_sources)

    gmxio.write_forcefield_header(
        ROOT / "solvent" / "forcefield_gaff_spce.itp",
        defaults_comment=defaults_comment,
        defaults_line=defaults_line,
        atomtype_rows=atomtype_rows,
        header_comments=[
            "; Minimal GAFF/SPC/E force-field header for this system.",
            "; Generated by scripts/task3/prepare_task3_inputs.py.",
            "; Keep this file included before any molecule .itp file.",
        ],
    )


def write_task3_mdp_files() -> None:
    """Create Task 3 MDP variants from the base MDP files.

    Inputs:
        ``mdp/nvt.mdp``, ``mdp/npt.mdp``, and ``mdp/nvt_prod.mdp``.
    Outputs:
        ``mdp/task3_nvt.mdp``, ``mdp/task3_npt.mdp``, and
        ``mdp/task3_prod.mdp``.

    The generated variants use a single ``System`` temperature-coupling group
    and add explicit velocity-generation or continuation settings.
    """

    replacements = {
        "tc-grps": "tc-grps              = System",
        "tau_t": "tau_t                = 0.1",
        "ref_t": "ref_t                = 300",
    }
    for source_name, target_name in [
        ("nvt.mdp", "task3_nvt.mdp"),
        ("npt.mdp", "task3_npt.mdp"),
        ("nvt_prod.mdp", "task3_prod.mdp"),
    ]:
        source_lines = (ROOT / "mdp" / source_name).read_text(encoding="utf-8").splitlines()
        target_lines: list[str] = []
        for line in source_lines:
            key = line.split("=", maxsplit=1)[0].strip() if "=" in line else ""
            target_lines.append(replacements.get(key, line))

        # NVT starts from minimized coordinates, so it needs generated
        # velocities.  Later stages should continue from the previous
        # checkpoint and keep the existing velocities.
        if target_name == "task3_nvt.mdp":
            target_lines.extend(["", "; Velocity generation", "gen_vel              = yes", "gen_temp             = 300", "gen_seed             = -1"])
        else:
            target_lines.extend(["", "; Continuation from previous checkpoint", "continuation          = yes", "gen_vel              = no"])

        (ROOT / "mdp" / target_name).write_text("\n".join(target_lines) + "\n", encoding="utf-8")


def main() -> None:
    """Generate all Task 3 helper files."""

    write_clean_ethanol_itp()
    write_clean_ethanol_coordinates()
    write_origin_centered_cls_pdb()
    write_water_pdb_template()
    write_forcefield_file()
    write_task3_mdp_files()

    print("Task 3 input files prepared.")
    print("Next: run packmol inside systems/water and systems/water_ethanol.")


if __name__ == "__main__":
    main()
