# Task 3 run commands

Run these commands in the WSL/conda environment where `packmol` and `gmx` are
available.  The molecule counts in each `topol.top` must match the PackMol input
and the converted `initial.gro`.

## Pure water system

```bash
cd /mnt/f/polymer_task/systems/water
packmol < packmol.inp
python ../../scripts/task3/packmol_pdb_to_gro.py

gmx grompp -f ../../mdp/em.mdp -c initial.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

gmx grompp -f ../../mdp/task3_nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

gmx grompp -f ../../mdp/task3_npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

gmx grompp -f ../../mdp/task3_prod.mdp -c npt.gro -t npt.cpt -p topol.top -o prod.tpr
gmx mdrun -deffnm prod
```

## Water/ethanol system

```bash
cd /mnt/f/polymer_task/systems/water_ethanol
packmol < packmol.inp
python ../../scripts/task3/packmol_pdb_to_gro.py

gmx grompp -f ../../mdp/em.mdp -c initial.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

gmx grompp -f ../../mdp/task3_nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

gmx grompp -f ../../mdp/task3_npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

gmx grompp -f ../../mdp/task3_prod.mdp -c npt.gro -t npt.cpt -p topol.top -o prod.tpr
gmx mdrun -deffnm prod
```
