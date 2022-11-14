# PAutodock - Parallelize AutoDock JOBs

This set of scripts aims to parallelize AutoDock Jobs
and run fast screening on several CPUs.

# License

PAutoDock is distributed under GPLv3 license.
To know more in detail how the license work,
please read the file "LICENSE" or go to
"http://www.gnu.org/licenses/gpl-3.0.en.html"

Copyright Giuseppe Marco Randazzo <gmrandazzo@gmail.com>

# Dependencies

- autodock
- autogrid 
- autodock-vina
- python2 
- MGLTools

# Changelog

- 2022: first time online

- 2017: initial release


# Usage:

```
python3 main.py --receptor receptor.pdb --ligand ligand.pdb --db multimol2_mol_to_screen.mol2 --mglpath <YourPath>/MGLTools-1.5.6/ --wdir="./" --smode=medium --out screening_results.csv
```


