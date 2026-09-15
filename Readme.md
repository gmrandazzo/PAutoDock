# PAutodock

PAutodock parallelizes molecular docking jobs: it screens a library of
small molecules against a protein receptor with AutoDock Vina and/or
AutoDock4, spreading the calculations over all the CPUs of the machine,
and collects the binding results in a single CSV table.

For every molecule of the input database, PAutodock prepares the ligand
with Open Babel (3D coordinates, Gasteiger partial charges, protonation
at a given pH), prepares the receptor, runs the selected docking engine
in a dedicated working directory and appends one row per molecule with
the binding energies and a geometric pose check. Molecules that already
have results are skipped, so an interrupted screening can simply be
re-run to continue where it stopped.

## Features

- Parallel docking across all available CPUs.
- Two docking engines, selectable independently: AutoDock Vina and
  AutoDock4 (AutoGrid4 + AutoDock4).
- Ligand protonation at a chosen pH with Open Babel (default pH 7.4).
- Receptor preparation with Open Babel by default; the MGLTools
  `prepare_receptor4.py` script remains available with `--mgl ON`.
- Screening of multi-mol2 databases or of a single ligand.
- Resumable screenings: molecules with existing results are skipped.
- Results table with the binding energies of every molecule and the
  distance between the expected centre and the baricentre of the best
  docking pose, as a sanity check.
- Companion command line tools: `pautodock-recover-output` rebuilds a
  results table from a working directory, and
  `pautodock-autogridmap2dx` converts AutoGrid maps to OpenDX format
  for visualization.

## Requirements

- Python >= 3.10
- Open Babel
- AutoDock Vina (to run Vina)
- AutoDock4 and AutoGrid4 (to run AutoDock4)
- MGLTools (only with `--mgl ON`)

## Installation

From PyPI:

    pip install pautodock

From source:

    git clone https://github.com/gmrandazzo/PAutoDock.git
    cd PAutoDock
    poetry install

## Usage

Prepare a receptor in PDB format and a database of ligands as a single
multi-mol2 file (3D coordinates, Gasteiger partial charges and a unique
name per molecule), then run for example:

    cd data/3EML
    pautodock --receptor rec.pdb --cx -9.06364 --cy -7.1446 --cz 55.8626 \
        --db dataset.mol2 --wdir example_calculation \
        --out screening_results.csv --vina ON --atd OFF

The ligands are protonated at pH 7.4 by default; use `--ph` to change
it. The option affects the ligands only, the receptor is never
protonated.

The main command line options:

| Option             | Default   | Description                                                        |
| ------------------ | --------- | ------------------------------------------------------------------ |
| `--receptor`       | required  | receptor PDB file                                                  |
| `--db`             |           | multi-mol2 database to screen                                      |
| `--ligand`         |           | single ligand (PDB or mol2); the grid centre is its baricentre     |
| `--cx --cy --cz`   |           | grid centre (required when `--ligand` is not given)                |
| `--gx --gy --gz`   | 30        | grid size in Angstrom                                              |
| `--vina`           | ON        | run AutoDock Vina                                                  |
| `--atd`            | OFF       | run AutoDock4                                                      |
| `--smode`          | fast      | AutoDock4 screening mode: fast, normal, thorough                   |
| `--exhaustiveness` | 32        | Vina exhaustiveness                                                |
| `--num_modes`      | 18        | number of Vina binding modes                                       |
| `--ph`             | 7.4       | protonation pH of the ligands                                      |
| `--mgl`            | OFF       | prepare the receptor with MGLTools instead of Open Babel           |
| `--out`            | output.txt| results table                                                      |
| `--wdir`           | required  | working directory                                                  |

## Output

The results table is a semicolon-delimited CSV with one row per
molecule: the AutoDock4 free energy terms when `--atd ON` (partition
function, free energy, internal energy, entropy and the cluster
averages), and the average, minimum and maximum Vina binding energy
together with the template-ligand baricentre distance of the best pose.

## Companion tools

Rebuild a lost results table from a working directory:

    pautodock-recover-output --wdir example_calculation --out results.csv

Convert an AutoGrid map to OpenDX (for example for PyMOL):

    pautodock-autogridmap2dx --map path/to/receptor_model.OA.map --dx map.dx

## Changelog

- 2026: version 1.1.0. Ligand protonation at pH 7.4 by default,
  receptor preparation with Open Babel (MGLTools optional), fixed
  coordinate handling and result parsing, strict type checking,
  automated releases.
- 2024: revamp in a more organized form.
- 2022: first release on PyPI.
- 2017: initial release.

## License

PAutodock is distributed under the GNU General Public License v3 or
later. Copyright (C) Giuseppe Marco Randazzo <gmrandazzo@gmail.com>

## Development

See [DEVELOPMENT.md](DEVELOPMENT.md).
