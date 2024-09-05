# PAutodock - Parallelize AutoDock JOBs

PAutodock is a powerful set of scripts designed to parallelize AutoDock jobs, enabling fast screening across multiple CPUs. This tool is particularly useful for researchers and scientists working in computational biology and drug discovery, allowing them to efficiently manage and execute docking simulations.

## Table of Contents

- [License](#license)
- [Dependencies](#dependencies)
- [Changelog](#changelog)
- [Usage](#usage)
- [Contributing](#contributing)
- [Contact](#contact)

## License

PAutodock is distributed under the GPLv3 license. For detailed information on how the license works, please refer to the file "LICENSE" or visit [GNU GPLv3 License](http://www.gnu.org/licenses/gpl-3.0.en.html).

Copyright © Giuseppe Marco Randazzo <gmrandazzo@gmail.com>

## Dependencies

To run PAutodock, you will need the following software installed:

- [AutoDock](http://autodock.scripps.edu/)
- [AutoGrid](http://autodock.scripps.edu/resources/autogrid)
- [AutoDock Vina](http://vina.scripps.edu/)

Ensure that these dependencies are properly installed and accessible from your command line.

## Changelog

- **2024**: Revamp in a more organized form
- **2022**: First time online
- **2017**: Initial release

## Usage

To use PAutodock for parallelizing AutoDock jobs, follow these steps:

1. **Prepare the Receptor and Ligand**:
   - Ensure you have a receptor file (in PDB format) that represents the target protein or enzyme.
   - Prepare a ligand file (also in PDB format) that contains the molecule you want to dock with the receptor.

2. **Create a Multimol2 File**:
   - Prepare a multimol2 file that includes all the ligands you wish to screen.
   - Each ligand in the file must have:
     - **Partial Charges**: Ensure that the ligands have Gasteiger partial charges assigned.
     - **3D Coordinates**: The ligands should be represented in 3D space.
     - **Unique Names**: Each molecule must have a unique name to avoid conflicts during the screening process.

3. **Execute the Command**:
   - Once your receptor and multimol2 file are ready, execute the following command in your terminal:

   ```bash
   python3 main.py --receptor receptor.pdb --ligand ligand.pdb --db multimol2_mol_to_screen.mol2 --wdir="./" --smode=medium --out screening_results.csv
