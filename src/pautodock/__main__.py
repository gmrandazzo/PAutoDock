#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""main.py

This file is part of PAutoDock.
Copyright (C) 2020 Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
PAutoDock is distributed under GPLv3 license.
To know more in detail how the license work,
please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"


Provides the basic main function to run autodock screening
in parallel.
"""

import argparse
import sys
from dataclasses import dataclass
from typing import Optional

from pautodock.adparallel import ADParallel


@dataclass
class DockingConfig:
    """Configuration class for molecular docking parameters."""

    receptor: str
    ligand: Optional[str]
    db: Optional[str]
    wdir: str
    center_x: Optional[float]
    center_y: Optional[float]
    center_z: Optional[float]
    grid_x: int = 30
    grid_y: int = 30
    grid_z: int = 30
    screening_mode: str = "fast"
    output_path: str = "output.txt"
    autodock_enabled: bool = False
    vina_enabled: bool = True
    vina_exhaustiveness: int = 32
    vina_num_modes: int = 18


def parse_arguments() -> DockingConfig:
    """
    Parse command line arguments for molecular docking configuration.

    Returns:
        DockingConfig: Configuration object with all docking parameters
    """
    parser = argparse.ArgumentParser(
        description="Molecular Docking Program",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required arguments group
    required = parser.add_argument_group("Required Arguments")
    required.add_argument("--receptor", type=str, help="Path to receptor PDB file")
    required.add_argument("--wdir", type=str, help="Working directory path")

    # Optional arguments
    parser.add_argument("--ligand", type=str, help="Path to ligand PDB file")
    parser.add_argument("--db", type=str, help="Path multimol2 to screen")

    # Grid center coordinates
    grid_group = parser.add_argument_group("Grid Configuration")
    grid_group.add_argument("--cx", type=float, help="Grid center X coordinate")
    grid_group.add_argument("--cy", type=float, help="Grid center Y coordinate")
    grid_group.add_argument("--cz", type=float, help="Grid center Z coordinate")
    grid_group.add_argument("--gx", type=int, default=30, help="Grid size X")
    grid_group.add_argument("--gy", type=int, default=30, help="Grid size Y")
    grid_group.add_argument("--gz", type=int, default=30, help="Grid size Z")

    # Docking configuration
    dock_group = parser.add_argument_group("Docking Configuration")
    dock_group.add_argument(
        "--smode",
        type=str,
        default="fast",
        choices=["fast", "normal", "thorough"],
        help="Screening mode",
    )
    dock_group.add_argument(
        "--out", type=str, default="output.txt", help="Output file path"
    )
    dock_group.add_argument(
        "--atd", type=str, default="OFF", choices=["ON", "OFF"], help="Enable Autodock"
    )
    dock_group.add_argument(
        "--vina", type=str, default="ON", choices=["ON", "OFF"], help="Enable Vina"
    )
    dock_group.add_argument(
        "--exhaustiveness", type=int, default=32, help="Vina exhaustiveness parameter"
    )
    dock_group.add_argument(
        "--num_modes", type=int, default=18, help="Number of binding modes to generate"
    )

    args = parser.parse_args(sys.argv[1:])

    # Validate required arguments
    if not args.receptor or not args.wdir:
        parser.print_help()
        sys.exit(1)

    # Validate that either ligand or database is provided
    if not args.ligand and not args.db:
        parser.error("Either --ligand or --db must be provided")

    # Validate grid center coordinates when no ligand is provided
    if not args.ligand and not all([args.cx, args.cy, args.cz]):
        parser.error(
            "Grid center coordinates (--cx, --cy, --cz) are required when no ligand is provided"
        )

    return DockingConfig(
        receptor=args.receptor,
        ligand=args.ligand,
        db=args.db,
        wdir=args.wdir,
        center_x=args.cx,
        center_y=args.cy,
        center_z=args.cz,
        grid_x=args.gx,
        grid_y=args.gy,
        grid_z=args.gz,
        screening_mode=args.smode,
        output_path=args.out,
        autodock_enabled=args.atd == "ON",
        vina_enabled=args.vina == "ON",
        vina_exhaustiveness=args.exhaustiveness,
        vina_num_modes=args.num_modes,
    )


def main() -> int:
    """
    Main function to run the poautodock program.

    Returns:
        int: Exit code (0 for success, non-zero for failure)
    """
    try:
        config = parse_arguments()

        # Initialize docking
        dock = ADParallel(
            receptor=config.receptor,
            ligand=config.ligand,
            db=config.db,
            wpath=config.wdir,
        )

        # Configure docking parameters
        if config.ligand is None:
            dock.cx = config.center_x
            dock.cy = config.center_y
            dock.cz = config.center_z

        dock.atd = config.autodock_enabled
        dock.vina = config.vina_enabled
        dock.speed = config.screening_mode
        dock.gsize_x = config.grid_x
        dock.gsize_y = config.grid_y
        dock.gsize_z = config.grid_z

        # Run virtual screening
        dock.virtual_screening(config.output_path)
        return 0

    except Exception as err:
        print(f"Error: {str(err)}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
