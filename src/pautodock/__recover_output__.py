#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""recover_output.py

This file is part of PAutoDock.
Copyright (C) 2020 Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
PAutoDock is distributed under GPLv3 license.
To know more in detail how the license work,
please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"

Provides a commandline script to recover the output result of a virtualscreening.
"""

import argparse
import os
import sys
from pathlib import Path

from pautodock.adparallel import ADParallel


def main():
    """
    main.py
    """
    p = argparse.ArgumentParser()
    p.add_argument("--wdir", default=None, type=str, help="work directory")
    p.add_argument("--out", default=None, type=str, help="screening output")
    p.add_argument("--ligand", default=None, type=str, help="ligand")
    args = p.parse_args(sys.argv[1:])

    if args.wdir is None or args.out is None:
        print("\nUsage: %s --receptor [input pdb]" % sys.argv[0])
        print("                --wdir [work path]")
        print("                --ligand [ligand PDB]")
        print("                --out [screening output]")
    else:
        dock = ADParallel("", None, None, args.ligand, None, args.wdir)
        vinalogout = []
        dpfout = []
        mnames = []
        for root, directories, filenames in os.walk(str(Path(args.wdir).absolute())):
            for filename in filenames:
                if "vina_log.txt" in filename:
                    vinalogout.append(os.path.join(root, filename))
                    mname = str(Path(vinalogout[-1]).parents[0].name)
                    if mname not in mnames:
                        mnames.append(mname)
                elif "ind.dlg" in filename:
                    dpfout.append(os.path.join(root, filename))
                    mname = str(Path(dpfout[-1]).parents[0].name)
                    if mname not in mnames:
                        mnames.append(mname)
                else:
                    continue
        dock.GenVSOutput(vinalogout, dpfout, mnames, args.out)
        return 0


if __name__ == "__main__":
    main()
