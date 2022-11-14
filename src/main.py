"""main.py
Provides the basic main function to run autodock screening
in parallel.
"""

#!/usr/bin/env python3

import sys
import argparse
from adparallel import ADParallel


def main():
    """
    main.py: with that you run the program.
    """
    p = argparse.ArgumentParser()
    p.add_argument('--receptor',
                   default=None, type=str,
                   help='receptor')
    p.add_argument('--ligand', default=None, type=str, help='ligand')
    p.add_argument('--cx', default=None, type=float, help='center x')
    p.add_argument('--cy', default=None, type=float, help='center y')
    p.add_argument('--cz', default=None, type=float, help='center z')
    p.add_argument('--gx', default=30, type=int, help='grid x size')
    p.add_argument('--gy', default=30, type=int, help='grid y size')
    p.add_argument('--gz', default=30, type=int, help='grid z size')
    p.add_argument('--db', default=None, type=str, help='multimol2 database')
    p.add_argument('--mglpath', default=None, type=str, help='MGLTools path')
    p.add_argument('--atdpath', default=None, type=str, help='Autodock path')
    p.add_argument('--wdir', default=None, type=str, help='work directory')
    p.add_argument('--smode', default="fast", type=str, help='screening mode')
    p.add_argument('--out', default=None, type=str, help='screening output')
    p.add_argument('--atd', default=True, type=str, help='Autodock ON')
    p.add_argument('--vina', default=True, type=str, help='Vina ON')
    args = p.parse_args(sys.argv[1:])

    if args.receptor is None or args.ligand is None and args.cx is None:
        print("\nUsage: %s  --receptor [input pdb]" % sys.argv[0])
        print("                --ligand [input pdb/mol2]")
        print("                --cx [center x]")
        print("                --cy [center y]")
        print("                --cz [center z]")
        print("                --mglpath [path mgltools]")
        print("                --atdpath [path autodock]")
        print("                --db [multi mol2]")
        print("                --wdir [work path]")
        print("                --smode=[slow,medium,fast]")
        print("                --gx=[grid x size]")
        print("                --gy=[grid y size]")
        print("                --gz=[grid z size]")
        print("                --out [screening output]")
        print("                --atd [ON;OFF]")
        print("                --vina [ON;OFF]")
    else:
        # Load import pdb; pdb.set_trace()#
        dock = ADParallel("", args.mglpath, args.receptor,
                          args.ligand, args.db, args.wdir)
        if args.ligand is None:
            dock.cx = args.cx
            dock.cy = args.cy
            dock.cz = args.cz

        if args.atd == "OFF":
            dock.atd = False

        if args.vina == "OFF":
            dock.vina = False

        dock.speed = args.smode
        dock.gsize_x = args.gx
        dock.gsize_y = args.gy
        dock.gsize_z = args.gz
        dock.VS(args.out)
    return 0


if __name__ == '__main__':
    main()
