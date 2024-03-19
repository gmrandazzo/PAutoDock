"""recover_output.py
Provides a commandline script to recover the output result of a virtualscreening.
"""

#!/usr/bin/env python3

import sys
import os
from pathlib import Path
import argparse
from adparallel import ADParallel


def main():
    """
    main.py
    """
    p = argparse.ArgumentParser()
    p.add_argument('--wdir', default=None, type=str, help='work directory')
    p.add_argument('--out', default=None, type=str, help='screening output')
    p.add_argument('--ligand', default=None, type=str, help='ligand')
    args = p.parse_args(sys.argv[1:])

    if args.wdir is None or args.out is None:
        print("\nUsage: %s --receptor [input pdb]" % sys.argv[0])
        print("                --wdir [work path]")
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


if __name__ == '__main__':
    main()
