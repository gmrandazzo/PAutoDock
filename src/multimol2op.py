"""multimol2op.py
Provides the basic operation for MOL2 files in parallel

"""

#!/usr/bim/env python3

import sys
import shutil
import os
from pathlib import Path


def ReadMol2Name(filemol2):
    f = open(filemol2, "r")
    molname = None
    next_is_name = False
    for line in f:
        if line.strip().lower() == "@<TRIPOS>MOLECULE".lower():
            next_is_name = True
        else:
            if next_is_name is True:
                molname = line.strip()
                next_is_name = False
    f.close()
    return molname


def GetMol2Name(filemol2):
    molname = ReadMol2Name("tmp_")
    filename = molname+".mol2"
    cc = 1
    while True:
        if os.path.isfile(filename) is False:
            break
        else:
            filename = molname+"_"+str(cc)+".mol2"
            cc = cc + 1
    return filename


def SplitMol2(mmol2, path="./"):
    mol2splitted = []
    fmol2 = open(mmol2, "r")
    ftmp = open("tmp_", "w")
    firstmol = True
    for line in fmol2:
        if line.strip().lower() == "@<TRIPOS>MOLECULE".lower():
            if ftmp.closed is False and firstmol is False:
                ftmp.close()
                filename = GetMol2Name("tmp_")
                mpath = str(Path(path+"/"+filename).absolute())
                if Path(mpath).exists() is False:
                    shutil.move("tmp_", mpath)
                    mol2splitted.append(mpath)
                else:
                    os.remove("tmp_")
                ftmp = open("tmp_", "w")
            firstmol = False
            ftmp.write(line)
        else:
            ftmp.write(line)

    if ftmp.closed is False:
        ftmp.close()
        filename = GetMol2Name("tmp_")
        mpath = str(Path(path+"/"+filename).absolute())
        if Path(mpath).exists() is False:
            shutil.move("tmp_", mpath)
            mol2splitted.append(mpath)
        else:
            os.remove ("tmp_")
    fmol2.close()
    return mol2splitted


def main():
    if len(sys.argv) < 2:
        print("\n Usage: ", sys.argv[0], "\t <Mulimol2>\n")
    else:
        SplitMol2(sys.argv[1])


if __name__ == "__main__":
    main()
