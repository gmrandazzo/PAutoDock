#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""multimol2op.py

This file is part of PAutoDock.
Copyright (C) 2020 Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
PAutoDock is distributed under GPLv3 license.
To know more in detail how the license work,
please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"

Provides the basic operation for MOL2 files in parallel

"""

import os
import shutil
from pathlib import Path


def read_molname(filemol2: str):
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


def get_mol2_name(filemol2: str):
    molname = read_molname("tmp_")
    filename = molname + ".mol2"
    cc = 1
    while True:
        if os.path.isfile(filename) is False:
            break
        else:
            filename = molname + "_" + str(cc) + ".mol2"
            cc = cc + 1
    return filename


def split_mol2(mmol2, path="./"):
    mol2splitted = []
    fmol2 = open(mmol2, "r")
    ftmp = open("tmp_", "w")
    firstmol = True
    for line in fmol2:
        if line.strip().lower() == "@<TRIPOS>MOLECULE".lower():
            if ftmp.closed is False and firstmol is False:
                ftmp.close()
                filename = get_mol2_name("tmp_")
                mpath = str(Path(path + "/" + filename).absolute())
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
        filename = get_mol2_name("tmp_")
        mpath = str(Path(path + "/" + filename).absolute())
        if Path(mpath).exists() is False:
            shutil.move("tmp_", mpath)
            mol2splitted.append(mpath)
        else:
            os.remove("tmp_")
    fmol2.close()
    return mol2splitted
