#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""molop.py

This file is part of PAutoDock.
Copyright (C) 2020 Giuseppe Marco Randazzo <gmrandazzo@gmail.com>
PAutoDock is distributed under GPLv3 license.
To know more in detail how the license work,
please read the file "LICENSE" or
go to "http://www.gnu.org/licenses/gpl-3.0.en.html"

Provides the basic operation for molecular files.

"""

import logging
import subprocess
from pathlib import Path


def nsplit(s, delim=None):
    return [x for x in s.split(delim) if x]


def get_mol_baricentre(mol: str) -> tuple:
    """
    Get molecular baricentre from a molecule
    """
    cc = [0.0, 0.0, 0.0]
    n = 0
    if ".pdbqt" in mol:
        j = -5
    else:
        j = -4
    with open(mol, "r", encoding="utf-8") as f:
        for line in f:
            if ("ATOM" in line or "HETATM" in line) and "REMARK" not in line:
                try:
                    # Get the results from back
                    items = list(filter(None, str.split(line.strip(), " ")))
                    cc[0] += float(items[j - 2])
                    cc[1] += float(items[j - 1])
                    cc[2] += float(items[j])
                    n += 1
                except IndexError as err:
                    logging.error("%s get_mol_baricentre problem with %s", err, line)
    return [cc[i] / float(n) for i in range(len(cc))]


class Receptor(object):
    def __init__(self, receptor, mglpath):
        self.receptor = receptor
        self.mglpath = str(Path(mglpath).resolve())

    def topdbqt(self):
        python_env = (
            "export LD_LIBRARY_PATH=\"%s/lib\"${LD_LIBRARY_PATH:+':'$LD_LIBRARY_PATH};"
            % (self.mglpath)
        )
        python_env += "%s/bin/python2" % (self.mglpath)
        prep_rec = self.mglpath
        prep_rec += "/MGLToolsPckgs/AutoDockTools/Utilities24/"
        prep_rec += "prepare_receptor4.py"
        pdbqt = self.receptor.replace(".pdb", ".pdbqt")
        cmd = "%s %s -r '%s' -o '%s'" % (python_env, prep_rec, self.receptor, pdbqt)
        subprocess.call([cmd], shell=True)
        return str(Path(pdbqt).resolve())


class Molecule(object):
    def __init__(self, molecule, mglpath):
        self.molecule = molecule
        self.mglpath = str(Path(mglpath).resolve())

    def topdbqt(self, tran0=[]):
        """
        tran0 is the vector of centre x,y,z where to translate the molecule
        """
        obabel = "/usr/bin/obabel"
        molname = self.molecule
        if ".mol2" in molname.lower():
            molname = molname.replace(".mol2", ".pdbqt")
        else:
            molname = molname.replace(".pdb", ".pdbqt")
        cmd = "%s -p gastaiger -imol2 '%s' -opdbqt -O '%s'" % (
            obabel,
            self.molecule,
            molname,
        )
        subprocess.call([cmd], shell=True)
        # Translate to the new center
        fpdbqt = str(Path(molname).resolve())
        if len(tran0) > 0:
            mem = []
            fi = open(fpdbqt, "r")
            for line in fi:
                if "ATOM" in line:
                    v = nsplit(line.strip(), " ")
                    x = None
                    y = None
                    z = None
                    copy_line = line
                    for i in range(len(v)):
                        if i == 6:
                            x = float(v[i]) + tran0[0]
                            copy_line.replace(v[i], str(x))
                        elif i == 7:
                            y = float(v[i]) + tran0[1]
                            copy_line.replace(v[i], str(y))
                        elif i == 8:
                            z = float(v[i]) + tran0[2]
                            copy_line.replace(v[i], str(z))
                        else:
                            continue
                    mem.append(copy_line)
                else:
                    mem.append(line)
            fi.close()

            fo = open(fpdbqt, "w")
            for line in mem:
                fo.write(line)
            fo.close()
        return fpdbqt
