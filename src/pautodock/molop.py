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

from pautodock.fileutils import get_bin_path


def nsplit(s, delim=None):
    return [x for x in s.split(delim) if x]


def extract_coordinates(line, ftype):
    if ftype == "pdb":
        if "HETATM" in line:
            x = float(line[29:39].strip())
            y = float(line[39:47].strip())
            z = float(line[47:55].strip())
            return [float(x), float(y), float(z)]
    elif ftype == "pdbqt":
        if "ATOM" in line:
            x = float(line[31:39].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            return [float(x), float(y), float(z)]
    else:
        return None


def get_mol_baricentre(mol: str) -> tuple:
    """
    Get molecular baricentre from a molecule
    """
    cc = [0.0, 0.0, 0.0]
    n = 0
    if mol.endswith(".pdbqt"):
        ftype = "pdbqt"
    elif mol.endswith(".pdb"):
        ftype = "pdb"
    else:
        raise ValueError(
            "Molecual format not supported {mol}. Supported formats: pdb or pdbqt"
        )

    with open(mol, "r", encoding="utf-8") as f:
        for line in f:
            if ("ATOM" in line or "HETATM" in line) and "REMARK" not in line:
                try:
                    ex_cc = extract_coordinates(line.strip(), ftype)
                    if ex_cc:
                        for i, val in enumerate(ex_cc):
                            cc[i] += val
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
        self.obabel_path = get_bin_path("obabel")

    def topdbqt(self, tran0=[]):
        """
        tran0 is the vector of centre x,y,z where to translate the molecule
        """
        obabel = f"{self.obabel_path}/obabel"
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
                    ex_cc = extract_coordinates(line.strip(), "pdbqt")
                    if ex_cc:
                        copy_line = line
                        x = float(ex_cc[0]) + tran0[0]
                        copy_line.replace(str(ex_cc[0]), str(x))
                        y = float(ex_cc[1]) + tran0[0]
                        copy_line.replace(str(ex_cc[1]), str(y))
                        z = float(ex_cc[2]) + tran0[0]
                        copy_line.replace(str(ex_cc[2]), str(z))
                        mem.append(copy_line)
                    else:
                        msg = "Molecule.topdbqt Error!\n"
                        msg += " X Y Z coordinates not found in "
                        msg += f"line {line.strip()}"
                        raise ValueError(msg)
                else:
                    mem.append(line)
            fi.close()

            fo = open(fpdbqt, "w")
            for line in mem:
                fo.write(line)
            fo.close()
        return fpdbqt
