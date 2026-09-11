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
    # Coordinates live in the fixed-width fields 31-38, 39-46, 47-54
    # (1-based PDB specification).
    if ftype == "pdb":
        if "ATOM" in line or "HETATM" in line:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            return [x, y, z]
    elif ftype == "pdbqt":
        if "ATOM" in line:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            return [x, y, z]
    return None


def get_mol_baricentre(mol: str) -> list:
    """
    Get the geometric centre (unweighted mean of atomic coordinates)
    of a molecule.
    """
    cc = [0.0, 0.0, 0.0]
    n = 0
    if mol.endswith(".pdbqt"):
        ftype = "pdbqt"
    elif mol.endswith(".pdb"):
        ftype = "pdb"
    else:
        raise ValueError(
            "Molecule format not supported %s. Supported formats: pdb or pdbqt"
            % (mol)
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
    if n == 0:
        raise ValueError("No atoms found in %s" % (mol))
    return [cc[i] / float(n) for i in range(len(cc))]


def get_first_pose_baricentre(mol: str) -> list:
    """
    Get the geometric centre of the first pose (up to ENDMDL) of a
    multimodel pdbqt file, e.g. vina docking poses.
    """
    cc = [0.0, 0.0, 0.0]
    n = 0
    if not mol.endswith(".pdbqt"):
        raise ValueError(
            "Molecule format not supported %s. Supported format: pdbqt" % (mol)
        )
    with open(mol, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if ("ATOM" in line or "HETATM" in line) and "REMARK" not in line:
                ex_cc = extract_coordinates(line.strip(), "pdbqt")
                if ex_cc:
                    for i, val in enumerate(ex_cc):
                        cc[i] += val
                    n += 1
    if n == 0:
        raise ValueError("No atoms found in %s" % (mol))
    return [cc[i] / float(n) for i in range(len(cc))]


def read_active_torsions(pdbqt: str) -> int:
    """
    Read the number of active torsions written by Open Babel in the
    REMARK of a pdbqt file. This is the value AutoDock expects for
    the torsdof parameter.
    """
    with open(pdbqt, "r", encoding="utf-8") as f:
        for line in f:
            if "active torsions:" in line:
                parts = line.replace(":", " ").split()
                return int(parts[parts.index("active") - 1])
    raise ValueError("Active torsions remark not found in %s" % (pdbqt))


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

    def topdbqt(self, center=None):
        """
        center is the x,y,z point where the molecule baricentre
        will be translated to.
        """
        obabel = f"{self.obabel_path}/obabel"
        molname = self.molecule
        if ".mol2" in molname.lower():
            molname = molname.replace(".mol2", ".pdbqt")
        else:
            molname = molname.replace(".pdb", ".pdbqt")
        cmd = "%s --partialcharge gasteiger -imol2 '%s' -opdbqt -O '%s'" % (
            obabel,
            self.molecule,
            molname,
        )
        subprocess.call([cmd], shell=True)
        # Translate the molecule so that its baricentre is at center
        fpdbqt = str(Path(molname).resolve())
        if center:
            cc_mol = get_mol_baricentre(fpdbqt)
            tran0 = [center[i] - cc_mol[i] for i in range(3)]
            mem = []
            fi = open(fpdbqt, "r", encoding="utf-8")
            for line in fi:
                if "ATOM" in line:
                    ex_cc = extract_coordinates(line.strip(), "pdbqt")
                    if ex_cc:
                        x = ex_cc[0] + tran0[0]
                        y = ex_cc[1] + tran0[1]
                        z = ex_cc[2] + tran0[2]
                        copy_line = "%s%8.3f%8.3f%8.3f%s" % (
                            line[:30],
                            x,
                            y,
                            z,
                            line[54:],
                        )
                        mem.append(copy_line)
                    else:
                        msg = "Molecule.topdbqt Error!\n"
                        msg += " X Y Z coordinates not found in "
                        msg += f"line {line.strip()}"
                        raise ValueError(msg)
                else:
                    mem.append(line)
            fi.close()

            fo = open(fpdbqt, "w", encoding="utf-8")
            for line in mem:
                fo.write(line)
            fo.close()
        return fpdbqt
