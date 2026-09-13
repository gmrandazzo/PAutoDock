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

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

from pautodock.fileutils import get_bin_path


def nsplit(s: str, delim: str | None = None) -> list[str]:
    return [x for x in s.split(delim) if x]


def extract_coordinates(line: str, ftype: str) -> list[float] | None:
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


def get_mol_baricentre(mol: str) -> list[float]:
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


def get_first_pose_baricentre(mol: str) -> list[float]:
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
    def __init__(self, receptor: str, method: str = "obabel") -> None:
        if method not in ("obabel", "mgltools"):
            raise ValueError("Unknown receptor preparation method %s" % (method))
        self.receptor = receptor
        self.method = method
        self.obabel_path = get_bin_path("obabel")
        self.mglpath = str(Path(f"{Path.home()}/.pautodock/MGLTools").resolve())

    def topdbqt(self) -> str:
        """
        Convert the receptor to pdbqt. With the default "obabel"
        method the hydrogens are added (the pdbqt writer keeps the
        polar ones, which is the AutoDock convention), Gasteiger
        charges are assigned and the water molecules are removed.
        With the "mgltools" method the MGLTools prepare_receptor4.py
        script is used instead.
        """
        if self.method == "mgltools":
            return self._topdbqt_mgltools()
        return self._topdbqt_obabel()

    def _topdbqt_obabel(self) -> str:
        obabel = Path(self.obabel_path) / "obabel"
        pdbqt = self.receptor.replace(".pdb", ".pdbqt")
        subprocess.run(
            [
                str(obabel),
                "-ipdb",
                self.receptor,
                "-opdbqt",
                "-h",
                "--partialcharge",
                "gasteiger",
                "-O",
                pdbqt,
            ],
            check=True,
        )
        self._remove_waters(pdbqt)
        return str(Path(pdbqt).resolve())

    def _topdbqt_mgltools(self) -> str:
        if not Path(self.mglpath).exists():
            msg = "MGLTools is not installed in %s. " % (self.mglpath)
            msg += "Run pautodock with --mgl ON to install it, or use "
            msg += "the default obabel receptor preparation."
            raise ValueError(msg)
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

    @staticmethod
    def _remove_waters(pdbqt: str) -> None:
        with open(pdbqt, "r", encoding="utf-8") as fi:
            lines = fi.readlines()
        with open(pdbqt, "w", encoding="utf-8") as fo:
            for line in lines:
                if ("ATOM" in line or "HETATM" in line) and (
                    "HOH" in line[17:20] or "WAT" in line[17:20]
                ):
                    continue
                fo.write(line)


class Molecule(object):
    def __init__(self, molecule: str) -> None:
        self.molecule = molecule
        self.obabel_path = get_bin_path("obabel")

    def topdbqt(
        self, center: list[float] | None = None, ph: float | None = None
    ) -> str:
        """
        center is the x,y,z point where the molecule baricentre
        will be translated to.

        ph, when set, protonates the molecule at the given pH
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
        if ph is not None:
            # Protonate the ligand at the given pH
            cmd += " -p %s" % (ph)
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
