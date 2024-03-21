"""molop.py
Provides the basic operation for molecular files.

"""

#!/usr/bin/env python3

import logging
from pathlib import Path
import subprocess


def nsplit(s, delim=None):
    return [x for x in s.split(delim) if x]

def getMolBaricentre(mol: str) -> tuple:
    """
    Get molecular baricentre from a molecule
    """
    cc = [0., 0., 0.]
    n = 0
    if ".pdbqt" in mol:
        j=-5
    else:
        j=-4
    with open(mol, 'r', encoding='utf-8') as f:
        for line in f:
            if (("ATOM" in line or "HETATM" in line) and not 'REMARK' in line):
                try:
                    # Get the results from back
                    items = list(filter(None, str.split(line.strip(), " ")))
                    cc[0] += float(items[j-2])
                    cc[1] += float(items[j-1])
                    cc[2] += float(items[j])
                    n += 1
                except IndexError as err:
                    logging.error('%s getMolBaricentre problem with %s', err, line)
    return [cc[i]/float(n) for i in range(len(cc))]

class Receptor(object):
    def __init__(self, receptor, mglpath):
        self.receptor = receptor
        self.mglpath = str(Path(mglpath).resolve())

    def topdbqt(self):
        python_env = 'export LD_LIBRARY_PATH="%s/lib"${LD_LIBRARY_PATH:+\':\'$LD_LIBRARY_PATH};' % (self.mglpath)
        python_env += "%s/bin/python2" % (self.mglpath)
        prep_rec = self.mglpath
        prep_rec += "/MGLToolsPckgs/AutoDockTools/Utilities24/"
        prep_rec += "prepare_receptor4.py"
        pdbqt = self.receptor.replace(".pdb", ".pdbqt")
        cmd = ("%s %s -r '%s' -o '%s'" %
               (python_env, prep_rec, self.receptor, pdbqt))
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
        obabel="/usr/bin/obabel"
        molname = self.molecule
        if ".mol2" in molname.lower():
            molname = molname.replace(".mol2", ".pdbqt")
        else:
            molname = molname.replace(".pdb", ".pdbqt")
        cmd = ("%s -p gastaiger -imol2 '%s' -opdbqt -O '%s'" %
               (obabel, self.molecule, molname))
        subprocess.call([cmd], shell=True)
        #Translate to the new center
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
                            x = (float(v[i])+tran0[0])
                            copy_line.replace(v[i], str(x))
                        elif i == 7:
                            y = (float(v[i])+tran0[1])
                            copy_line.replace(v[i], str(y))
                        elif i == 8:
                            z = (float(v[i])+tran0[2])
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
    
    def topdbqt_old(self, tran0=[]):
        python_env = 'export LD_LIBRARY_PATH="%s/lib"${LD_LIBRARY_PATH:+\':\'$LD_LIBRARY_PATH};' % (self.mglpath)
        python_env += "%s/bin/python2" % (self.mglpath)
        prep_lig = self.mglpath
        prep_lig += "/MGLToolsPckgs/AutoDockTools/Utilities24/"
        prep_lig += "prepare_ligand4.py"
        molname = self.molecule
        if ".mol2" in molname.lower():
            molname = molname.replace(".mol2", ".pdbqt")
        else:
            molname = molname.replace(".pdb", ".pdbqt")
        cmd = ("%s %s -l '%s' -o '%s'" %
               (python_env, prep_lig, self.molecule, molname))
        subprocess.call([cmd], shell=True)
        #Translate to the new center
        fpdbqt = str(Path(molname).resolve())
        if len(tran0) > 0:
            mem = []
            fi = open(fpdbqt, "r")
            for line in fi:
                if "ATOM" in line:
                    v = nsplit(line.strip(), " ")
                    for i in range(len(v)):
                        if i == 5:
                            v[i] = (float(v[i])+tran0[0])
                        elif i == 6:
                            v[i] = (float(v[i])+tran0[1])
                        elif i == 7:
                            v[i] = (float(v[i])+tran0[2])
                        else:
                            continue
                    mem.append("%s%7s  %-4s%s %s%16.3f%8.3f%8.3f%6.2f%6.2f%10.3f %s\n" %
                             (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], float(v[8]), float(v[9]), float(v[10]), v[11]))
                else:
                    mem.append(line)
            fi.close()

            fo = open(fpdbqt, "w")
            for line in mem:
                fo.write(line)
            fo.close()
        return fpdbqt

