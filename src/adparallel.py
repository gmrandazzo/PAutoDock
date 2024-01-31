"""adparallel.py

Provides the basic code to run autodock in parallel on a machine.

"""

#!/usr/bin/env python3

import molop
import multimol2op
import os
import tempfile
import shutil
from pathlib import Path
import multiprocessing


class ADParallel(object):
    def __init__(
        self,
        atdpath: str,
        mglpath: str,
        receptor: str,
        ligand: str,
        db: str,
        wpath: str):
        if len(atdpath) == 0:
            self.atdpath = "/usr/bin/"
        else:
            self.atdpath = atdpath
        self.mglpath = mglpath
        self.receptor = receptor
        self.ligand = ligand
        self.db = db
        self.wpath = wpath
        self.results = []
        self.cx = 0.
        self.cy = 0.
        self.cz = 0.
        self.gsize_x = 30
        self.gsize_y = 30
        self.gsize_z = 30
        self.speed = "slow"
        self.atd = True
        self.vina = True
        self.smina_variant = True

    def readAtomTypes(self, rec_mol):
        # Read the atom types
        atlst = []
        f = open(rec_mol, "r")
        for line in f:
            if "ATOM" in line:
                at = line.strip().split(" ")[-1]
                if at in atlst:
                    continue
                else:
                    if len(at) > 0:
                        atlst.append(at)
                    else:
                        continue
        f.close()
        atlst = list(set(atlst))
        atypes_str = "%s" % atlst[0]
        for i in range(1, len(atlst)):
            atypes_str += " %s" % atlst[i]
        return atypes_str, atlst

    def WriteParamFiles(self, path, rec_pdbqt, mol_pdbqt, cc):
        path_ = Path(path).absolute()
        rat_str, rat_lst = self.readAtomTypes(rec_pdbqt)
        lat_str, lat_lst = self.readAtomTypes(path + "/" + mol_pdbqt)
        # write the GPF
        f = open(path + "/grid.gpf", "w")
        f.write("npts %d %d %d # num.grid points in xyz\n" %
                (self.gsize_x, self.gsize_y, self.gsize_z))
        f.write("gridfld %s/receptor_model.maps.fld     # grid_data_file\n" %
                (path_))
        f.write("spacing 0.33 # spacing(A)\n")
        f.write("receptor_types %s # receptor atom types\n" %
                (rat_str))
        f.write("ligand_types %s # ligand atom types\n" %
                (lat_str))
        f.write("receptor %s    # macromolecule\n" % (rec_pdbqt))
        f.write("gridcenter %.3f %.3f %.3f       # xyz-coordinates or auto\n" %
                (cc[0], cc[1], cc[2]))
        f.write("smooth 0.5  # store minimum energy w/in rad(A)\n")
        for at in lat_lst:
            f.write("map %s/receptor_model.%s.map # atom-specific aff. map\n" %
                    (path_, at))
        f.write("elecmap %s/receptor_model.e.map # El. Pot potential map\n" %
                (path_))
        f.write("dsolvmap %s/receptor_model.d.map # desolv. potential map\n" %
                (path_))
        f.write("dielectric -0.1465 # <0, AD4 dist-dep.diel;>0, constant\n")
        f.close()

        # Write autodock file
        f = open(path + "/ind.dpf", "w")
        f.write("autodock_parameter_version 4.2")
        f.write("# used by autodock to validate parameter set\n")
        f.write("outlev 1 # diagnostic output level\n")
        f.write("intelec # calculate internal electrostatics\n")
        f.write("seed time pid # seeds for random generator\n")
        f.write("ligand_types %s # atoms types in ligand\n" % (lat_str))
        f.write("fld %s/receptor_model.maps.fld      # grid_data_file\n" %
                (path_))
        for at in lat_lst:
            f.write("map %s/receptor_model.%s.map # atom-specific aff. map\n" %
                    (path_, at))
        f.write("elecmap %s/receptor_model.e.map    # electrostatics map\n" %
                (path_))
        f.write("desolvmap %s/receptor_model.d.map   # desolvation map\n" %
                (path_))
        f.write("move %s/%s                    # small molecule\n" %
                (path_, mol_pdbqt))
        f.write("about %.4f %.4f %.4f       # small molecule center\n" %
                (cc[0], cc[1], cc[2]))
        f.write("tran0 %.4f %.4f %.4f " % (cc[0], cc[1], cc[2]))
        f.write("# initial coordinates/A or random\n")
        f.write("quaternion0 random ")
        f.write("# initial orientation\n")
        f.write("dihe0 random ")
        f.write("# initial dihedrals (relative) or random\n")
        f.write("torsdof 5 ")
        f.write("# torsional degrees of freedom\n")
        f.write("rmstol 2.0  # cluster_tolerance/A\n")
        f.write("extnrg 1000.0  # external grid energy\n")
        f.write("e0max 0.0 10000 ")
        f.write("# max initial energy; max number of retries\n")
        f.write("ga_pop_size 150 ")
        f.write("# number of individuals in population\n")
        if self.speed == "slow":
            f.write("ga_num_evals 25000000 ")
        elif self.speed == "medium":
            f.write("ga_num_evals 2500000 ")
        else:
            f.write("ga_num_evals 250000 ")
        f.write("               # maximum number of energy evaluations\n")
        f.write("ga_num_generations 27000 ")
        f.write("# maximum number of generations\n")
        f.write("ga_elitism 1 ")
        f.write("# number of top individuals to survive to next generation\n")
        f.write("ga_mutation_rate 0.02 # rate of gene mutation\n")
        f.write("ga_crossover_rate 0.8 # rate of crossover\n")
        f.write("ga_window_size 10 # \n")
        f.write("ga_cauchy_alpha 0.0 ")
        f.write("# Alpha parameter of Cauchy distribution\n")
        f.write("ga_cauchy_beta 1.0 ")
        f.write("# Beta parameter Cauchy distribution\n")
        f.write("set_ga")
        f.write(" # set the above parameters for GA or LGA\n")
        f.write("sw_max_its 300")
        f.write("# iterations of Solis & Wets local search\n")
        f.write("sw_max_succ 4 ")
        f.write(" # consecutive successes before changing rho\n")
        f.write("sw_max_fail 4 ")
        f.write(" # consecutive failures before changing rho\n")
        f.write("sw_rho 1.0 ")
        f.write("# size of local search space to sample\n")
        f.write("sw_lb_rho 0.01                       # lower bound on rho\n")
        f.write("ls_search_freq 0.06 ")
        f.write("# probability of performing local search on individual\n")
        f.write("set_psw1 ")
        f.write("# set the above pseudo-Solis & Wets parameters\n")
        f.write("unbound_model bound ")
        f.write("# state of unbound ligand\n")
        f.write("ga_run 10 ")
        f.write("# do this many hybrid GA-LS runs\n")
        f.write("analysis ")
        f.write("# perform a ranked cluster analysis\n")
        f.close()
        grid_path = Path(path + "/grid.gpf").absolute()
        ind_path = Path(path + "/ind.dpf").absolute()
        return grid_path, ind_path

    def WriteVinaParams(self, path, cc, ss):
        f = open(path + "/vina_conf.txt", "w")
        f.write("center_x = %.4f\n" % (cc[0]))
        f.write("center_y = %.4f\n" % (cc[1]))
        f.write("center_z = %.4f\n" % (cc[2]))
        f.write("size_x = %d\n" % (ss[0]))
        f.write("size_y = %d\n" % (ss[1]))
        f.write("size_z = %d\n" % (ss[2]))
        f.write("num_modes = 9\n")
        f.close()
        return Path(path + "/vina_conf.txt").absolute()

    def RunAutoGrid(self, cmd):
        atg_path = str(Path("%s/autogrid4" % (self.atdpath)).absolute())
        return os.system("%s %s" % (atg_path, cmd))

    def RunAutoDock(self, cmd):
        atd_path = str(Path("%s/autodock4" % (self.atdpath)).absolute())
        return os.system("%s %s" % (atd_path, cmd))

    def RunVina(self, cmd):
        bin_path = None
        if self.smina_variant is True:
            bin_path = f'{Path(__file__).parent.parent}/third_party/smina.static'
        else:
            bin_path = str(Path("%s/vina" % (self.atdpath)).absolute())
        return os.system("%s %s" % (bin_path, cmd))

    def ReadOutput(self, ofile):
        r = []
        header = []
        benergy = []
        c_rmsd = []
        r_rmsd = []
        f = open(str(Path(ofile).absolute()), "r")
        for line in f:
            if "Partition function, Q =" in line:
                header.append("Part. Func.")
                r.append(molop.nsplit(line.strip(), " ")[4])
            elif "Free energy,        A ~" in line:
                r.append(molop.nsplit(line.strip(), " ")[4])
                header.append("Free Energy")
            elif "Internal energy,    U =" in line:
                header.append("Internal Energy")
                r.append(molop.nsplit(line.strip(), " ")[4])
            elif "Entropy,            S =" in line:
                header.append("Entropy")
                r.append(molop.nsplit(line.strip(), " ")[3])
            elif "RANKING" in line:
                v = molop.nsplit(line.strip(), " ")
                benergy.append(float(v[3]))
                c_rmsd.append(float(v[4]))
                r_rmsd.append(float(v[5]))
        f.close()
        header.append("Binding Energy Average")
        r.append(round(sum(benergy) / float(len(benergy)),3))
        header.append("Cluster RMSD Average")
        r.append(round(sum(c_rmsd) / float(len(c_rmsd)),3))
        header.append("Ref. RMSD Average")
        r.append(round(sum(r_rmsd) / float(len(r_rmsd)),3))
        return header, r

    def ReadVinaOutput(self, ofile):
        benergy = []
        f = open(str(Path(ofile).absolute()), "r")
        getres = False
        for line in f:
            if getres is True:
                if "Writing output ... done." in line:
                    getres = False
                elif 'Refine time' in line:
                    getres = False
                else:
                    v = molop.nsplit(line.strip(), " ")
                    benergy.append(float(v[1]))
            else:
                if "-----+------------+----------+----------" in line:
                    getres = True
                else:
                    continue
        f.close()
        if len(benergy) > 0:
            return round(sum(benergy) / float(len(benergy)), 3), min(benergy), max(benergy)
        else:
            return 9999., 9999., 9999.

    def GenVSOutput(self, vinalogout, dpfout, mnames, otab):
        # Collect the vina results
        vbind = []
        for vout in vinalogout:
            avg_b, min_b, max_b = self.ReadVinaOutput(vout)
            vbind.append([avg_b, min_b, max_b])
        # Collect results
        fo = open(otab, "w")
        firstline = True
        for i in range(len(mnames)):
            h = []
            r = []
            if len(dpfout) > 0:
                if Path(dpfout[i]).is_file():
                    h, r = self.ReadOutput(dpfout[i])
            if firstline is True:
                firstline = False
                fo.write("Molname;")
                for j in range(len(h)):
                    fo.write("%s;" % (h[j]))
                fo.write("Avg. vina Binding Energy;")
                fo.write("Min vina Binding Energy;Max vina Binding Energy\n")
            fo.write("%s;" % (mnames[i]))
            for j in range(len(r)):
                fo.write("%s;" % (r[j]))
            fo.write("%f;%f;%f\n" % (vbind[i][0], vbind[i][1], vbind[i][2]))
        fo.close()

    def VS(self, otab):
        # Prepare the receptor
        rec = molop.Receptor(self.receptor, self.mglpath)
        rec_pdbqt = rec.topdbqt()
        cc = None
        if self.ligand is None:
            cc = [self.cx, self.cy, self.cz]
        else:
            # Prepare the ligand
            lig = molop.Molecule(self.ligand, self.mglpath)
            cc = lig.getMolCentre()
        # Prepare the database
        # split multi mol2
        tmppath = tempfile.mkdtemp()
        mol2lst = multimol2op.SplitMol2(self.db, tmppath)
        agcmdlst = []
        adcmdlst = []
        vinacmdlst = []
        vinalogout = []
        dpfout = []
        mnames = []
        # Create a directory with the name of the mol2 molecule
        # and copy the receptor and itself
        for mol2 in mol2lst:
            molname_ext = str(Path(mol2).resolve().name)
            molname = molname_ext.replace(".mol2", "")
            mnames.append(molname)
            mpath = str(Path(self.wpath + "/" + molname).absolute())
            if Path(mpath).exists() is False:
                os.makedirs(mpath)
            if Path(mpath + "/" + molname_ext).exists() is False:
                shutil.move(str(Path(mol2).resolve()), mpath)

            mol = molop.Molecule(str(Path(mpath + "/" + molname_ext).absolute()), self.mglpath)
            mol_pdbqt = mol.topdbqt(cc)
            mol_pdbqt_name = str(Path(mol_pdbqt).resolve().name)
            if self.atd is True:
                gpf_path, dpf_path = self.WriteParamFiles(mpath,
                                                rec_pdbqt,
                                                mol_pdbqt_name,
                                                cc)
                ag = "-p \"%s\" -l \"%s.glg\"" % (str(gpf_path), str(gpf_path).replace(".gpf", ""))
                agcmdlst.append(ag)
                dpfout.append(str(dpf_path).replace(".dpf", ".dlg"))
                ad = "-p \"%s\" -l \"%s\"" % (str(dpf_path), dpfout[-1])
                adcmdlst.append(ad)

            if self.vina is True:
                ss = [self.gsize_x, self.gsize_y, self.gsize_z]
                vconf_path = self.WriteVinaParams(mpath, cc, ss)
                vc = "--config \"%s\" --receptor \"%s\" --ligand \"%s\"" % (vconf_path, rec_pdbqt, mol_pdbqt)
                vinalogout.append("%s/vina_log.txt" % (mpath))
                vc += " --out \"%s/dock_confs_%s.pdbqt\" >> \"%s\"" % (mpath, molname, vinalogout[-1])
                vinacmdlst.append(vc)

        shutil.rmtree(tmppath)
        # Run AutoGrid
        ncpu = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(ncpu)
        pool.map(self.RunAutoGrid, agcmdlst)

        # RunAutodock
        pool = multiprocessing.Pool(ncpu)
        pool.map(self.RunAutoDock, adcmdlst)

        # RunVina or Smina
        pool = multiprocessing.Pool(ncpu)
        pool.map(self.RunVina, vinacmdlst)

        # Write the output table
        self.GenVSOutput(vinalogout, dpfout, mnames, otab)

        return 1
