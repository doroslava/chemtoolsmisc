
import MDAnalysis as mda
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.rms import rmsd
#reformat files
#Add chain information, remove water
path="/mdspace/dsribar/TLR8/MD/DESMOND/ssRNA/noions/Analysis/"
os.mkdir(path + "pycontact/")
for i in range (0,3):
        structureFile=path+"test"+str(i)+".psf"
        trajectoryFile=path+"test"+str(i)+".dcd"

        u=mda.Universe(structureFile,trajectoryFile)

        u.atoms[0:12155].segids = "A"
        u.atoms[24382:24389].segids = "A"
        u.atoms[24395:24398].segids = "A"
        u.atoms[12155:24190].segids = "B"
        u.atoms[24389:24395].segids = "B"
        u.atoms[24398].segid = "B"
#       u.atoms[24229:24295].resnames = "LIG"

        A = u.selectAtoms("not resname SPC")

        A.write(path+"pycontact/pycontact"+str(i)+".pdb")

        with mda.Writer(path+"pycontact/pycontact"+str(i)+".dcd", A.n_atoms) as W:
            for ts in u.trajectory:
                W.write(A)

#RMSD plots
path="/mdspace/dsribar/TLR8/MD/DESMOND/CL097/Holo/Analysis/"
for i in range (0,3):
        structureFile=path+"dyno"+str(i)+".pdb"
        trajectoryFile=path+"dyno"+str(i)+".dcd"

        u = mda.Universe(structureFile, trajectoryFile)
        domains = {
        '1': u.select_atoms("resname LIG"),
        }
        colors = {'1': 'black'}

        u.trajectory[0]  # rewind trajectory
        xref0 = dict((name, g.positions - g.center_of_mass()) for name, g in domains.items())

        nframes = len(u.trajectory)
        results = dict((name, np.zeros((nframes, 2), dtype=np.float64)) for name in domains)

        for iframe, ts in enumerate(u.trajectory):
                for name, g in domains.items():
                        results[name][iframe, :] = (u.trajectory.time,
                                            rmsd(g.positions, xref0[name],
                                                 center=True, superposition=True))




        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        for name in ("1"):
                data = results[name]
                ax.plot(data[:,0], data[:,1], linestyle="-", color=colors[name], lw=2, label=name)
        ax.legend(loc="best")
        ax.set_xlabel(r"time  $t$ (ps)")
        ax.set_ylabel(r"C$_\alpha$ RMSD from $t=0$, $\rho_{\mathrm{C}_\alpha}$ ($\AA$)")

        for ext in ('svg', 'pdf', 'png'):
                fig.savefig("AdK_domain_rigidity.{0}".format(ext))