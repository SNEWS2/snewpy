import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.special import lpmv, gamma
import argparse

parser = argparse.ArgumentParser(description="Read neutrino energy files")
parser.add_argument("--mass",help="progenitor mass to use as input, default 10",type=str,default='10')
parser.add_argument("--neutrino",help="neutrino species to use as input, default 0",type=str,default='0')
args = parser.parse_args()

#read luminosity angle-averaged data, spectra data, and full luminosity data, respectively
lum_dat = np.loadtxt("lum_{}M_inu{}.dat".format((args.mass),(args.neutrino)))
spec_dat = np.loadtxt("lum_spec_{}M_inu{}.dat".format((args.mass),(args.neutrino)))
lum_file = h5py.File("lum_spec_{}M.h5".format(args.mass),"r")

grid_file= h5py.File("grid.h5","r")
phic = np.array(grid_file["phic"])
thetac = np.array(grid_file["thc"])
dOmega = np.array(grid_file["dOmega"])

def fact(n):
    """
    Note: math.factorial returns an int, but we want to work only with floats
    """
    return gamma(n + 1.)

def real_sph_harm(l, m, phi=phic, theta=thetac):
    """
    Computes the orthonormalized real spherical harmonics Y_lm
    """
    if m < 0:
        norm = np.sqrt((2*l + 1.)/(2*np.pi)*fact(l + m)/fact(l - m))
        return norm*np.outer(lpmv(-m, l, np.cos(theta)), np.sin(-m*phi))
    elif m == 0:
        norm = np.sqrt((2*l + 1.)/(4*np.pi))
        return norm*np.outer(lpmv(0, l, np.cos(theta)), np.ones_like(phi))
    else:
        norm = np.sqrt((2*l + 1.)/(2*np.pi)*fact(l - m)/fact(l + m))
        return norm*np.outer(lpmv(m, l, np.cos(theta)), np.cos(m*phi))

#read time and luminosity
time = lum_dat[:,0]
lum_avg = lum_dat[:,1]

#plot angle-averaged Luminosity
plt.figure()
plt.plot(time,lum_avg/100.)
plt.xlabel(r'Time after bounce [s]')
plt.ylabel(r'Luminosity [10$^{52}$ erg s$^{-1}$]')
plt.savefig(r'Lum_avg.pdf')
plt.close()

#plot nu0 luminosity spectra at a given time (index 300 here)
plt.figure()
plt.plot(spec_dat[300][1:13],spec_dat[300][13:25],'o')
plt.figtext(0.6,0.6,'Time = {} s'.format(spec_dat[300][0]))
plt.xlabel(r'Energy [MeV]')
plt.ylabel(r'$dL_{\nu_e}/d\varepsilon$ [10$^{50}$ erg s$^{-1}$ MeV$^{-1}$]')
plt.xlim([0,60])
plt.savefig(r'spec.pdf')
plt.close()

#plot mean neutrino energy, averaged over all bins and over solid angle
time = lum_file["nu{}".format(args.neutrino)]["g0"].attrs["time"]
eave = np.array(lum_file["nu{}".format(args.neutrino)]["eave"]["l=0 m=0"])/np.sqrt(4*np.pi)
plt.figure()
plt.plot(time,eave)
plt.xlabel(r'Time after bounce [s]')
plt.ylabel(r'Average Neutrino Energy [10$^{50}$ erg s$^{-1}$]')
plt.savefig(r'eave_avg.pdf')
plt.close()

#define total, angle-dependent luminosity Ltot_angle
dum_idx=1
Ltot_angle=np.zeros((len(lum_avg),len(thetac),len(phic)))
for l in range(3):
    for m in range(-l,l+1):
        print(l,m)
        Ltot_angle += real_sph_harm(l,m)*lum_dat[:,dum_idx][:,np.newaxis,np.newaxis]
        dum_idx+=1

PH, TH = np.meshgrid(phic, thetac)
cf = plt.contourf(PH, TH, Ltot_angle[0])
cb = plt.colorbar(cf)
plt.show()
