#%%

import os
import subprocess
import glob
import datetime
import numpy as np
import pandas as pd

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker
get_ipython().run_line_magic("matplotlib","qt")
matplotlib.rcParams['mathtext.default'] = 'rm'

MPI = 0.14                          # pion mass = 0.14 GeV
MSIGMA_MIN, MSIGMA_MAX = 2*MPI, 2   # sigma mass, spectral function considered up to 2 GeV
NMASSSAMPLES = 300                  # 300 masses sampled
# parentdir = "data/examplescenario"
# parentdir = "data/AuAu200_cc1"
parentdir = "data/XeXe544_cc1"
# parentdir = "data/PbPb276_cc1"

#%%
####
####    PARAMETERS FOR PRIMARY (SIGMA) SPECTRUM COMPUTATION
####

dccinitdata = glob.glob(parentdir+"/initdata.csv")[0]
freezeoutdata = glob.glob(parentdir+"/freezeout.csv")[0]
masses = np.linspace(MSIGMA_MIN, MSIGMA_MAX,NMASSSAMPLES+2)[1:-1]
qTmax = 4
NqT = 400
epsrel = 1e-5
iterations=10000
primaryspecfolder = "primaryspec"

#%%
####
####    COMPUTE PRIMARY SPECTRA
####

# prevent accidental execution of notebook cell
confirm = input("start computation? (y/n):").strip().lower()
if confirm != 'y':
    raise RuntimeError("aborted")

for (i,m) in enumerate(masses):
    foldername = primaryspecfolder+"/spec_{:%Y%m%d_%H%M%S}_".format(datetime.datetime.now())+str(i).zfill(int(np.ceil(np.log10(len(masses)))))
    print("save to", foldername)
    result = subprocess.run(args=[
        "./bin/specV2",
        "--m=%f"%(masses[i]),
        "--pTmax=%f"%(qTmax),
        "--NpT=%d"%(NqT),
        "--epsabs=0",
        "--epsrel=%f"%(epsrel),
        "--iter=%d"%(iterations),
        "--parentdir=%s"%(parentdir),
        "--initpath=%s"%(dccinitdata),
        "--foldername=%s"%(foldername),
        "--freezeoutpath=%s"%(freezeoutdata)
    ])
    print(result)

#%%
####
####    VISUALIZE PRIMARY SPECTRA
####

TICKLABELSIZE=20
FIGSIZE = (7,7)
AXISLABELSIZE = 20
LINEWIDTH = 2

CMAP = LinearSegmentedColormap.from_list("custom", ["blue","red"])
CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]
CMAP_LABELSIZE = 15
CMAP_TICKSIZE = 15
LC_LABEL = r"$m\ [GeV]$"

SPEC_XLABEL = r"$q_T\ [GeV]$"
SPEC_YLABEL = r"$(2\pi q_T)^{-1}dN_{coherent}/(dq_Td\eta_q)\ [GeV^{-2}]$"

fig, ax = plt.subplots(figsize=(7,7))
lines = []

specfiles = sorted(glob.glob(parentdir+"/"+primaryspecfolder+"/*/*spectr.txt"))

masses = np.zeros(len(specfiles))
for (i,file) in enumerate(specfiles):
    with open(file) as f:
        mylines = f.readlines()
        masses[i] = float(mylines[3].replace("# particle mass:\t",""))

for (i, file) in enumerate(specfiles):
    df = pd.read_csv(file,comment="#")

    pTs, spec = df.to_numpy().T
    lines.append(np.column_stack((pTs, spec)))
    
linecol = LineCollection(lines,array=masses,cmap=CMAP,lw=LINEWIDTH)
ax.add_collection(linecol)
ax.set_yscale("log")
ax.set_xlim(0,qTmax)
ax.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)
ax.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)

ax.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.grid(False)

cax = ax.inset_axes(CMAP_LBWH)
cbar = fig.colorbar(linecol, cax=cax)
cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)

fig.tight_layout()
fig.show()

#%%
####
####    PARAMETERS FOR DECAY (PION) SPECTRUM COMPUTATION
####

decayspecfolder = "decayspec"
pTmax = 1
NpT = 50
epsrel=1e-5
iterations = 100000
primespecfiles = sorted(glob.glob(parentdir+"/"+primaryspecfolder+"/*/*spectr.txt"))
B=1
Q=1

#%%
####
####    COMPUTE DECAY SPECTRA
####

# prevent accidental execution of notebook cell
confirm = input("start computation? (y/n):").strip().lower()
if confirm != 'y':
    raise RuntimeError("aborted")

masses = np.zeros(len(primespecfiles))
for (i,file) in enumerate(primespecfiles):
    with open(file) as f:
        mylines = f.readlines()
        masses[i] = float(mylines[3].replace("# particle mass:\t",""))

for (i,MSIGMA) in enumerate(masses):
    foldername = decayspecfolder+"/decay_{:%Y%m%d_%H%M%S}_".format(datetime.datetime.now())+str(i).zfill(int(np.ceil(np.log10(len(masses)))))
    print("save to", foldername)
    result = subprocess.run(args=[
        "./bin/decayV2",
        "--ma=%f"%(MSIGMA),
        "--mb=%f"%(MPI),
        "--mc=%f"%(MPI),
        "--pTmax=%f"%(pTmax),
        "--NpT=%d"%(NpT),
        "--epsabs=0",
        "--epsrel=%f"%(epsrel),
        "--iter=%d"%(iterations),
        "--primespecpath=%s"%(primespecfiles[i]),
        "--parentdir=%s"%(parentdir),
        "--foldername=%s"%(foldername),
        "--B=%f"%(B),
        "--Q=%f"%(Q)
    ])
    print(result)

#%%
####
####    VISUALIZE DECAY SPECTRA
####

SCALE = 0.25

### \begin{equation}
###     \frac{1}{2\pi p_{\text{T}}}\frac{\dt N}{\dt p_{\text{T}}\dt\eta_p}\Big\vert_{\pi^+}\sim B\,\vert\underbrace{\mathcal{A}_\sigma}_{\text{condensate amplitude}}\vert^2=\underbrace{\tilde{B}}_{\to 1}\ \vert\underbrace{\tilde{\mathcal{A}}_\sigma}_{\to\SI{0.1}{\GeV}}\vert^2\cdot\text{scale}\quad\implies\quad\mathcal{A}_\sigma=\tilde{\mathcal{A}}_\sigma\sqrt{\text{scale}/B}
### \end{equation}

### B... branching ratio, A... sigma condensate amplitude
### in computation we instead use for simplicity
###     B~ = 1
###     A~ = 0.1 GeV
### and adjust the scale in the plotting script.
###
### ->  spec \propto B * |A|^2 = B~ * |A~|^2 * scale        =>      A = A~ * sqrt(scale/B)
###

TICKLABELSIZE=20
FIGSIZE = (7,7)
AXISLABELSIZE = 20
LINEWIDTH = 2
MARKERSIZE = 5

CMAP = LinearSegmentedColormap.from_list("custom", ["blue","red"])
CMAP_LBWH = [0.025, 0.025, 0.05, 0.45]
CMAP_LABELSIZE = 15
CMAP_TICKSIZE = 15
LC_LABEL = r"$m\ [GeV]$"

SPEC_XLABEL = r"$p_T\ [GeV]$"
SPEC_YLABEL = r"$(2\pi p_T)^{-1}dN_{\pi^+}/(dp_Td\eta_p)\ [GeV^{-2}]$"
LEGEND= r"$\sigma_{DCC}\to\pi\,\pi$"

fig, ax = plt.subplots(figsize=(7,7))
fig_full, ax_full = plt.subplots(figsize=(7,7))

lines = []
pT_full, spec_full = np.zeros(shape=(2,1))

decayfiles = sorted(glob.glob(parentdir+"/"+decayspecfolder+"/*/*decayspec.txt"))

### COMPUTE SPECTRAL FUNCTION AND WEIGHTS
Mpole = 0.4
Gpole = 2* 0.2

msigma = np.sqrt(1/4 * (16 * MPI**2 + 
                        np.sqrt(16 * Gpole**2 * Mpole**2 + 
                                (-16 * MPI**2 - Gpole**2 + 4*Mpole**2)**2)))
Gam = np.sqrt(1/2 * (16 * MPI**2 + Gpole**2 - 
     4*  Mpole**2 + np.sqrt(16*Gpole**2 * Mpole**2 + (-16 * MPI**2 - Gpole**2 + 4*Mpole**2)**2)))

def Delta(s):
    return 1/(s-msigma**2+1j*Gam*np.sqrt(s-(2*MPI)**2))

def S(k):
    return -1/np.pi*np.imag(Delta(k**2))

masses = np.zeros(len(decayfiles))
for (i,file) in enumerate(decayfiles):
    with open(file) as f:
        mylines = f.readlines()
        masses[i] = float(mylines[3].replace("# ma:\t",""))

weights = 2*masses*S(masses)*np.ptp(masses)/len(masses)
weights /= np.sum(weights)
###

for (i, file) in enumerate(decayfiles):
    df = pd.read_csv(file,comment="#")

    pTs, spec = df.to_numpy().T
    spec = SCALE * spec
    lines.append(np.column_stack((pTs, spec)))

    ### avoid gaps due to failed integral evaluations that lead to inf
    exclude = np.isinf(spec)
    idcs = np.where(1-exclude)
    spec = np.exp(np.poly1d(np.polyfit(pTs[idcs],np.log(spec[idcs]),10))(pTs))


    pT_full = pTs
    spec_full = spec_full + weights[i] * spec
    
### INDIVIDUAL DECAYSPEC PLOT
linecol = LineCollection(lines,array=masses,cmap=CMAP,lw=LINEWIDTH)
ax.add_collection(linecol)
ax.set_yscale("log")
ax.set_xlim(0,pTmax)
ax.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)
ax.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)

ax.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.grid(False)

cax = ax.inset_axes(CMAP_LBWH)
cbar = fig.colorbar(linecol, cax=cax)
cbar.set_label(LC_LABEL, fontsize=CMAP_LABELSIZE)
cbar.ax.tick_params(labelsize=CMAP_TICKSIZE)

fig.tight_layout()
fig.show()
###

### FULLSPEC PLOT
# PROCESS EXPERIMENTAL AND FLUIDUM DATA
df_exp = pd.read_csv(parentdir+"/experimentaldata.csv",comment="#")
pTs_exp, spec_exp, spec_exp_err = df_exp.to_numpy().T

df_fluid = pd.read_csv(parentdir+"/fluidumdata.csv",comment="#")
pTs_fluid, spec_fluid = df_fluid.to_numpy().T

#

def mylogpolyfit(datax, datay):
    popt = np.polyfit(datax,np.log(datay),10)
    return np.poly1d(popt)

spec_fluidum_loginterp = mylogpolyfit(pTs_fluid, spec_fluid)
spec_fluidum_interp = np.exp(spec_fluidum_loginterp(pT_full))
spec_full += spec_fluidum_interp

COL = (1,0,0)
ax_full.plot(pT_full, spec_fluidum_interp,lw=LINEWIDTH,c="b",label=r"Fluid$\mathdefault{u}$m")
ax_full.fill_between(pT_full,spec_fluidum_interp,spec_full,facecolor=(*COL,0.2),edgecolor=COL,lw=LINEWIDTH,label=LEGEND)
ax_full.errorbar(pTs_exp, spec_exp, spec_exp_err,label="experiment",c="b",fmt="o",markersize=MARKERSIZE,lw=LINEWIDTH)

ax_full.set_yscale("log")
ax_full.set_ylim(3e1,2e3)
ax_full.set_xlim(0,pTmax)
ax_full.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)
ax_full.set_ylabel(SPEC_YLABEL, fontsize=AXISLABELSIZE)

ax_full.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax_full.xaxis.set_ticks_position("bottom")
ax_full.yaxis.set_ticks_position("left")
ax_full.grid(False)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8))
ax_full.yaxis.set_minor_locator(locmin)
ax_full.yaxis.set_minor_formatter(matplotlib.ticker.LogFormatterSciNotation(base=10,labelOnlyBase=False,minor_thresholds=(5,2.5))) # means: the data spans ~5 decades and we want to see all minor ticks if we zoom in on a region of 2.5 decades
ax_full.legend()

fig_full.tight_layout()
fig_full.show()

### DATA/MODEL PLOT

pTs_grid = np.linspace(0,1,50)
spec_fluidum_grid = np.exp(spec_fluidum_loginterp(pTs_grid))
spec_exp_loginterp = mylogpolyfit(pTs_exp, spec_exp)
spec_exp_grid = np.exp(spec_exp_loginterp(pTs_grid))

fig, ax = plt.subplots(figsize=(7,7))
ax.plot(pTs_grid, spec_exp_grid/spec_fluidum_grid)

ax.set_xlim(0,pTmax)
ax.set_xlabel(SPEC_XLABEL, fontsize=AXISLABELSIZE)
ax.set_ylabel("data/model", fontsize=AXISLABELSIZE)

ax.tick_params(axis="both",labelsize=TICKLABELSIZE)
ax.xaxis.set_ticks_position("bottom")
ax.yaxis.set_ticks_position("left")
ax.grid()

fig.show()

###
# %%
