#python3
#This program plots the SEDs of the objects in the root directory given. It also plots the IRAS LRS spectra and SDSS optical spectra if available.
#install colour package
#need root.inp and root_phot.csv on computer in same folder as this program
#run as follows:
#python phot_plotter.py root.inp
#ex .inp file: kastner08.inp
#this output is saved in root directory
#Major history upgrades
#8/28/21: Added proper motion data from Gaia in both root_coords.csv and dist-lum.c
#8/27/21: Included upper and lower limits of provided distances and thus the calculations of upper and lower limits of luminosity
#08/26/21: Changed logic so that the code will compute luminosity from input distances if needed
#08/25/21: Added flags in the .inp for when qout.dat needs to be used and if distances are given in the coords file
#06/29/21: changed logic of radial-vel distance computation for the case when d is given/not given in root.inp file
# TODOs: modify output root_dist_lum.csv file when d is given so that it does not list radial velocities and near/far distance, but just lists the input distance 
#09/17/20: Added interpolations and exterpolations for teff and phot_ioniz to calculate teff and nuFnu expected (added values to dist_lum.csv file)
#08/31/20: Added option to replot certain plots with different limits on the axes
#08/18/20: Added postage stamp images
#08/11/20: Added distance/luminosity calculations
#08/07/20: Added calculation for finding total flux using Simpsons rule and the trapezoidal rule
#08/06/20: Added downwards arrows on upper limit IRAS points
#08/03/20: Separated SED plotter and color-color/kmeans plotter into different programs. Color-color/kmeans program is named color-color_plotter.py
import numpy as np
import statistics
import matplotlib.pyplot as plt
import matplotlib
from matplotlib  import cm
import csv
import sys
import os
from pathlib import Path
import time
import random
import matplotlib.lines as mlines
from colour import Color
import mpld3
import astropy.units as u, astropy.constants as c
from astropy.coordinates import SkyCoord
from astropy.visualization import (MinMaxInterval, SqrtStretch, AsinhStretch, ImageNormalize)
from astropy.visualization import astropy_mpl_style
from urllib.parse import quote
from astropy.io import fits
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
#from kd import pdf_kd

start_time = time.time()

root = sys.argv[1][:-4]
coordfile = root + "_coords.csv"
fg = 0
ccfg = 0
count = 0
psfg =0
formats = []
hips_list = []
source = []
xmin = []
xmax = []
ymin = []
ymax = []
reps=[]
sscand = []
notsscand = []
lrad = []
con =0
con2 = 0
d = 0
c=2.9979e10
pi=3.14159
sunl=3.826e33
pc=3.085678e18
t_e=1e4
with open(sys.argv[1], 'r') as inp:
    lines = inp.readlines()
    for line in lines:
        if line[:1] == '#':
            continue
        else:
            index = line.find('=')
            if 'sep' in line[:index]:
                sep = str(line[index+1:].strip())
            if 'c_format' in line[:index]:
                form = str(line[index+1:].strip())
            if 'spectrumflag' in line[:index]:
                fg = int(line[index+1:].strip())
            if 'imageflag' in line[:index]:
                psfg = int(line[index+1:].strip())
            if 'hips_' in line[:index]:
                if line[index+1:].strip() == '1':
                    hips_list.append(str(line[5:index].strip()))
            if '.' in line[:index]:
                if line[index+1:].strip() == 'y':
                    formats.append(str(line[:index].strip()))
            if 'inp_cerr' in line[:index]:
                inp_cerr = float(line[index+1:].strip())
            if 'freq1_GHz' in line[:index]:
                freqone = float(line[index+1:].strip())
            if 'freq2_GHz' in line[:index]:
                freqtwo = float(line[index+1:].strip())
#             if 'source' in line[:index]:
#                 x = line[index+1:].strip().split(',')
#                 source.append(x[0])
#                 xmin.append(float(x[1]))
#                 xmax.append(float(x[2]))
#                 ymin.append(float(x[3]))
#                 ymax.append(float(x[4]))
            if 'distance' in line[:index]:
                d = float(line[index+1:].strip())
            if 'dist_each_source' in line[:index]:
                distfg = int(line[index+1:].strip())
            if 'massivestar' in line[:index]:
                qout = int(line[index+1:].strip())
path = str(Path(__file__).parent.absolute())
root = root+'-'+str(inp_cerr)+"asec"
dir = path + '/' + str(root)
ty = []
distances = []
with open(coordfile, 'r') as c:
    lines = c.readlines()
    for line in lines:

        if line[:1] == "#":
            continue
        else:
            x = line.split(',')
            if distfg == 1:
                distances.append(x[3:])
            #ty.append(str(x[3]))
#postage stamp images
if psfg ==1:

    objects = []
    raa = []
    decc = []
    with open(coordfile, 'r') as c:
        lines = c.readlines()
        for line in lines:

            if line[:1] == "#":
                continue
            else:
                x = line.split(sep)
                objects.append(str(x[0].strip()))
                if form == 'd':
                    raa.append(float(x[1]))
                    decc.append(float(x[2]))
                elif form =='s':
                    c = SkyCoord(str(x[1])+' '+str(x[2]), unit=(u.hourangle, u.deg))
                    raa.append(float(c.ra.degree))
                    decc.append(float(c.dec.degree))

    width = 150
    height = 150
    fov = 0.1

    nb_obj = len(objects)
    nb_hips = len(hips_list)

    fig, axs = plt.subplots(nb_obj, nb_hips, figsize=(4 * nb_hips, 3* nb_obj), facecolor='w', edgecolor='k')
    fig.subplots_adjust(hspace = .2, wspace=.001)

    axs = axs.ravel()

    subdir = dir + '/fitsfiles'
    os.makedirs(os.path.dirname(subdir+"/files.txt"), exist_ok=True)
    i = 0
    for j in range(len(objects)):
        for hips in hips_list:
            axs[i].set_axis_off()

            axs[i].set_title('{} - {}'.format(objects[j], hips))
            #sc = SkyCoord.from_name(objects[i])
            ra = raa[j]
            dec = decc[j]
            url = 'http://alasky.u-strasbg.fr/hips-image-services/hips2fits?hips={}&width={}&height={}&fov={}&projection=TAN&coordsys=icrs&ra={}&dec={}'.format(quote(hips), width, height, fov, ra, dec)

            hdu = fits.open(url)

            file_name = '{}-{}.fits'.format(objects[j].replace(' ',"_"), hips.replace('/', '_'))
            hdu.writeto(subdir + "/" + file_name, overwrite=True)
            print('Saving {}'.format(file_name))
            im = hdu[0].data
            norm = ImageNormalize(im, interval=MinMaxInterval(),
                          stretch=AsinhStretch())
            axs[i].imshow(im, cmap='magma', norm=norm, origin='lower')

            i += 1


data = dir + "/" + root + "_phot.csv"


#print(result)
integ = dir + '/' + root + "_integresults.txt"
with open(integ,"w") as r:
    r.write('Integration Results \n')
r.close()
info = []
with open(data, mode='r') as phot:
    phot_reader = csv.reader(phot, delimiter=',')
    #next(phot_reader)
    for row in phot_reader:
        info.append(row)

lsten = []
temp = []
ioniz_fl = []

if qout==1:
    with open('qout.dat', 'r') as dat: #this is a file of Teff, Lum vs #ionzing photons from models of massive stars, don't need it if we input sources that are not massive stars
        lines = dat.readlines()[1:]
        for line in lines:
            x1 = line.split(' ')
            x = []
            for i in range(len(x1)):
                if x1[i]!='':
                    x.append(x1[i])
            if x[1].strip() == '4.00':

                temp.append(float(x[0].strip()))
                lsten.append(float(x[3].strip()))
                flu = float(x[7].strip())+float(x[8].strip())
                ioniz_fl.append(flu*10**(-48))

def findind(st):
    if st == 'Landolt':
        return lab.index('LANDOLT')
    elif st == 'Sky':
        return lab.index('SKYMAPPER')
    elif st == "Gaia":
        return lab.index('GAIA/GAIA2')
    elif st == 'Cousins':
        return lab.index('COUSINS')
    elif st == 'Johnson':
        return lab.index('JOHNSON')
    elif st == 'glimpse':
        return lab.index('glimpse(spitzer)')
    elif st == 'WISE':
        return lab.index('WISE')
    elif st == 'ISO':
        return lab.index('ISO')
    elif st == 'mips':
        return lab.index('mipsgal(spitzer)')
    elif st == 'ppsc':
        return lab.index('PACS')
    elif st == "Herschel":
        return lab.index('HERSCHEL')
    elif st == 'spsc':
        return lab.index('SPIRE')
    elif st == 'ATLASGAL':
        return lab.index('atlasgal(spitzer)')
    elif st == 'HI-GAl':
        return lab.index('higal(spitzer)')
    elif st == 'other':
        return len(lab)-1
    else:
        for p in range(len(lab)):
            if st in lab[p]:
                return p

def sortSecond(val):
    return val[0]
r0 = 8.0
v0 = 220.0
def dist(l, b, vr):
# kd calculates kinematic distances and kinematic distance uncertainties
# installed as follows from https://github.com/tvwenger/kd
#  pip3 install git+https://github.com/tvwenger/pyqt-fit.git
#  pip3 install git+https://github.com/tvwenger/kd.git
    from kd import pdf_kd
    glong = float(l) # Galactic longitude, degrees
    velo = float(vr) # measured LSR velocity, km/s
    velo_err = 0 # measured LSR velocity uncertainty, km/s
    rotcurve = 'reid14_rotcurve' # the name of the script containing the rotation curve
    num_samples = 2000 # number of re-samples
    glat = float(b) # Galactic latitude, degrees
    d = pdf_kd.pdf_kd(glong, glat, velo, velo_err=velo_err, rotcurve=rotcurve, num_samples=num_samples)
    return d['near'], d['far']
def lum(lglum):
    alum=np.exp(2.302585*lglum)
    return alum
def myExpFunc(x, a, b):
    return a * np.power(x, b)
def coeffind(a, x, y):
    return y/(x**a)
def tef(lllum,dkpc):
    c=2.9979e10
    pi=3.14159
    sunl=3.826e33
    pc=3.085678e18
    t_e=1e4

    teff=np.empty(3)
    teff[0]=22600.0
    teff[1]=20500.0
    teff[2]=17900.0

    # luminosities log (L [lsun])
    alum_lg=np.empty(3)
    alum_lg[0]=3.72
    alum_lg[1]=3.46
    alum_lg[2]=3.02

    # these are in units of 1e48 /s
    phot_ioniz=np.empty(3)
    phot_ioniz[0]=0.0033
    phot_ioniz[1]=4.47e-4
    phot_ioniz[2]=4.898e-5

    stlin = interp1d(lsten, temp, kind='linear')
    st2lin = interp1d(lsten, ioniz_fl, kind = 'linear')

    lums = []
    for i in range(len(alum_lg)):
        lums.append(lum(alum_lg[i]))

    popt2, pcov2 = curve_fit(myExpFunc, lums, teff)
    cof = coeffind(popt2[1], lsten[0], temp[0])

    popt, pcov = curve_fit(myExpFunc, lums, phot_ioniz)

    cof2 = coeffind(popt[1], lsten[0], ioniz_fl[0])
    dcm=1.e3*float(dkpc)*pc
    fbol = float(lllum)*(sunl/dcm)/(4*pi*dcm)
    print(min(lums),min(lsten),max(lsten))
    if float(lllum) >= min(lums) and float(lllum) <= min(lsten):
        fion = float(cof2)*float(lllum)**float(popt[1])
    elif float(lllum) > min(lsten) and float(lllum)<= max(lsten):
        fion = float(st2lin(float(lllum)))
    else:
        fion = '-999'

    if float(lllum) >= min(lums) and float(lllum) <= min(lsten):
        #fion = cof2*llum**popt[1]
        teff = float(cof)*float(lllum)**float(popt2[1])

    elif float(lllum) > min(lsten) and float(lllum)<= max(lsten):
        #fion = float(st2lin(llum))
        teff = float(stlin(float(lllum)))
    else:
        #fion = '-999'
        teff = '-999'
    return(teff, fion)

def exprat(dkpc,fnu,llum):


    # luminosities log (L [lsun])
    alum_lg=np.empty(3)
    alum_lg[0]=3.72
    alum_lg[1]=3.46
    alum_lg[2]=3.02

    # these are in units of 1e48 /s
    phot_ioniz=np.empty(3)
    phot_ioniz[0]=0.0033
    phot_ioniz[1]=4.47e-4
    phot_ioniz[2]=4.898e-5

    st2lin = interp1d(lsten, ioniz_fl, kind = 'linear')

    lums = []
    for i in range(len(alum_lg)):
        lums.append(lum(alum_lg[i]))

    popt, pcov = curve_fit(myExpFunc, lums, phot_ioniz)

    cof2 = coeffind(popt[1], lsten[0], ioniz_fl[0])

    # Aversa+11 (AJ, 141, 125), eqn 2
    dcm=1.e3*float(dkpc)*pc

    fbol = float(llum)*(sunl/dcm)/(4*pi*dcm)
    if float(llum) >= min(lums) and float(llum) <= min(lsten):
        fion = float(cof2)*float(llum)**float(popt[1])
    elif float(llum) > min(lsten) and float(llum)<= max(lsten):
        fion = float(st2lin(float(llum)))
    else:
        fion = '-999'

    if fion != '-999':
        fradio=((fion/6.3e4)*(1e27/dcm)/(4.*pi*dcm))/((t_e/1e4)**(-0.45)*(float(fnu)**0.1))
        fradio_mjy=fradio*1e26
        freq_flux = float(fnu)*1e9*fradio
        ratio = freq_flux/fbol
    else:
        ratio = '-999'

    return(ratio)
#set color key
lab = ["GALEX", "HIP", "LANDOLT","SKYMAPPER", "APASS(sdss)","SDSS", "POSS", "GAIA/GAIA2", "HST/ACS","INT", "PAN-STARRS", "XMM-OT", "COUSINS", "DENIS", "JOHNSON", "VISTA", "2MASS", "UKIDSS", "UKIRT", 'OASIS', 'EIC', "DIRBE", "SAGE", "CATWISE","glimpse(spitzer)", 'CELOBJ(HALL74)', "WISE", "ALLWISE", "ISO", "MSX", "ISOGAL", 'HII-REGIONS(GIVEON05)', "mipsgal(spitzer)", "1RXS","IRAS", 'IRDC', "AKARI", "PACS", 'HERSCHEL', "PN-HIIregions(ANDERSON12)", "HI-GAL-CORNISH(HIIREGIONS)", "higal(spitzer)", "SPIRE", "BLAST", 'SCUBA', "JPS", "atlasgal(spitzer)", 'W31(BEUTHER13)', "PLANCK", "BOLOCAM", 'SEST', "SPECFIND", "GLOSTAR", "RMS-Survey", "3FGL-ATCA-VLA(SCHINZEL17)","AT20G","CRATES","VLA-INNERGAL(BECKER94)",'MAGPIS', "CORNISH", 'VLA-OUTERGAL(FICH86)', "THOR", "FIRST-NVSS(VAR)",  "SUMSS",  "other"]
print(len(lab))
labwav = []
red = Color("red")
blue = Color("blue")
v = list(blue.range_to(red, len(lab)))
for i in range(len(v)):
    labwav.append(v[i].rgb)
mark = []
for i in range(len(lab)):
    if (i % 5) == 0:
        mark.append('s')
    elif (i%5) == 1:
        mark.append('v')
    elif (i%5) == 2:
        mark.append('P')
    elif (i%5) == 3:
        mark.append('^')
    else:
        mark.append('o')
def marks(k):
    if (k % 5) == 0:
        st = 's'
        return st
    elif (k%5) == 1:
        st = 'v'
        return st
    elif (k%5) == 2:
        st = 'P'
        return st
    elif (k%5) == 3:
        st = '^'
        return st
    else:
        st = 'o'
        return st
style = []
for i in range(len(lab)):
    style.append(mlines.Line2D([], [], markeredgecolor=labwav[i], markerfacecolor = 'None', marker=mark[i], linestyle='None',
                          markersize=5, label=lab[i]))

lums = []
rats = []
hlen = 0
freqlist = []
fils = []
ob = []
t = []
#print(len(ty), len(info), 'hi')
for i in range(1,len(info)):
#for i in range(1,25):
    vels = []
    inte = []
    inte.append('\n')
    inte.append(info[i][0] + "\n")
    o = info[i][0]
    ob.append(o)
    obj = info[i][0].replace(" ", "_")
    ra = info[i][1]
    dec = info[i][2]
    type = info[i][3]
    if distfg == 1:
        dist_kpc = float(distances[i-1][0])
        dist_kpc_lower = float(distances[i-1][1])
        dist_kpc_higher = float(distances[i-1][2])
        pm = float(distances[i-1][3])
        pnra = float(distances[i-1][4])
        e_pmra = float(distances[i-1][5])
        pmde = float(distances[i-1][6])
        e_pmde = float(distances[i-1][7])
        excess_noise = float(distances[i-1][8])
        excess_noise_sig = float(distances[i-1][9])
        ipd_gof_ha = float(distances[i-1][10])
        ruwe = float(distances[i-1][11])
    #type = ty[i-1]
    #t.append(type)
    specfile_op = -999
    specfile_iras = '-999'
    wavelength = []
    fnu = []
    errorfnu = []
    filter = []
    cata = []
    for j in range(len(info[i])):
        if str(info[0][j][-18:]) == "wavelength(micron)" and "Lin" not in str(info[0][j]) and "Sp:" not in str(info[0][j]) and 'Lin@' not in str(info[0][j]):
            #print(info[i][j])
            #print(obj)
            if float(info[i][j]) != -999:
                wavelength.append(float(info[i][j]))
                fnu.append(float(info[i][j+1]))
                filter.append(info[0][j][:-19])
                if info[i][j+2] == "--":
                    errorfnu.append(0)
                else:
                    errorfnu.append(float(info[i][j+2]))
        elif "Lin" in str(info[0][j]) or 'Lin@' in str(info[0][j]):
            if 'catalogue' in str(info[0][j]):
                cat = info[i][j]
            if '_VLSR(km/s)' in str(info[0][j]) or '_VPEAK(km/s)' in str(info[0][j]) or '_vel(km/s)' in str(info[0][j]) or '_pkvel(km/s)' in str(info[0][j]) or '_minvel(km/s)' in str(info[0][j]) or "_maxvel(km/s)" in str(info[0][j]) or '_vpeak(km/s)' in str(info[0][j]) or '_vpeak2(km/s)' in str(info[0][j]):
                if str(info[i][j]) != '-999' and str(info[i][j]) != '--':
                    vels.append(str(info[i][j]) + '-' + str(info[0][j]))
                    cata.append(cat)
                    if str(info[0][j]) not in fils:
                        fils.append(str(info[0][j]))
            elif '_VL(km/s)' in str(info[0][j]):
                if str(info[i][j]) != '-999' and str(info[i][j]) != '--' and str(info[i][j+1]) != '-999' and str(info[i][j+1]) != '--' and "'" not in str(info[i][j+1]) and "'" not in str(info[i][j]):
                    avg = (float(info[i][j]) + float(info[i][j+1]))/2
                    ind = info[0][j].find('_')
                    vels.append(str(avg) + '-' + info[0][j][0:ind] + '_vavg(km/s)')
                    cata.append(cat)
                    if info[0][j][0:ind] + '_vavg(km/s)' not in fils:
                        fils.append(info[0][j][0:ind] + '_vavg(km/s)')
                elif str(info[i][j]) != '-999' and str(info[i][j]) != '--':
                    if str(info[i][j+1]) == '-999' or str(info[i][j+1]) == '--':
                        vels.append(str(info[i][j]) + '-' + str(info[0][j]))
                        cata.append(cat)
                        if str(info[0][j]) not in fils:
                            fils.append(str(info[0][j]))
                elif str(info[i][j+1]) != '-999' and str(info[i][j+1]) != '--':
                    if str(info[i][j]) == '-999' or str(info[i][j]) == '--':
                        vels.append(str(info[i][j+1]) + '-' + str(info[0][j+1]))
                        cata.append(cat)
                        if str(info[0][j+1]) not in fils:
                            fils.append(str(info[0][j+1]))
        elif "Sp:IRAS" in str(info[0][j]):
            specfile_iras = str(info[i][j])
        elif "Sp:SDSS_optical" in str(info[0][j]):
            specfile_op = str(info[i][j])

    print("This is fnu:" + str(fnu))
    
    nufnu = []
    errornufnu = []
    nu = []
    for k in range(len(fnu)):
        nu.append((2.99*10**(14))/wavelength[k])
        nufnu.append(nu[k]*fnu[k]*10**(-23))
        errornufnu.append(nu[k]*errorfnu[k]*10**(-23))

    num = []
    col = []
    markstyle = []
    for i in range(len(filter)):
        if 'GALEX' in filter[i]:
            st = 'GALEX'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'HIP' in filter[i]:
            st = 'HIP'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'Landolt' in filter[i]:
            st = 'Landolt'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'Sky' in filter[i]:
            st = 'Sky'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif "APASS" in filter[i]:
            st = 'APASS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'SDSS' in filter[i]:
            st = 'SDSS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'POSS' in filter[i]:
            st = 'POSS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'Gaia' in filter[i] or 'GAIA' in filter[i]:
            st = 'Gaia'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'HST/ACS' in filter[i]:
            st = 'HST/ACS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'INT' in filter[i]:
            st = 'INT'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'PAN' in filter[i]:
            st = 'PAN'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'XMM-OT' in filter[i]:
            st = 'XMM-OT'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'Cousins' in filter[i]:
            st = 'Cousins'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'DENIS' in filter[i]:
            st = 'DENIS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'Johnson' in filter[i]:
            st = 'Johnson'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'VIST' in filter[i]:
            st = 'VIST'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif '2MASS' in filter[i]:
            st = '2MASS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'UKIDSS' in filter[i]:
            st = 'UKIDSS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'UKIR' in filter[i]:
            st = 'UKIR'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'OASIS' in filter[i]:
            st = 'OASIS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'EIC' in filter[i]:
            st = 'EIC'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'DIRBE' in filter[i]:
            st = 'DIRBE'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'SAGE' in filter[i]:
            st = 'SAGE'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'CATWISE' in filter[i]:
            st = 'CATWISE'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'glimpse' in filter[i]:
            st = 'glimpse'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'HALL74' in filter[i]:
            st = 'HALL74'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'WISE' in filter[i] and "ALLWISE" not in filter[i]:
            st = 'WISE'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'ALLWISE' in filter[i]:
            st = 'ALLWISE'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'ISO' in filter[i] and 'GAL' not in filter[i]:
            st = 'ISO'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'MSX' in filter[i]:
            st = 'MSX'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'ISOGAL' in filter[i]:
            st = 'ISOGAL'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'GIVEON05' in filter[i]:
            st = 'GIVEON05'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'mips' in filter[i]:
            st = 'mips'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif '1RXS' in filter[i]:
            st = '1RXS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'IRAS' in filter[i]:
            st = 'IRAS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'IRDC' in filter[i]:
            st = 'IRDC'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'AKAR' in filter[i]:
            st = 'AKAR'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'ppsc' in filter[i] or 'PACS' in filter[i]:
            st = 'ppsc'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'Herschel' in filter[i]:
            st = 'Herschel'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'PN-HIIREGIONS(ANDERSON12)' in filter[i]:
            st = 'PN-HIIregions(ANDERSON12)'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'HI-GAL-CORNISH' in filter[i]:
            st = 'HI-GAL-CORNISH'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'HI-GAL' in filter[i] and 'CORNISH' not in filter[i]:
            st = 'HI-GAL'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'spsc' in filter[i] or 'SPIRE' in filter[i]:
            st = 'spsc'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'BLAST' in filter[i]:
            st = 'BLAST'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'SCUBA' in filter[i]:
            st = 'SCUBA'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'JPS' in filter[i]:
            st = 'JPS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'ATLASGAL' in filter[i]:
            st = 'ATLASGAL'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'W31' in filter[i]:
            st = 'W31'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'PLANCK' in filter[i]:
            st = 'PLANCK'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'BOLOCAM' in filter[i]:
            st = 'BOLOCAM'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'SEST' in filter[i]:
            st = 'SEST'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'SPECFIND' in filter[i]:
            st = 'SPECFIND'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'GLOSTAR' in filter[i]:
            st = 'GLOSTAR'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'RMS-Survey' in filter[i]:
            st = 'RMS-Survey'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif '3FGL' in filter[i]:
            st = '3FGL'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'AT20G' in filter[i]:
            st = 'AT20G'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'CRATES' in filter[i]:
            st = 'CRATES'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'VLA-INNERGAL' in filter[i]:
            st = 'VLA-INNERGAL'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'MAGPIS' in filter[i]:
            st = 'MAGPIS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'CORNISH' in filter[i]:
            st = 'CORNISH'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'VLA-OUTERGAL(FICH86)' in filter[i]:
            st = 'VLA-OUTERGAL'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'THOR' in filter[i]:
            st = 'THOR'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'FIRST-NVSS' in filter[i]:
            st = 'FIRST-NVSS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        elif 'SUMSS' in filter[i]:
            st = 'SUMSS'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
        else:
            st = 'other'
            k = findind(st)
            col.append(labwav[k])
            markstyle.append(marks(k))
            if k not in num:
                num.append(k)
    num = sorted(num)
    specstyle = []

    #label = []
    #labelwav = []
    for i in range(len(num)):
        #label.append(lab[num[i]])
        #labelwav.append(labwav[num[i]])
        specstyle.append(style[num[i]])

    #iras spectra data
    if str(specfile_iras) != "-999" and str(specfile_iras) != '':
        lam = []
        fl = []
        lam2 = []
        fl2 = []
        flag = 0
        with open(specfile_iras, mode='r') as sp:
            reader = csv.reader(sp, delimiter=',')
            next(reader)
            #print(reader[1])
            for row in reader:
                #print(row)
                if len(lam) ==0:
                    lam.append(float(row[0]))
                    fl.append(float(row[1]))
                elif float(row[0]) > lam[-1] and flag ==0:
                    lam.append(float(row[0]))
                    fl.append(float(row[1]))
                elif float(row[0]) < lam[-1] and flag ==0:
                    flag = 1
                    lam2.append(float(row[0]))
                    fl2.append(float(row[1]))
                if flag ==1:
                    lam2.append(float(row[0]))
                    fl2.append(float(row[1]))
        sp.close()
        freq = []
        freq2 = []
        flcor = []
        nfl = []
        flcor2 = []
        nfl2 = []
        for m in range(len(lam)):
            freq.append(2.9979*10**(14)/lam[m])
            factors = 10**(26)*lam[m]*lam[m]*10**(-4)/(2.9979*10**(10))
            flcor.append(fl[m]*factors)
            nfl.append(freq[m]*flcor[m]*10**(-23))
        for m in range(len(lam2)):
            freq2.append(2.9979*10**(14)/lam2[m])
            factors = 10**(26)*lam2[m]*lam2[m]*10**(-4)/(2.9979*10**(10))
            flcor2.append(fl2[m]*factors)
            nfl2.append(freq2[m]*flcor2[m]*10**(-23))

    #optical spectra data
    if str(specfile_op) != "-999":
        lamb = []
        fluxop = []
        with open(specfile_op, mode='r') as s:
            read = csv.reader(s, delimiter=',')
            next(read)
            for row in read:
                lamb.append(float(row[0]))
                fluxop.append(float(row[1]))
        s.close()
        opfl = []
        for m in range(len(lamb)):
            lamang = lamb[m]/0.0001
            factors = 10**(26)*lamb[m]*lamb[m]*10**(-8)/(2.9979*10**(10))
            opfl.append(factors*fluxop[m])

        oplam = lamb

    if o in source:
        ind = source.index(o)
        plt.figure()
        ax = plt.axes()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(xmin[ind],xmax[ind])
        ax.set_ylim(ymin[ind],ymax[ind])
        ax.grid(True, alpha=0.3)
        for j in range(len(wavelength)):
            ax.errorbar(wavelength[j], fnu[j], yerr=errorfnu[j], capsize=5, ecolor = col[j], linestyle="None", elinewidth=1)
            ax.scatter(wavelength[j],fnu[j], edgecolors=col[j], s=20, marker = markstyle[j], facecolors='none', linewidth=.8)
            if '-uplim' in filter[j]:
                #ax.scatter(wavelength[j] ,fnu[j], c='purple',marker= r'$\downarrow$' ,s=40, label='arrow' )
                annotate = ax.annotate('', xy=(wavelength[j], fnu[j]), xytext=(0, -25), textcoords='offset points', arrowprops=dict(arrowstyle="<|-"), label = 'upperlimit')
        leg = ax.legend(handles=specstyle, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})
        #ax.legend(handles = [annotate], handler_map={type(annotate) : AnnotationHandler(5)})
        #legar =ax.legend(loc='center left', bbox_to_anchor=(1, 0.2), prop={'size': 8})
        ax.add_artist(leg)
        if str(specfile_iras) != "-999":
            ax.plot(lam,flcor, '-g.', markersize=2, linewidth = 1, label = 'IRAS LRS Spectra')
            ax.plot(lam2,flcor2, '-r.', markersize=2, linewidth = 1, label = 'IRAS LRS Spectra')
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.2), prop={'size': 8})
        plt.xlabel(r'Wavelength$\,(\mu m)$', size=14)
        plt.ylabel(r'F$(\nu)\,$(Jy)', size=14)
        plt.title(obj + ' (' + ra + ", " + dec + ")")
        plt.tight_layout()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width, box.height])
        # Put a legend to the right of the current axis

        #plt.legend(recs,lab,loc='best', prop={'size': 5})
        for j in range(len(formats)):
            plt.savefig(dir + '/' + obj + '_fnu-wav_limitaxes' + formats[j])
        plt.close()

        plt.figure()
        ax = plt.axes()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(xmin[ind],xmax[ind])
        ax.set_ylim(ymin[ind],ymax[ind])
        ax.grid(True, alpha=0.3)
        for j in range(len(wavelength)):
            ax.errorbar(wavelength[j], nufnu[j], yerr=errornufnu[j], capsize=5, ecolor = col[j], linestyle="None", elinewidth=1)
            ax.scatter(wavelength[j],nufnu[j],s=20,edgecolors=col[j], marker = markstyle[j], facecolors='none', linewidth=.8)
            if '-uplim' in filter[j]:
                ax.annotate('', xy=(wavelength[j], nufnu[j]), xytext=(0, -25), textcoords='offset points', arrowprops=dict(arrowstyle="<|-"))
        leg1 = ax.legend(handles=specstyle, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})
        ax.add_artist(leg1)
        if str(specfile_iras) != "-999":
            grsp = ax.plot(lam,nfl, '-g.', markersize=2, linewidth = 1, label = 'IRAS LRS Spectra')
            rsp = ax.plot(lam2,nfl2, '-r.', markersize=2, linewidth = 1, label = 'IRAS LRS Spectra')
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.2), prop={'size': 8})
            # Manually add the first legend back
        plt.xlabel(r'Wavelength$\,(\mu m)$', size=14)
        plt.ylabel(r'$\nu$F$(\nu)\,$(erg/cm$^2$/s)', size=14)
        plt.title(obj + ' (' + ra + ", " + dec + ")")
        plt.tight_layout()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width , box.height])
        for j in range(len(formats)):
            plt.savefig(dir + '/' + obj + '_nufnu-wav_limitaxes' +formats[j])
        plt.close()

    plt.figure()
    ax = plt.axes()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    for j in range(len(wavelength)):
        ax.errorbar(wavelength[j], fnu[j], yerr=errorfnu[j], capsize=5, ecolor = col[j], linestyle="None", elinewidth=1)
        ax.scatter(wavelength[j],fnu[j], edgecolors=col[j], s=20, marker = markstyle[j], facecolors='none', linewidth=.8)
        if '-uplim' in filter[j]:
            #ax.scatter(wavelength[j] ,fnu[j], c='purple',marker= r'$\downarrow$' ,s=40, label='arrow' )
            annotate = ax.annotate('', xy=(wavelength[j], fnu[j]), xytext=(0, -25), textcoords='offset points', arrowprops=dict(arrowstyle="<|-"), label = 'upperlimit')
    leg = ax.legend(handles=specstyle, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})
    #ax.legend(handles = [annotate], handler_map={type(annotate) : AnnotationHandler(5)})
    #legar =ax.legend(loc='center left', bbox_to_anchor=(1, 0.2), prop={'size': 8})
    ax.add_artist(leg)
    if str(specfile_iras) != "-999":
        ax.plot(lam,flcor, '-g.', markersize=2, linewidth = 1, label = 'IRAS LRS Spectra')
        ax.plot(lam2,flcor2, '-r.', markersize=2, linewidth = 1, label = 'IRAS LRS Spectra')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.2), prop={'size': 8})
    plt.xlabel(r'Wavelength$\,(\mu m)$', size=14)
    plt.ylabel(r'F$(\nu)\,$(Jy)', size=14)
    plt.title(obj + ' (' + ra + ", " + dec + ")")
    plt.tight_layout()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height])
    # Put a legend to the right of the current axis

    #plt.legend(recs,lab,loc='best', prop={'size': 5})
    for j in range(len(formats)):
        plt.savefig(dir + '/' + obj + '_fnu-wav' + formats[j])
    plt.close()

    #nu flux vs wavelength
    plt.figure()
    ax = plt.axes()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    for j in range(len(wavelength)):
        ax.errorbar(wavelength[j], nufnu[j], yerr=errornufnu[j], capsize=5, ecolor = col[j], linestyle="None", elinewidth=1)
        ax.scatter(wavelength[j],nufnu[j],s=20,edgecolors=col[j], marker = markstyle[j], facecolors='none', linewidth=.8)
        if '-uplim' in filter[j]:
            ax.annotate('', xy=(wavelength[j], nufnu[j]), xytext=(0, -25), textcoords='offset points', arrowprops=dict(arrowstyle="<|-"))
    leg1 = ax.legend(handles=specstyle, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8})
    ax.add_artist(leg1)
    if str(specfile_iras) != "-999":
        grsp = ax.plot(lam,nfl, '-g.', markersize=2, linewidth = 1, label = 'IRAS LRS Spectra')
        rsp = ax.plot(lam2,nfl2, '-r.', markersize=2, linewidth = 1, label = 'IRAS LRS Spectra')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.2), prop={'size': 8})
        # Manually add the first legend back
    plt.xlabel(r'Wavelength$\,(\mu m)$', size=14)
    plt.ylabel(r'$\nu$F$(\nu)\,$(erg/cm$^2$/s)', size=14)
    plt.title(obj + ' (' + ra + ", " + dec + ")")
    plt.tight_layout()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width , box.height])
    for j in range(len(formats)):
        plt.savefig(dir + '/' + obj + '_nufnu-wav' +formats[j])
    plt.close()

    #spectra plot separately
    if str(specfile_iras) != "-999" and fg ==1:

        plt.figure()
        ax = plt.axes()
        #ax.set_xscale("log")
        #ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        #ax.scatter(lam,fl, s=20)
        ax.plot(lam,flcor, '-g.')
        ax.plot(lam2,flcor2, '-r.')
        plt.xlabel(r'Wavelength$\,(\mu m)$', size=14)
        plt.ylabel(r'Flux (Jy)', size=14)
        plt.title(obj + ' (' + ra + ", " + dec + ") Spectrum")
        for j in range(len(formats)):
            plt.savefig(dir + '/' + obj + '_spectrum' +formats[j])
        plt.close()
    if str(specfile_op) != "-999" and fg ==1:

        plt.figure()
        ax = plt.axes()
        #ax.set_xscale("log")
        #ax.set_yscale("log")
        ax.grid(True, alpha=0.3)
        #ax.scatter(lam,fl, s=20)
        ax.plot(oplam,opfl, '-b.')
        plt.xlabel(r'Wavelength$\,(\mu m)$', size=14)
        plt.ylabel(r'Flux (Jy)', size=14)
        plt.title(obj + ' (' + ra + ", " + dec + ") Optical Spectrum")
        for j in range(len(formats)):
            plt.savefig(dir + '/' + obj + '_op_spectrum' +formats[j])
        plt.close()

    #integrate under the curve
    #f = nu
    #fden = fnu

    from scipy.integrate import simps

    if len(fnu) != 0:
        f = []
        fden =[]
        nufden = []
        for j in range(len(nu)):
            if nu[j] not in f and 'uplim' not in filter[j]:
                f.append(nu[j])
                fden.append(fnu[j])
                nufden.append(nufnu[j])

        #print("This is fden:" + str(fden))
        
        fsort = (sorted(f))

        fdensort = [x for _,x in sorted(zip(f,fden))]
        nufdensort = [x for _,x in sorted(zip(f,nufden))]
        
        print("This is fdensort:" + str(fdensort))
        #print("This is nufdensort:" + str(nufdensort))
        print("This is fsort:" + str(fsort))

        simpson = simps(fdensort, fsort, even='avg')
        trap = np.trapz(fdensort,fsort)

        fsort1 = []
        fdensort1 = []
        fsort2 = []
        fdensort2 = []

        fsortmaxind = nufdensort.index(max(nufdensort))
        for j in range(fsortmaxind+1):
            fsort1.append(fsort[j])
            fdensort1.append(fdensort[j])
        for j in range(fsortmaxind,len(fsort)):
            fsort2.append(fsort[j])
            fdensort2.append(fdensort[j])

        robsim = simps(fdensort1, fsort1, even='first') + simps(fdensort2, fsort2, even='last')

        print("This is trap:" + str(trap))
        print("This is robsim:" + str(robsim))
        print("This is simpson:" + str(simpson))
        
        #if robsim < 0:
        #    finflux = trap
        #else:
        #    finflux = robsim
        if trap > 0:
            finflux = trap
        elif robsim > 0:
            finflux = robsim
        elif simpson > 0:
            finflux = simpson
        else:
            finflux = -999

        if finflux==-999:
            finflux = -999
        else:
            finflux = finflux *(10**(-23))
        simpson = simpson*(10**(-23))
        robsim = robsim *(10**(-23))
        trap = trap* (10**(-23))
    else:
        simpson = -999
        trap = -999
        robsim =-999
        finflux = -999

    units = 'erg/cm^2/s'

    inte.append('Robust Simpsons rule: ' + str(robsim) + " " + units + "\n")
    inte.append("Simpsons rule: " + str(simpson) + " " + units + "\n")
    inte.append("Trapezoidal rule: " + str(trap) + " " + units + "\n")
    inte.append("Total Flux: " + str(finflux) + " " + units + "\n")

    with open(integ, "a") as i:
        for j in range(len(inte)):
            i.write(inte[j])
    i.close()
    #finding distance and luminosities using flux
    r = []
        #cols: source name, coords (both), bol flux, vr, 2 distances, 2 lum
    r.append(o)
    r.append(type)
    r.append(ra)
    r.append(dec)

    c_icrs = SkyCoord(ra=float(ra)*u.degree, dec=float(dec)*u.degree, frame='icrs')
    c_gal =c_icrs.galactic
    l = c_gal.l.degree
    b = c_gal.b.degree

    r.append(l)
    r.append(b)
    r.append(finflux)
    print(finflux)
    av = []
    av2=[]
    vs = []
    a2 = '-999'
    if distfg == 0:
        for j in range(len(vels)):
            print('# radial velocities ',vels[j])
            if vels[j][0] == '-':
                #print('y')
                index = vels[j].find('-', 2)
            else:
                index = vels[j].find('-', 1)
            vel = vels[j][:index]
            print(vel)
            if vel[-1] == ')':
                vel = vel[:-1]
            print(vel)
            vs.append(float(vel))
            r.append(cata[j])
            r.append(vels[j]) 
            if d !=0:
                dnear = d
                dfar = d
            else:
                dnear, dfar = dist(l, b, vel)
                if 'nan' in str(dnear):
                    dnear = '-999'
                if 'nan' in str(dfar):
                    dfar = '-999'

            r.append(dnear)
            r.append(dfar)
            if dnear != '-999':
                dnearcm = dnear * 3.086e+21
            else:
                dnearcm = '-999'
            if dfar != '-999':
                dfarcm = dfar * 3.086e+21
            else:
                dfarcm = '-999'
            if finflux >0:
                if dnear!='-999':
                    lumnear = 4*np.pi*dnearcm**2*finflux
                    lumnearsol = lumnear / 3.839e+33

                else:
                    lumnear = -999
                    lumnearsol = -999
                if dfar!='-999':
                    lumfar = 4*np.pi*dfarcm**2*finflux
                    lumfarsol = lumfar / 3.839e+33
                else:
                    lumfar = -999
                    lumfarsol = -999
            else:
                lumnear = -999
                lumfar = -999
                lumnearsol = -999
                lumfarsol = -999


            av.append(lumnearsol)
            av.append(lumfarsol)
            r.append(lumnearsol)
            r.append(lumfarsol) 
            #    if len(vs) == 0 and d!=0: possibly incorrect logic, we want to use distance given in inp file, whether or not there is radial vel info (stored in vs array)
            if d!=0:
        #        print('yy')
                r.append('-999')
                r.append('-999')
                dnear = d
                dfar = d
                r.append(dnear)
                r.append(dfar)
                dnearcm = dnear * 3.086e+21
                dfarcm = dfar * 3.086e+21
                if finflux >0:
                    lumnear = 4*np.pi*dnearcm**2*finflux
                    lumfar = 4*np.pi*dfarcm**2*finflux
                    lumnearsol = lumnear / 3.839e+33
                    lumfarsol = lumfar / 3.839e+33
                else:
                    lumnear = -999
                    lumfar = -999
                    lumnearsol = -999
                    lumfarsol = -999
                av.append(lumnearsol)
                av.append(lumfarsol)
                #print(obj, av, finflux)
                a2 = d
                r.append(lumnearsol)
                r.append(lumfarsol)
        

    #print('Source, AllLumVals ',obj,av)
#    if len(vs) >0: possibly incorrect logic, we want to conpute distance from rv if no distance given in input file
        if len(vs)>0 and d==0:
            print(vs)
            med = statistics.median(vs)
            dnear, dfar = dist(l, b, med)
            if 'nan' not in str(dnear) and 'nan' not in str(dfar):
                a2 = (dnear+dfar)/2
            elif 'nan' in str(dnear) and 'nan' not in str(dfar):
                a2 = dfar
            elif 'nan' in str(dfar) and 'nan' not in str(dnear):
                a2 = dnear
            else:
                a2 = -999
        if len(av)>0 and str(a2)!='-999':
            count = 0
            a=0
            for lval in av:
                if lval >0:
                    count+=1
                    a+=float(lval)
            if count ==0:
                a = 0
            else:
                a = a / count
            t_eff,iofl = tef(a,a2)

        else:
            a = '-999'
            t_eff = '-999'
            iofl = '-999'


        #print(a2)
        #print(av2)
        r.append(a)
        r.append(a2)
        r.append(t_eff)
        r.append(iofl)
        lums.append(r)
    elif distfg == 1:
        #a = '-999'
        dist_cm = dist_kpc * 3.086e+21
        luminosity = 4*np.pi*dist_cm**2*finflux
        luminosity_sol = luminosity / 3.839e+33
        
        dist_cm_lower = dist_kpc_lower * 3.086e+21
        luminosity_lower = 4*np.pi*dist_cm_lower**2*finflux
        luminosity_sol_lower = luminosity_lower / 3.839e+33
        
        dist_cm_higher = dist_kpc_higher * 3.086e+21
        luminosity_higher = 4*np.pi*dist_cm_higher**2*finflux
        luminosity_sol_higher = luminosity_higher / 3.839e+33
        
        t_eff = '-999'
        iofl = '-999'
        
        a = luminosity_sol
        a_lower = luminosity_sol_lower
        a_higher = luminosity_sol_higher
        a2 = dist_kpc
        a2_lower = dist_kpc_lower
        a2_higher = dist_kpc_higher
        t_eff = '-999'
        iofl = '-999'
        r.append(a)
        r.append(a_lower)
        r.append(a_higher)
        r.append(a2)
        r.append(a2_lower)
        r.append(a2_higher)
        r.append(pm)
        r.append(pnra)
        r.append(e_pmra)
        r.append(pmde)
        r.append(e_pmde)
        r.append(excess_noise)
        r.append(excess_noise_sig)
        r.append(ipd_gof_ha)
        r.append(ruwe)
        r.append(t_eff)
        r.append(iofl)
        lums.append(r)

    print(r)
    ct = 0
    f = []
    n = []
    flagrad = 0
    for j in range(len(nu)):
        if nu[j] not in n:
            n.append(nu[j])
    for j in range(len(n)):
        if n[j]*10**(-9) > freqone and n[j]*10**(-9) < freqtwo:
            #print(obj)
            if n[j]*10**(-9) not in reps:
                reps.append(n[j]*10**(-9))
                reps.append(1)
            else:
                indx = reps.index(n[j]*10**(-9))
                reps[indx+1] +=1
    for j in range(len(nu)):

        if nu[j]*10**(-9) > freqone and nu[j]*10**(-9) < freqtwo:
            flagrad =1

            if str(nu[j]*10**(-9))[:6] not in freqlist:
                freqlist.append(str(nu[j]*10**(-9))[:6])
            ct+=1
            rw = []
            rw.append(nu[j]*10**(-9))


            if str(finflux) != '-999':
                ratio = nufnu[j]/finflux
            else:
                ratio = -999
            if str(finflux)!='-999' and len(av)>0:
                exrat = exprat(a2,nu[j]*10**(-9),a)
            else:
                exrat = '-999'

            nfnmjy = (nufnu[j]/(nu[j]*10**(-9)))*10**(17)

            print(obj, finflux, av, exrat, ratio)
            rw.append(ratio)
            rw.append(exrat)
            rw.append(nfnmjy)


            ratobsexp = -999
            if str(ratio)!='-999':
                con2+=1
            if str(ratio)!='-999' and str(exrat)!='-999':
                con+=1
                ratobsexp = float(ratio)/float(exrat)
            rw.append(ratobsexp)
            if ratobsexp < 0.1 and ratobsexp>0:
                sscand.append(obj)
                print(obj,nu[j])
            else:
                notsscand.append(obj)
            f.append(rw)
        else:
            notsscand.append(obj)
    if ct >  hlen:
        hlen = ct
    f.sort(key = sortSecond, reverse = True)
    rats.append(f)
    if flagrad ==1:
        lrad.append(obj)
high = 0
fhigh = 0
for i in range(len(reps)):
    if i%2==1 and reps[i]>high:
        high = reps[i]
        fhigh = reps[i-1]


l2 = []
for i in range(len(freqlist)):
    l2.append(float(freqlist[i]))
l2.sort(reverse=True)
freqlist=[]
for i in range(len(l2)):
    freqlist.append(str(l2[i]))
ssc = []
nssc = []
for i in sscand:
    if i not in ssc:
        ssc.append(i)
for i in notsscand:
    if i not in ssc and i not in nssc:
        nssc.append(i)

if distfg == 0:
    h = ['source object', 'type','ra(deg)', 'dec(deg)', 'l(deg)', 'b(deg)', 'bolometric flux (erg/cm^2/s)']
    for i in range(len(fils)):
        fin = fils[i].find('(')
        f = fils[i][:fin]
        fin2 = fils[i].find('_')
        f2 = fils[i][:fin2]
        h.append(f2+ '_catalogue')
        h.append(f + '(km/s, vlsr)')
        h.append(f2+'_dnear' +  '(kpc)')
        h.append(f2+ '_dfar' +  '(kpc)')
        h.append(f2+'_Lumnear' +  '(Lsun)')
        h.append(f2+'_Lumfar' +  '(Lsun)')
    h.append('avg_Lum(Lsun)')
    h.append('avg_dist(kpc)')
    h.append('Teff(K)')
    h.append("Ionizing_flux(1e48/s)")
    for i in range(len(freqlist)):
        h.append(freqlist[i]+"_freq(GHz)")
        h.append(freqlist[i]+"_nuFnu/Fbol_observed")
        h.append(freqlist[i]+"_nuFnu/Fbol_expected")
        h.append(freqlist[i]+"_nuFnu(mjy)")
        h.append(freqlist[i]+"_nuFnu/Fbol_obs/exp")

    rows = []
    for i in range(len(lums)):
        r = []
        r.append(lums[i][0])
        r.append(lums[i][1])
        r.append(lums[i][2])
        r.append(lums[i][3])
        r.append(lums[i][4])
        r.append(lums[i][5])
        r.append(lums[i][6])
        for j in range(len(fils)):
            flag = 0
            for k in range(len(lums[i])):
                if '(km/s)' in str(lums[i][k]):
                    if fils[j] in lums[i][k]:
                        if str(lums[i][k])[0] == '-':
                            #print('y')
                            ind = str(lums[i][k]).find('-', 2)
                        else:
                            ind = str(lums[i][k]).find('-', 1)

                        r.append(lums[i][k-1])
                        r.append(lums[i][k][:ind])
                        r.append(lums[i][k+1])
                        r.append(lums[i][k+2])
                        r.append(lums[i][k+3])
                        r.append(lums[i][k+4])
                        flag = 1
            if flag ==0:
                r.append(-999)
                r.append(-999)
                r.append(-999)
                r.append(-999)
                r.append(-999)
                r.append(-999)
        r.append(lums[i][-4])
        r.append(lums[i][-3])
        r.append(lums[i][-2])
        r.append(lums[i][-1])
        rows.append(r)

    for i in range(len(rats)):
        for p in range(len(freqlist)):
            flg = 0
            for j in range(len(rats[i])):
                for k in range(len(rats[i][j])):
                    print(len(rats[i][j]),k)
                    if freqlist[p] == str(rats[i][j][k])[:6] and k ==0:
                        ro = rats[i][j][k]
                        ro2 = rats[i][j][k+1]
                        ro3 = rats[i][j][k+2]
                        ro4 = rats[i][j][k+3]
                        #if len(rats[i][j]>k+3):
                        ro5 = rats[i][j][k+4]
                        flg = 1

            if flg ==1:
                rows[i].append(ro)
                rows[i].append(ro2)
                rows[i].append(ro3)
                rows[i].append(ro4)
                rows[i].append(ro5)

            else:
                rows[i].append(-999)
                rows[i].append(-999)
                rows[i].append(-999)
                rows[i].append(-999)
                rows[i].append(-999)
        while len(rows[i]) < (hlen*5+len(fils)*6+6):
            rows[i].append(-999)
elif distfg == 1:
    h = ['source object', 'type','ra(deg)', 'dec(deg)', 'l(deg)', 'b(deg)', 'bolometric flux (erg/cm^2/s)', 'avg_Lum(Lsun)', 'avg_Lum_lower(Lsun)', 'avg_Lum_upper(Lsun)','avg_dist(kpc)','avg_dist_lower(kpc)','avg_dist_upper(kpc)', 'PM(mas/yr)', 'pmRA(mas/yr)', 'e_pmRA(mas/yr)', 'pmDE(mas/yr)', 'e_pmDE(mas/yr)', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'ipd_gof_harmonic_amplitude', 'RUWE', 'Teff(K)', "Ionizing_flux(1e48/s)"]

with open(dir + "/" + root+"_dist_lum.csv", 'w') as l:
    l_writer = csv.writer(l, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    l_writer.writerow(h)
    if distfg == 0:
        for i in range(len(rows)):
            l_writer.writerow(rows[i])
    elif distfg == 1:
        for i in range(len(lums)):
            l_writer.writerow(lums[i])
l.close()

with open(dir + "/" + root+"_swollen_star_candidates.txt", 'w') as l:
    l.write("Swollen star candidates:\n")
    for i in ssc:
        l.write(i+'\n')
    l.write("\nOther:\n")
    for i in nssc:
        l.write(i+'\n')
l.close()

with open(dir + "/" + root+"_obj_radcontinuum.txt", 'w') as l:
    for i in lrad:
        l.write(i+'\n')
l.close()

rthd = dir + "/" + root + "_hdr_column_distlum.csv"
with open(rthd, 'w') as he:
    h_writer = csv.writer(he, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for i in range(len(h)):
        h_writer.writerow([h[i], " " + str(i+1)])

print(str(len(ob)) + ' objects done in ' "--- %s seconds ---" % (time.time() - start_time))
print(con)
print(con2)
