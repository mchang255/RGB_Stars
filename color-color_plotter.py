#python3
#This program plots the color-color plots and k-means clustering plots of the objects in the root directory.
#need root_colors.inp, root_coords.csv, and root_phot.csv on computer in same folder as this program
#run as follows:
#python color-color_plotter.py root_colors.inp
#ex .inp file: kastner08_colors.inp
#this output is saved in root_colors directory
#Major history upgrades
#7/15/21: Added the errors of the medians to the plottable .csv file. Also, created plots that would take into the account the medians.
#7/13/21: I edited the portion on retrieving values to collect all values for a certain filter and then take the median. I added near the end (section called TAKING MEDIAN INTO ACCOUNT) for the .csv files to be outputted.
#7/12/21: Added a section where magnitudes and the names of stars (that ARE NOT plottable) are outputted to a .csv file. There were some values that were infinity, so I modified code to add those stars with infinity values to the reject list.
#07/09/21: For 3D plots, same issue as 2D plots: Copied and pasted a section from k-means to the section creating the main color plot. Originally, main color plot outputted blank plots. Not the case anymore.
#07/08/21: Added a section where the magnitudes of the colors and the names of stars (that are plottable) are outputted to a .csv file.
#07/07/21: Copied and pasted a section from k-means to the section creating the main color plot. Originally, main color plot outputted blank plots. Not the case anymore.
#09/01/20: User input to rotate 3d plots and save
#08/21/20: Colors are now user inputs
#08/03/20: Separated SED plotter and color-color/kmeans plotter into different programs. SED program is named phot_plotter.py. Added 3D color-color and kmeans plots.
#IMPORTANT NOTES from a previous student Rabia, 7/12/21:
# 1) The color-color program will reject a source if it is missing any of the four fluxes
# necessary to plot the source. For example if you are plotting J-K vs H-A, if a source is
# missing A, then it will not be plotted.

# 2) if a list of objects had object type info from the source website (as opposed to just
# querying Simbad and taking the type listed there), and this info is the 3rd col in the
# coords file, then if the spfg is 1 it will take the type from that 3rd column in the
# coords file instead of the phot file column.

# 3) When there are multiple instances of a flux for a certain filter, it may be taking the
# first instance of that

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib  import cm
from matplotlib import colors as mcolors
import math
import csv
import sys
import os
from pathlib import Path
import time
import random
import matplotlib.lines as mlines
import mpld3
import astropy.units as u, astropy.constants as c
from mpl_toolkits import mplot3d
import pandas as pd
from shutil import copy
from mpl_toolkits.mplot3d import Axes3D
from astropy.stats import median_absolute_deviation

start_time = time.time()

root = sys.argv[1][:-4]
count = 0
count2 = 0
count3 = 0
count4 = 0
formats = []
nclus = 0
coordfile = []
objfg = []
rootobj = []
ccplots = []
azim = []
elev = []
with open(sys.argv[1], 'r') as inp:
    lines = inp.readlines()
    for line in lines:
        if line[:1] == '#':
            continue
        else:
            index = line.find('=')
            if 'dirnam' in line[:index]:
                rootobj.append(str(line[index+1:].strip()))
            if 'spfg' in line[:index]:
                objfg.append(str(line[index+1:].strip()))
            if 'coordinates' in line[:index]:
                coordfile.append(str(line[index+1:].strip()))
            if '.' in line[:index]:
                if line[index+1:].strip() == 'y':
                    formats.append(str(line[:index].strip()))
            if 'nclusters' in line[:index]:
                nclus = int(line[index+1:].strip())
            if 'azim' in line[:index]:
                azim.append(float(line[index+1:].strip()))
            if 'elev' in line[:index]:
                elev.append(float(line[index+1:].strip()))
            if 'ccax' in line[:index]: #which colors to plot as read from the _colors.inp file
                right = str(line[index+1:].strip())
                n = right.count(',')
                ind = 0
                ccplots.append([])
                for i in range(n+1):
                    ind2 = right.find('-',ind)
                    ind3 = right.find(',', ind+1)
                    if i == n:
                        l = [right[ind:ind2],right[ind2+1:]]
                    else:
                        l = [right[ind:ind2],right[ind2+1:ind3]]
                    ind = ind3+1
                    ccplots[-1].append(l)

print(ccplots)

path = str(Path(__file__).parent.absolute())
dir = path + '/' + str(root)
d = path + '/'

totlen = 0
types = []
for i in range(len(objfg)):
    if objfg[i] == '1':
        types.append(['-'])
        with open(coordfile[i], 'r') as c:
            lines = c.readlines()
            for line in lines:

                if line[:1] == "#":
                    continue
                else:
                    xx = line.split(",")
                    print(line)
                    types[i].append(str(xx[3].strip()))
    else:
        types.append([])

#print(result)
# radconobj = []
# for i in range(len(rootobj)):
#     with open(d + str(rootobj[i])+ '-2.5asec'+"/" + rootobj[i]+ '-2.5asec'+"_obj_radcontinuum.txt") as rad:

#         lines=rad.readlines()
#         for line in lines:
#             radconobj.append(line.strip())
info = []
for i in range(len(coordfile)):
    info.append([])
    with open(d + str(rootobj[i])+ '-2.5asec'+"/" + rootobj[i] + '-2.5asec'+"_phot.csv", mode='r') as phot:
        phot_reader = csv.reader(phot, delimiter=',')
        #next(phot_reader)
        for row in phot_reader:
            info[i].append(row)
axes = []
rejects = []
cccol1 = []
lab = []
for i in range(len(ccplots)):
    axes.append([])
    rejects.append([])
    cccol1.append([])
    lab.append([])
    for j in range(len(ccplots[i])):
        axes[i].append([])
        rejects[i].append([])

filter_list = []
for x in range(len(ccplots)):
    for y in range(len(ccplots[x])):
        for z in range(len(ccplots[x][y])):
            filter_list.append(ccplots[x][y][z])

filter_list = list(dict.fromkeys(filter_list))
#print(filter_list)

        
output = []
typelist = []
lablist = []
obs = []
medians = []
medians_error = []
for l in range(len(info)):
    obs.append([])
    for i in range(1,len(info[l])):
        obs[-1].append(info[l][i][0].strip())
        obj = info[l][i][0].replace(" ", "_")
        ra = info[l][i][1]
        dec = info[l][i][2]
        if objfg[l] == '1':
            type = types[l][i]
        else:
            type = info[l][i][3]
        wavelength = []
        fnu = []
        errorfnu = []
        filter = []
        for k in range(len(info[l][i])):
            if str(info[l][0][k][-18:]) == "wavelength(micron)" and "Lin:" not in str(info[l][0][k]) and "Sp:" not in str(info[l][0][k]):
                if float(info[l][i][k]) != -999:
                    wavelength.append(float(info[l][i][k]))
                    fnu.append(float(info[l][i][k+1]))
                    filter.append(info[l][0][k][:-19])
                    if info[l][i][k+2] == "--":
                        errorfnu.append(0)
                    else:
                        errorfnu.append(float(info[l][i][k+2]))

        nufnu = []
        errornufnu = []
        nu = []
        for p in range(len(fnu)):
            nu.append((2.99*10**(14))/wavelength[p])
            nufnu.append(nu[p]*fnu[p]*10**(-23))
            errornufnu.append(nu[p]*errorfnu[p]*10**(-23))
        #color-color plots
        totlen += 1
        #first, convert flux into magnitudes
        flux = []
        mags = []
        miss = []
        fils = []
        #print(filter)
        for k in range(len(filter)):
            if "(" in filter[k]:
                index = filter[k].find("(")
                fil = filter[k][:index]
                if ":" in fil:
                    index2 = fil.find(":")
                    fils.append(fil[index2+1:])
                    flux.append(fnu[k])
                    miss.append(fil)
            else:
                if ":" in filter[k]:
                    index2 = filter[k].find(":")
                    fils.append(filter[k][index2+1:])
                    miss.append(filter[k])
                    flux.append(fnu[k])
        fils2 = []
        miss2 = []
        out = 0
        def mag(f0, flux): #calculating magnitude from flux
            coef = -1/0.4
            log = np.log10(flux/f0)
            return coef*log
        obj = {}
        for i in filter_list:
            obj[i] = []
        median_star = []
        median_error = []
        for k in range(len(miss)):
            #print(miss[k])
            if "MSX" in miss[k]:
                if "A" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    out = flux[k]
                    mags.append(mag(55.81,flux[k]))
                    #print(str(mag(55.81,flux[k]))+',second round')
                    if flux[k] != -999:
                        if fils[k] in obj.keys(): #if band is part of input
                            obj[fils[k]].append(mag(55.81,flux[k])) #append magnitude to list
                if "C" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(26.432,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(26.432,flux[k]))
                if "D" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(18.267,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(18.267,flux[k]))
                if "E" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(8.666,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(8.666,flux[k]))
            if "2MASS" in miss[k]:
                if "J" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(1594,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(1594,flux[k]))
                if "H" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(1024,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(1024,flux[k]))
                if "K" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(666.7,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(666.7,flux[k]))
            if "Johnson" in miss[k]:
                if "J" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(1636.6,flux[k]))
                    #print(str(mag(1636.6,flux[k]))+'Johnson')
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(1636.6,flux[k]))
                if "H" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(1049.5,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(1049.5,flux[k]))
                if "K" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(653.2,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(653.2,flux[k]))
            if "AKARI" in miss[k]:
                if "9W" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(56.26,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(56.26,flux[k]))
                if "18W" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(12.00,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(12.00,flux[k]))
            if "WISE" in miss[k]:
                if "W1" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(306.68,flux[k]))
                    #print(mag(306.68,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(306.68,flux[k]))
                if "W2" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(170.66,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(170.66,flux[k]))
                if "W3" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(29.045,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(29.045,flux[k]))
                if "W4" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(8.2839,flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(8.2839,flux[k]))
            if "GALEX" in miss[k]:
                if "NUV" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(3.365*10**(-5),flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(3.365*10**(-5),flux[k]))
                if "FUV" in fils[k]:
                    miss2.append(miss[k])
                    fils2.append(fils[k])
                    mags.append(mag(1.079*10**(-4),flux[k]))
                    if flux[k] != -999:
                        if fils[k] in obj.keys():
                            obj[fils[k]].append(mag(1.079*10**(-4),flux[k]))
                

        for i in filter_list:
            median_star.append(np.median(obj[i]))
            median_error.append(median_absolute_deviation(obj[i])) #calculating error for each band
        medians.append(median_star)
        medians_error.append(median_error)

        #cccol = 'str'
        label = 'str'
        if 'GMV' in type:
            cccol='black'
            label='GMV'
        elif 'RSG' in type:
            cccol='dodgerblue'
            label= 'RSG'
        elif 'H II' in type or "HII" in type:
            cccol='forestgreen'
            label='HII region'
        elif 'O AGB' in type:
            cccol='darkblue'
            label='O AGB'
        elif 'C AGB' in type:
            cccol='red'
            label='C AGB'
        else:
            #cccol='blue'
            if '?' in type:
                place = type.index('?')
                label=type[:place]
            else:
                label = type
                cccol = 'red'
        if label not in lablist:
            lablist.append(label)
        if type not in typelist:
            typelist.append(type)
        for g in range(len(ccplots)): #the data to plot for each color is being extracted
            placeh = []
            for p in range(len(ccplots[g])):

                for f in range(len(ccplots[g][p])):
                    placeh.append(-99) #if a star does not have a value for a certain filter, -99 will be in place instead

                    for m in range(len(fils2)):
                        #print (ccplots[g][p][f])
                        if ccplots[g][p][f] =='K' or ccplots[g][p][f] =='J' or ccplots[g][p][f]=='H':
                            if ccplots[g][p][f] in fils2[m] and '2MASS' in miss2[m]:
                                placeh[-1] = mags[m]
                        else:
                            if ccplots[g][p][f] in fils2[m]:
                                placeh[-1] = mags[m]
            flg = 0
            #print(placeh)
            #print(mags)
            for t in range(len(placeh)):
                if placeh[t] ==-99:
                    flg = 1
            if flg ==0:
                if len(ccplots[g]) == 2: #if 2D
                    if (np.isfinite(placeh[0]-placeh[1]) & np.isfinite(placeh[2]-placeh[3])):
                        axes[g][0].append(placeh[0]-placeh[1]) #x-axis color
                        axes[g][1].append(placeh[2]-placeh[3]) #y-axis color
                        axes[g].append([])
                        axes[g][2].append(obs[-1][-1]) #name of star
                        cccol1[g].append(cccol)
                        lab[g].append(label)
                    else:
                        rejects[g][0].append(placeh[0]) #first mag
                        rejects[g][1].append(placeh[1]) #second mag
                        rejects[g].append([])
                        if ((ccplots[g][0][0] == ccplots[g][1][0])|(ccplots[g][0][1] == ccplots[g][1][0])): #ex: if I have J-H,H-K. There is an H in the two middle columns
                            rejects[g][2].append(placeh[3]) #third mag
                            rejects[g].append([])
                            rejects[g][3].append(obs[-1][-1]) #name
                        elif ((ccplots[g][0][0] == ccplots[g][1][1])|(ccplots[g][0][1] == ccplots[g][1][1])): #ex: if I have H-K,J-H. There is an H in the first and fourth column
                            rejects[g][2].append(placeh[2]) #third mag
                            rejects[g].append([])
                            rejects[g][3].append(obs[-1][-1]) #name
                        else:
                            rejects[g][2].append(placeh[2])
                            rejects[g].append([])
                            rejects[g][3].append(placeh[3])
                            rejects[g].append([])
                            rejects[g][4].append(obs[-1][-1]) #name
                elif len(ccplots[g]) == 3: #if 3D


                    axes[g][0].append(placeh[0]-placeh[1]) #x-axis color
                    axes[g][1].append(placeh[2]-placeh[3]) #y-axis color
                    axes[g][2].append(placeh[4]-placeh[5]) #z-axis color
                    axes[g].append([])
                    axes[g][3].append(obs[-1][-1]) #name of star
                    cccol1[g].append(cccol)
                    lab[g].append(label)
            else: #the stars are not plottable are being compiled
                if len(ccplots[g]) == 2: #if 2D
                    rejects[g][0].append(placeh[0]) #first mag
                    rejects[g][1].append(placeh[1]) #second mag
                    rejects[g].append([])
                    if ((ccplots[g][0][0] == ccplots[g][1][0])|(ccplots[g][0][1] == ccplots[g][1][0])): #ex: if I have J-H,H-K. There is an H in the two middle columns
                        rejects[g][2].append(placeh[3]) #third mag
                        rejects[g].append([])
                        rejects[g][3].append(obs[-1][-1]) #name
                    elif ((ccplots[g][0][0] == ccplots[g][1][1])|(ccplots[g][0][1] == ccplots[g][1][1])): #ex: if I have H-K,J-H. There is an H in the first and fourth column
                        rejects[g][2].append(placeh[2]) #third mag
                        rejects[g].append([])
                        rejects[g][3].append(obs[-1][-1]) #name
                    else:
                        rejects[g][2].append(placeh[2])
                        rejects[g].append([])
                        rejects[g][3].append(placeh[3])
                        rejects[g].append([])
                        rejects[g][4].append(obs[-1][-1]) #name


os.makedirs(os.path.dirname(dir + '/' + root+'_J-KvsK-A.pdf'), exist_ok=True)
copy(sys.argv[1], dir)
def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)
output = []
for j in range(len(axes)):
    
    if len(ccplots[j]) == 2:

        output.append('#centroid,'+ccplots[j][0][0]+'-'+ccplots[j][0][1]+','+ccplots[j][1][0]+'-'+ccplots[j][1][1]+'\n')
        X = []
        for i in range(len(axes[j][0])):
            new = []
            new.append(axes[j][0][i])
            new.append(axes[j][1][i])
            X.append(new)
        X = np.array(X)
        import seaborn as sns; sns.set()  # for plot styling
        from sklearn.cluster import KMeans
        kmeans2 = KMeans(n_clusters=nclus)
        #print(X)
        Y = []
        for k in range(len(X)):
            #print(X[k][0])
            if np.isfinite(float(X[k][0])) and np.isfinite(float(X[k][1])):
                Y.append(X[k])
        X=np.array(Y)
        
    if len(ccplots[j]) == 2:

        lablist = []
        for i in lab[j]:
            if i not in lablist:
                lablist.append(i)
        n=len(lablist)


        fig = plt.figure()
        ax = plt.axes()

        """cmap = plt.get_cmap('RdYlBu')

        colors = cmap(np.linspace(0, 1, n))
        #n = 0
        map = []
        for i in range(len(lablist)):
            map.append([lablist[i],colors[i]])
        col = []
        for i in lab[j]:
            for k in map:
                if k[0]==i:
                    col.append(k[1])"""


        plt.xlabel(ccplots[j][0][0]+'-'+ccplots[j][0][1])
        plt.ylabel(ccplots[j][1][0]+'-'+ccplots[j][1][1])
        #cm = plt.get_cmap('gist_rainbow')
        #ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
        for i in range(len(axes[j][0])):
            #if axes[j][2][i] not in radconobj:
            if lab[j][i]=='Evolved star' or lab[j][i] == 'HII region':

                ax.plot(axes[j][0][i],axes[j][1][i], '.', label = lab[j][i],color=cccol1[j][i])#color = col[i])
        handles, labels = ax.get_legend_handles_labels()
        handle_list, label_list = [], []
        for handle, label in zip(handles, labels):
            if label not in label_list:
                handle_list.append(handle)
                label_list.append(label)
        lgd = plt.legend(handle_list, label_list, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}, ncol = 1)
        #for i in range(len(axes[j][0])):
        #    if axes[j][2][i] in radconobj:
        #        if lab[j][i]=='HII region' or lab[j][i] == 'Evolved star':
        #            ax.plot(axes[j][0][i],axes[j][1][i], 's', label = lab[j][i], color=col[i])

        if ccplots[j][1][0]=='J' and ccplots[j][1][1]=='K' and ccplots[j][0][0]=='K':
            plt.xlim(0,11)
            plt.ylim(-1,8)
            #RSG
            ax.plot([0.75,2.7,2.7,0.75, 0.75],[0.7,0.7,1.75,1.75, 0.7],',-', color='dodgerblue', linewidth = 1)
            #HII
            ax.plot([6.25,9.0,9,6.25, 6.25],[0,0,2.5,2.5, 0],',-', color='forestgreen', linewidth = 1)
            #GMV
            ax.plot([0.5,2.5,2.5,0.5, 0.5],[1.1,1.6,1.9,1.4, 1.1],',-', color='black', linewidth = 1)
            #O AGB
            ax.plot([1.95,3.78,5.5,5.5, 1.95],[2,2,3.2,6.4, 2],',-', color='darkblue', linewidth = 1)
            #C AGB
            ax.plot([2.39,8.702,9.505,6.2067, 2.39],[2.4,2.4,3.5,6.5, 2.4],',-', color='red', linewidth = 1)
        elif ccplots[j][1][0]=='H' and ccplots[j][1][1]=='K' and ccplots[j][0][0]=='K':
            plt.xlim(0,11)
            plt.ylim(-1,4)
            #RSG
            ax.plot([0.75,2.7,2.7,0.75, 0.75],[0.1,0.1,0.7,0.7, 0.1],',-', color='dodgerblue', linewidth = 1)
            #HII
            ax.plot([6.25,9.0,9,6.25, 6.25],[-0.25,-0.25,1.5,1.5,-0.25],',-', color='forestgreen', linewidth = 1)
            #GMV
            ax.plot([0.5,2.5,2.5,0.5, 0.5],[0.2,0.6,0.9,0.5, 0.2],',-', color='black', linewidth = 1)
            #O AGB
            ax.plot([1.784,3.81,5.5,5.5, 1.784],[0.7,0.7,1.495,2.595, 0.7],',-', color='darkblue', linewidth = 1)
            #C AGB
            ax.plot([2.582,6.945,9.66,7.379, 2.582],[1.092,1.485,2.3,3.398, 1.092],',-', color='red', linewidth = 1)
        plt.scatter(X[:, 0], X[:, 1], s=20, cmap='viridis')
        plt.title(root)# + "(" +str(len(axes[j][0]))+'/' +str(totlen)+")" )
        plt.tight_layout()
        box = ax.get_position()
        xlab = ccplots[j][0][0]+'-'+ccplots[j][0][1]
        ylab = ccplots[j][1][0]+'-'+ccplots[j][1][1]
        ax.set_position([box.x0, box.y0, box.width , box.height])
        fig.set_size_inches(10,5)
        #plt.show()
        plt.savefig(dir + '/' + root+'_'+ylab+'vs'+xlab+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close()

    elif len(ccplots[j]) == 3:
        for y in range(len(azim)):
            lablist = []
            for i in lab[j]:
                if i not in lablist:
                    lablist.append(i)
            n=len(lablist)
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            cmap = plt.get_cmap('RdYlBu')

            colors = cmap(np.linspace(0, 1, n))
            map = []
            for i in range(len(lablist)):
                map.append([lablist[i],colors[i]])
            col = []
            for i in lab[j]:
                for k in map:
                    if k[0]==i:
                        col.append(k[1])
            output.append('#centroid,'+ccplots[j][0][0]+'-'+ccplots[j][0][1]+','+ccplots[j][1][0]+'-'+ccplots[j][1][1]+','+ccplots[j][2][0]+'-'+ccplots[j][2][1]+'\n')

            X = []
            for i in range(len(axes[j][0])): #getting data points to plot
                new = []
                new.append(axes[j][0][i])
                new.append(axes[j][1][i])
                new.append(axes[j][2][i])
                X.append(new)
            X = np.array(X)
            ax.set_xlabel(ccplots[j][0][0]+'-'+ccplots[j][0][1])
            ax.set_ylabel(ccplots[j][1][0]+'-'+ccplots[j][1][1])
            ax.set_zlabel(ccplots[j][2][0]+'-'+ccplots[j][2][1])
            ax.scatter(X[:, 0], X[:, 1], X[:, 2], s=20, cmap='viridis') #plotting data points
            for i in range(len(axes[j][0])):
                #if axes[j][3][i] not in radconobj:
                if lab[j][i]=='Evolved star' or lab[j][i] == 'HII region':
                    ax.scatter3D(axes[j][0][i], axes[j][1][i], axes[j][2][i], s=20,  label = lab[j][i], depthshade=True, color=cccol1[j][i])
            for handle, label in zip(handles, labels):
                if label not in label_list:
                    handle_list.append(handle)
                    label_list.append(label)
            lgd =plt.legend(handle_list, label_list, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}, ncol = 1)
            #for i in range(len(axes[j][0])):
            #    if axes[j][3][i] in radconobj:
            #        if lab[j][i]=='HII region' or lab[j][i] == 'Evolved star':
            #            ax.scatter3D(axes[j][0][i], axes[j][1][i], axes[j][2][i], s=20,  label = lab[j][i], depthshade=True, color = col[i])

            ax.set_title(root)# + "(" +str(len(axes[j][0]))+'/' +str(totlen)+")", pad=40)
            #plt.title(root + "(" + str(len(axes[j][0]))+'/' +str(totlen)+ ")")
            plt.tight_layout()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width , box.height])
            xlab = ccplots[j][0][0]+'-'+ccplots[j][0][1]
            ylab = ccplots[j][1][0]+'-'+ccplots[j][1][1]
            zlab = ccplots[j][2][0]+'-'+ccplots[j][2][1]
            ax.view_init(azim=azim[y], elev = elev[y])
            fig.set_size_inches(10,5)
            plt.savefig(dir + '/' + root+'_'+zlab+'vs'+ylab+'vs'+xlab+'_'+str(azim[y])+','+str(elev[y])+'.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
            #plt.show()
            plt.close()


#Kmeans
output = []
for j in range(len(axes)):

    if len(ccplots[j]) == 2:

        output.append('#centroid,'+ccplots[j][0][0]+'-'+ccplots[j][0][1]+','+ccplots[j][1][0]+'-'+ccplots[j][1][1]+'\n')
        X = []
        for i in range(len(axes[j][0])):
            new = []
            new.append(axes[j][0][i])
            new.append(axes[j][1][i])
            X.append(new)
        X = np.array(X)
        import seaborn as sns; sns.set()  # for plot styling
        from sklearn.cluster import KMeans
        kmeans2 = KMeans(n_clusters=nclus)
        #print(X)
        Y = []
        for k in range(len(X)):
            #print(X[k][0])
            if np.isfinite(float(X[k][0])) and np.isfinite(float(X[k][1])):
                Y.append(X[k])
        X=np.array(Y)
        #print(X) #this shows up in the terminal
        kmeans2.fit(X)
        y_kmeans = kmeans2.predict(X)

        fig = plt.figure()
        ax = plt.axes()
        plt.xlabel(ccplots[j][0][0]+'-'+ccplots[j][0][1])
        plt.ylabel(ccplots[j][1][0]+'-'+ccplots[j][1][1])
        if ccplots[j][1][0]=='J' and ccplots[j][1][1]=='K' and ccplots[j][0][0]=='K':
            plt.xlim(0,11)
            plt.ylim(-1,8)
            #RSG
            ax.plot([0.75,2.7,2.7,0.75, 0.75],[0.7,0.7,1.75,1.75, 0.7],',-', color='dodgerblue', linewidth = 1)
            #HII
            ax.plot([6.25,9.0,9,6.25, 6.25],[0,0,2.5,2.5, 0],',-', color='forestgreen', linewidth = 1)
            #GMV
            ax.plot([0.5,2.5,2.5,0.5, 0.5],[1.1,1.6,1.9,1.4, 1.1],',-', color='black', linewidth = 1)
            #O AGB
            ax.plot([1.95,3.78,5.5,5.5, 1.95],[2,2,3.2,6.4, 2],',-', color='darkblue', linewidth = 1)
            #C AGB
            ax.plot([2.39,8.702,9.505,6.2067, 2.39],[2.4,2.4,3.5,6.5, 2.4],',-', color='red', linewidth = 1)
        elif ccplots[j][1][0]=='H' and ccplots[j][1][1]=='K' and ccplots[j][0][0]=='K':
            plt.xlim(0,11)
            plt.ylim(-1,4)
            #RSG
            ax.plot([0.75,2.7,2.7,0.75, 0.75],[0.1,0.1,0.7,0.7, 0.1],',-', color='dodgerblue', linewidth = 1)
            #HII
            ax.plot([6.25,9.0,9,6.25, 6.25],[-0.25,-0.25,1.5,1.5,-0.25],',-', color='forestgreen', linewidth = 1)
            #GMV
            ax.plot([0.5,2.5,2.5,0.5, 0.5],[0.2,0.6,0.9,0.5, 0.2],',-', color='black', linewidth = 1)
            #O AGB
            ax.plot([1.784,3.81,5.5,5.5, 1.784],[0.7,0.7,1.495,2.595, 0.7],',-', color='darkblue', linewidth = 1)
            #C AGB
            ax.plot([2.582,6.945,9.66,7.379, 2.582],[1.092,1.485,2.3,3.398, 1.092],',-', color='red', linewidth = 1)
        centers = kmeans2.cluster_centers_
        sc = plt.scatter(X[:, 0], X[:, 1], c=y_kmeans, s=20, cmap='viridis')
        labs = []
        for i in range(len(y_kmeans)):
            ind = y_kmeans[i]
            labs.append(str(ind+1)+": (" + str(centers[ind,0])[:6] +',' + str(centers[ind,1])[:6] + ')')
        clset = set(zip(y_kmeans, labs))
        handles = [plt.plot([],color=sc.get_cmap()(sc.norm(c)),ls="", marker="o")[0] for c,l in clset ]
        labels = [l for c,l in clset]
        ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1, 0.5), prop={'size': 8})

        plt.scatter(centers[:, 0], centers[:, 1], c='black', s=100, alpha=0.5)
        for i in range(nclus):
            plt.text(centers[i,0], centers[i,1], i+1)
            output.append(str(i+1) + ',' + str(centers[i,0]) +',' + str(centers[i,1]) + '\n')
        plt.title(root + "(" +str(len(axes[j][0]))+'/' +str(totlen)+")" )
        #plt.legend( loc="lower center", bbox_to_anchor=(0.5, -0.3))
        plt.tight_layout()

        #box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width , box.height])
        plt.savefig(dir + '/' + root+'_'+ylab+'vs'+xlab+'_kmeans.pdf')
        plt.close()

    elif len(ccplots[j]) == 3:
        print('Y')

#3d kmeans
        print(azim)
        for y in range(len(azim)):
            print('YES')
            output.append('#centroid,'+ccplots[j][0][0]+'-'+ccplots[j][0][1]+','+ccplots[j][1][0]+'-'+ccplots[j][1][1]+','+ccplots[j][2][0]+'-'+ccplots[j][2][1]+'\n')

            X = []
            for i in range(len(axes[j][0])):
                new = []
                new.append(axes[j][0][i])
                new.append(axes[j][1][i])
                new.append(axes[j][2][i])
                X.append(new)
            X = np.array(X)
            from sklearn.cluster import KMeans
            kmeans2 = KMeans(n_clusters=nclus)
            Y = []
            for k in range(len(X)):
                #print(X[k][0])
                if np.isfinite(float(X[k][0])) and np.isfinite(float(X[k][1])):
                    Y.append(X[k])
            X = np.array(Y)
            print(X)
            kmeans2.fit(X)
            y_kmeans = kmeans2.predict(X)

            #matplotlib.use('WebAgg')

            fig = plt.figure()
            ax = fig.add_subplot(111, projection = '3d')
            ax.set_xlabel(ccplots[j][0][0]+'-'+ccplots[j][0][1])
            ax.set_ylabel(ccplots[j][1][0]+'-'+ccplots[j][1][1])
            ax.set_zlabel(ccplots[j][2][0]+'-'+ccplots[j][2][1])
            sc = ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=y_kmeans, s=20, cmap='viridis')
            centers = kmeans2.cluster_centers_

            labs = []
            for i in range(len(y_kmeans)):
                ind = y_kmeans[i]
                labs.append(str(ind+1)+": (" + str(centers[ind,0])[:6] +',' + str(centers[ind,1])[:6] + ',' + str(centers[ind,2])[:6] +')')
            clset = set(zip(y_kmeans, labs))
            handles = [plt.plot([],[],color=sc.get_cmap()(sc.norm(c)),ls="", marker="o")[0] for c,l in clset ]
            labels = [l for c,l in clset]
            ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1, 0.5), prop={'size': 8})

            ax.scatter(centers[:, 0], centers[:, 1], centers[:, 2], c='black', s=100, alpha=0.5)
            for i in range(nclus):
                ax.text(centers[i,0], centers[i,1], centers[i,2], i+1)
                output.append(str(i+1) + ',' + str(centers[i,0]) +',' + str(centers[i,1]) + ','+ str(centers[i,2])+ '\n')
            plt.title(root+ "(" +str(len(axes[j][0]))+'/' +str(totlen)+")" )
            plt.tight_layout()
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width, box.height])
            xlab = ccplots[j][0][0]+'-'+ccplots[j][0][1]
            ylab = ccplots[j][1][0]+'-'+ccplots[j][1][1]
            zlab = ccplots[j][2][0]+'-'+ccplots[j][2][1]
            print(zlab)
            ax.view_init(azim=azim[y], elev = elev[y])
            plt.savefig(dir + '/' + root+'_'+zlab+'vs'+ylab+'vs'+xlab+'_'+str(azim[y])+','+str(elev[y])+'_kmeans3d.pdf')
            #plt.show()
            plt.close()

    #higher dim kmeans
    elif len(ccplots[j]) >3:

        s = '#centroid'
        for i in range(len(ccplots[j])):
            s = s+','
            for k in range(len(ccplots[j][i])):
                if k !=0:
                    s=s+'-'
                s = s+ccplots[j][i][k]
        s=s+'\n'
        output.append(s)

        X = []
        for i in range(len(axes[j][0])):
            new = []
            for k in range(len(ccplots[j])):
                new.append(axes[j][k][i])
            X.append(new)
        X = np.array(X)
        from sklearn.cluster import KMeans
        kmeans2 = KMeans(n_clusters=nclus)
        kmeans2.fit(X)
        y_kmeans = kmeans2.predict(X)
        centers = kmeans2.cluster_centers_

        for i in range(nclus):
            s = str(i+1)+ ','
            for k in range(len(ccplots[j][i])):
                s=s+str(centers[i,0])
                if k !=len(ccplots[j][i])-1:
                    s=s+','
            s = s+'\n'
            output.append(s)

        st = dir + "/" + root + "_centroid_"
        for i in range(len(ccplots[j])):
            if i !=0:
                st = st+'vs'
            for k in range(len(ccplots[j][i])):
                st = st+ccplots[j][i][k]
                if k !=len(ccplots[j][i])-1:
                    st=st+'-'
        st = st+'.csv'

    if len(ccplots[j]) == 3:
        filnam = dir + "/" + root + "_centroid_"+zlab+'vs'+ylab+'vs'+xlab+".csv"
    elif len(ccplots[j]) ==2:
        filnam=dir + "/" + root + "_centroid_"+ylab+'vs'+xlab+".csv"
    elif len(ccplots[j]) > 3:
        filnam=st

    with open(filnam, 'w') as r:
        for i in range(len(output)):
            r.write(output[i])
    r.close()
    output = []
    
for j in range(len(ccplots)): #creating the list of stars (that have plottable colors) along with the magnitudes of colors
    first_color = axes[j][0]
    second_color = axes[j][1]
    name = axes[j][2]
    xlab = ccplots[j][0][0]+'-'+ccplots[j][0][1]
    ylab = ccplots[j][1][0]+'-'+ccplots[j][1][1]
    
    data = list(zip(name, first_color, second_color))
    f_out = (dir + '/' + root+'_mag_'+ylab+'vs'+xlab+'.csv')
    
    top = 'Name, ' +xlab+', '+ylab
    
    np.savetxt(f_out, data, fmt='%s', delimiter=",", header=top)

for j in range(len(ccplots)): #creating the list of stars are not plottable along with the magnitudes viewed in a specific band (ex: J, H, K)
    xlab = ccplots[j][0][0]+'-'+ccplots[j][0][1]
    ylab = ccplots[j][1][0]+'-'+ccplots[j][1][1]
    
    f_out = (dir + '/' + root+'_rejects_'+ylab+'vs'+xlab+'.csv')
    
    if ((ccplots[j][0][0] == ccplots[j][1][0])|(ccplots[j][0][1] == ccplots[j][1][0])):
        reject_mag1 = rejects[j][0]
        reject_mag2 = rejects[j][1]
        reject_mag3 = rejects[j][2]
        name = rejects[j][3]
        top = 'Name, ' +ccplots[j][0][0]+', '+ccplots[j][0][1]+', '+ccplots[j][1][1]
        data = list(zip(name, reject_mag1, reject_mag2, reject_mag3))
    elif ((ccplots[j][0][0] == ccplots[j][1][1])|(ccplots[j][0][1] == ccplots[j][1][1])):
        reject_mag1 = rejects[j][0]
        reject_mag2 = rejects[j][1]
        reject_mag3 = rejects[j][2]
        name = rejects[j][3]
        data = list(zip(name, reject_mag1, reject_mag2, reject_mag3))
        top = 'Name, ' +ccplots[j][0][0]+', '+ccplots[j][0][1]+', '+ccplots[j][1][0]
    else:
        reject_mag1 = rejects[j][0]
        reject_mag2 = rejects[j][1]
        reject_mag3 = rejects[j][2]
        reject_mag4 = rejects[j][3]
        name = rejects[j][4]
        data = list(zip(name, reject_mag1, reject_mag2, reject_mag3, reject_mag4))
        top = 'Name, ' +ccplots[j][0][0]+', '+ccplots[j][0][1]+', '+ccplots[j][1][0]+', '+ccplots[j][1][1]
    
    np.savetxt(f_out, data, fmt='%s', delimiter=",", header=top)

#TAKING MEDIAN INTO ACCOUNT
for i in range(len(medians)):
    for j in range(len(medians[i])):
        if np.isnan(medians[i][j]):
            medians[i][j] = -999 #replacing nan's in list with -999

master_first_color = []
master_second_color = []
for b in range(len(ccplots)): #outputting medians of plottable and non-plottable stars to .csv files
    for x in range(len(obs)):
        first_colors = []
        second_colors = []
        error_first_colors = []
        error_second_colors = []
        name = []
        rejects = []
        rejects_name = []
        top_rejects = 'Name'
        xlab = ccplots[b][0][0]+'-'+ccplots[b][0][1]
        ylab = ccplots[b][1][0]+'-'+ccplots[b][1][1]
        for y in range(len(obs[x])):
            temp_list = []
            if medians[y][filter_list.index(ccplots[b][0][0])] == -999 or medians[y][filter_list.index(ccplots[b][0][1])] == -999 or medians[y][filter_list.index(ccplots[b][1][0])] == -999 or medians[y][filter_list.index(ccplots[b][1][1])] == -999:
                temp_list.append(obs[x][y])
                for a in filter_list: #the rejects being appended to their lists
                    if a in ccplots[b][0] or a in ccplots[b][1]:
                        temp_list.append(medians[y][filter_list.index(a)])
                rejects.append(temp_list)
            else: #the plottables and their errors being appended to their lists
                first_color = medians[y][filter_list.index(ccplots[b][0][0])] - medians[y][filter_list.index(ccplots[b][0][1])]
                second_color = medians[y][filter_list.index(ccplots[b][1][0])] - medians[y][filter_list.index(ccplots[b][1][1])]
                error_first_color = medians_error[y][filter_list.index(ccplots[b][0][0])] - medians_error[y][filter_list.index(ccplots[b][0][1])]
                error_second_color = medians_error[y][filter_list.index(ccplots[b][1][0])] - medians_error[y][filter_list.index(ccplots[b][1][1])]
                
                first_colors.append(first_color)
                second_colors.append(second_color)
                error_first_colors.append(error_first_color)
                error_second_colors.append(error_second_color)
                name.append(obs[x][y])

        #outputting the rejects to their .csv file
        f_out_rejects = (dir + '/' + root+'_median_rejects_'+ylab+'vs'+xlab+'.csv')
        for s in filter_list:
            if s in ccplots[b][0] or s in ccplots[b][1]:
                top_rejects = top_rejects+','+s

        np.savetxt(f_out_rejects, rejects, fmt='%s', delimiter=",", header=top_rejects)

        #for plotting later
        master_first_color.append(first_colors)
        master_second_color.append(second_colors)
        
        #outputting the plottables to their .csv file
        data = list(zip(name, first_colors, error_first_colors, second_colors, error_second_colors))
        f_out = (dir + '/' + root+'_median_mag_'+ylab+'vs'+xlab+'.csv')
        top = 'Name,'+xlab+',Error (±) of '+xlab+','+ylab+',Error (±) of '+ylab
        np.savetxt(f_out, data, fmt='%s', delimiter=",", header=top)

for c in range(len(master_first_color)): #plotting the plottables
    plt.figure()
    xlab = ccplots[c][0][0]+'-'+ccplots[c][0][1]
    ylab = ccplots[c][1][0]+'-'+ccplots[c][1][1]
    plt.scatter(master_first_color[c],master_second_color[c])
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(ylab+'vs.'+xlab)
    plt.savefig(dir + '/' + root+'_median_'+ylab+'vs'+xlab+'.png', dpi = 150)
    
print(typelist)
print("plots created in --- %s seconds ---" % (time.time() - start_time))
