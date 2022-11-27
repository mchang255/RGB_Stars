#python3
#This program queries various online databases (Vizier, IRSA) and retreives photometry data for objects in a given input coordinate file. Then, it compiles this data into a dictionary-readable format.
#install requests package, PyAstronomy package
#have root_inp.txt, root_coords.csv, optical spectrum file (optional), and mol line files (5 total, optional) on computer in same folder as this program
#mol line files are named as follows: root_molline_data.csv where molline can be h2o, oh, sio, met1, met2
#optical spectra file is named as follows: root_optical.csv
#run as follows:
#python get-dict-phot.py root.inp
#ex .inp file: kastner08.inp
#this output is saved in root directory
#if wanting to query vizier sed tool, use code between keywords "begin query vizier sed" and "end query vizier sed" and modify the points in url to get the data in a 2D list. Then extract whatever information you need.
#if wanting to query vizier catalogue, use code between keywords "begin query vizier cat" and "end query vizier cat" and modify the catalogue in url to get the data in a 2D list. Then extract whatever information you need.
#if wanting to query irsa catalogue, use code between keywords "begin query irsa cat" and "end query irsa cat" and modify the catalogue in url to get the data in a 2D list. Then extract whatever information you need.
#if wanting to query mol lins,
    #1) add catalogue to for loop under comment "mollin #1" in same format as other entries
    #2) add catalogue to for loop under comment "mollin #2" with the cols that you want in the phot file (can be found by inspecting catalogue in vizier if theres velocity, peaks, etc)
    #3) add catalogue as an if statement under comment "mollin #3" following format as other catalogues, but make sure that you change the columns to whatever entires you wanted to add
#Major updates
#8/17/21: Changed .dat file to be raw5.csv, and it also contains nuFnu and nuFnu error. The final .dat file is what would be run into MoD. That is three columns: filter name and appropriate flux value and magnitude (along with three hashtag/pound signs at the bottom)
#8/9/21: Added section that converts flux (mJy) and flux error (mJy) to mag for filters that require it. This outputs to the final .dat file.
#8/6/21: Added section that takes median average of those with the same filter name (section '#taking median averages'). raw3 deletes identical duplicates between filter name, flux, flux error, and wavelength (this includes commented lines). raw4 takes median average of lines that have duplicates filters.
#8/3/21: Added a section 'SED-RAW' that would output the filter name, flux, flux error, wavelength, offset, and catalog name into a .dat for each star. raw1 is directly taking the filter names from _phot_order.csv. raw2 is substituting the filter names as specified from MoD.
#6/25/21: Made it so .inp file would be the first one read
#11/03/20: Fixed error reading certain VOTables
#09/02/20: Added UNWISE
#08/31/20: Added option of inputting sexagesimal coords and having a different separator in coord file, created file with phot ordered by wavelength
#08/12/20: Added .csv file that gives header names in col 1, column number in phot file in col2
#08/06/20: Added upper limit flags to IRAS filters when applicable
#08/05/20: Added IRAS II/126/sources catalogue from VizieR
#07/31/20: Added CATWISE catalogue from IRSA database
import urllib.request
import urllib.parse
import requests
import csv
from astropy.io.votable import parse
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip
from collections import Counter
import time
import sys
from pathlib import Path
import os
import math
from PyAstronomy import pyasl
import wget
import numpy as np
from shutil import copy
from subprocess import call
import weightedstats as ws

start_time = time.time()

#read in .csv file with object coordinates and create lists with names and coords of objects
obj = []
ra = []
dec = []
acc = []
tab2 = []
num = 0
mission = []
amac = []
amacmission = []

rad = 0
flag = 0
root = sys.argv[1][:-4]
coordfile = root + "_coords.csv"
inpfile = root + ".inp"
sp3flag = 0
sp2flag = 0
sp1flag = 0
sep = ','
form = 'd'
with open(inpfile, 'r') as inp:
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
            if 'radius' in line[:index]:
                rad = int(line[index+1:].strip())
            elif 'fluxdiff_werr' in line[:index]:
                fluxdiff_werr = float(line[index+1:].strip())
            elif 'fluxdiff_frac' in line[:index]:
                fluxdiff_frac = float(line[index+1:].strip())
            elif 'inp_cerr' in line[:index]:
                print('y')
                inp_cerr = float(line[index+1:].strip())
            elif 'xmlflag' in line[:index]:
                flag = int(line[index+1:].strip())
            elif 'amac' in line[:index]:
                index2 = line.find('_amac')
                amac.append(float(line[index+1:].strip()))
                amacmission.append(str(line[:index2].strip()))


path = str(Path(__file__).parent.absolute())
root = root+'-'+str(inp_cerr)+"asec"
dir = path + '/' + str(root)
#with open("coords_ohir.csv") as c:
with open(coordfile, 'r') as c:
    lines = c.readlines()
    for line in lines:

        if line[:1] == "#":
            continue
        else:
            x = line.split(sep)
            obj.append(str(x[0].strip()))
            if form == 'd':
                ra.append(float(x[1]))
                dec.append(float(x[2]))
            elif form =='s':
                print(x[1])
                c = SkyCoord(str(x[1])+' '+str(x[2]), unit=(u.hourangle, u.deg))
                ra.append(float(c.ra.degree))
                dec.append(float(c.dec.degree))
            elif form =='b':
                if x[3].strip() =='d':
                    ra.append(float(x[1]))
                    dec.append(float(x[2]))
                elif x[3].strip()=='h':
                    c = SkyCoord(str(x[1])+' '+str(x[2]), unit=(u.hourangle, u.deg))
                    ra.append(float(c.ra.degree))
                    dec.append(float(c.dec.degree))

counter = 1
downloading_phot_time = 0
for i in range(len(obj)):
    if obj[i] == "":
        obj[i] = "noname" + str(counter)
        counter+=1

noheader = dir + "/" + root + "_phot_noheader.csv"
rejects = dir + '/' + root + "_rejects.txt"
os.makedirs(os.path.dirname(noheader), exist_ok=True)
subdir = dir + '/xmlfiles'
if flag == '1':
    os.makedirs(os.path.dirname(subdir+"/files.txt"), exist_ok=True)
photometry = open(noheader, "w")
with open(rejects,"w") as r:
    r.write('input coord error is: ' + str(inp_cerr) + ' arcsec \n')
r.close()
for i in range(len(ra)):
    lins = []
    headers = {"User-Agent": "Mozilla/5.0 (Windows NT 6.1; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/72.0.3626.119 Safari/537.36"}
    points = str(ra[i]) + ',' + str(dec[i])
    print('Querying Vizier')
#VIZIER
    print(obj[i])
    #get photometry of object from vizier
    #begin query vizier sed
    url = 'https://vizier.u-strasbg.fr/viz-bin/sed?-c=' + urllib.parse.quote(points) + '&-c.rs=' + str(rad)
    print(url)
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)

    #put photometry in a votable and extract info
    if flag == '1':
        xmlfile = subdir + "/" + str(obj[i]) + "_vizier_phot.xml"
    else:
        xmlfile = dir + "/" + "phot.xml"
    p = open(xmlfile, "wb")
    p.write(req.content)
    exfg = 0
    with open(xmlfile,'r') as x:
        lines = x.readlines()
        for line in lines:
            if 'COOSYS' in line:
                exfg =1
    x.close()

    if exfg ==1:
        votable = parse(xmlfile)
        table = votable.get_first_table()
        data = table.array
        p.close()
    else:
        data = []
        p.close()
    #end query vizier sed
    #put individual information in lists
    fGhz = []
    fluxJy = []
    efluxJy = []
    filter = []
    righ = []
    dela = []
    tab = []

    ind = []
    newnam = []


    #get rid of duplicates too
    for j in range(len(data)):
        print(data[j])
        x = str(data[j]).split(", ")
        print(x)
        if len(x)!= 10:
            riga, decl, tabn, i_d, freq, fl, eflux, _filter = str(data[j]).split(", ")
        else:
            riga = x[0]
            decl=x[1]
            tabn=x[2]
            i_d=x[3]
            freq = x[6]
            fl = x[7]
            eflux = x[8]
            _filter = x[9]
        if fl !='--':
            if str(_filter)[1:-2] not in filter:
                righ.append(float(riga[1:]))
                dela.append(float(decl))
                fGhz.append(float(freq))
                fluxJy.append(float(fl))
                tab.append(str(tabn)[2:-1])
                if eflux == '--':
                    efluxJy.append(str(eflux))
                else:
                    efluxJy.append(float(eflux))
                filter.append(str(_filter)[1:-2])
            else:
                indices = [m for m, x in enumerate(filter) if x == str(_filter)[2:-2]]
                flg = []
                wh = 0
                for h in range(len(indices)):
                    if str(eflux) != '--' and str(efluxJy[indices[h]]) != '--':
                        if abs(float(fl)-float(fluxJy[indices[h]])) >= fluxdiff_werr*(float(eflux)**2+float(efluxJy[indices[h]])**2)**(1/2):
                            wh = 1
                        else:
                            flg.append(1)
                    else:
                        ave = (float(fl)+float(fluxJy[indices[h]]))/2
                        if ave != 0 and abs(float(fl)-float(fluxJy[indices[h]]))/ave >= fluxdiff_frac:
                            wh = 1
                        else:
                            flg.append(1)
                if len(flg) == len(indices):
                    righ.append(float(riga[1:]))
                    dela.append(float(decl))
                    fGhz.append(float(freq))
                    fluxJy.append(float(fl))
                    tab.append(str(tabn)[2:-1])
                    if eflux == '--':
                        efluxJy.append(str(eflux))
                    else:
                        efluxJy.append(float(eflux))
                    filter.append(str(_filter)[2:-2])
                elif efluxJy[filter.index(str(_filter)[2:-2])] != '--':
                    continue
                elif efluxJy[filter.index(str(_filter)[2:-2])] == '--' and eflux != '--':
                    efluxJy[filter.index(str(_filter)[2:-2])] = eflux
                elif efluxJy[filter.index(str(_filter)[2:-2])] == '--' and eflux == '--':
                    continue

    #fill wavelengths list from frequency
    wmicron = []
    for l in range(len(fGhz)):
        wmicron.append(299792.458/fGhz[l])

    #SPLASH
    #begin query vizier cat
    url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=J/ApJS/247/5/table1&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    p = open(xmlfile, "wb")
    p.write(req.content)
    votable = parse(xmlfile)
    table = votable.get_first_table()
    data = table.array
    p.close()
    #end query vizier cat

    if len(data)>0:
        sp3flag = 1
        x = str(data[0]).split(", ")

        hmsra = x[4][2:-1]
        hmsdec = x[5][2:-1]
        c = SkyCoord(str(hmsra)+' '+str(hmsdec), unit=(u.hourangle, u.deg))
        righ.append(float(c.ra.degree))
        dela.append(float(c.dec.degree))
        fluxJy.append(float(x[6]))
        freq = float(x[3][-6:-2])
        fGhz.append(float(freq))
        wmicron.append(299792458/freq)
        efluxJy.append('--')
        filter.append('SPLASH3')
        tab.append('J/ApJS/247/5/table1')

    url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=J/ApJS/239/15/table2&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    p = open(xmlfile, "wb")
    p.write(req.content)
    votable = parse(xmlfile)
    table = votable.get_first_table()
    data = table.array
    p.close()

    if len(data)>0:
        sp2flag = 1
        x = str(data[0]).split(", ")

        hmsra = x[4][2:-1]
        hmsdec = x[5][2:-1]
        c = SkyCoord(str(hmsra)+' '+str(hmsdec), unit=(u.hourangle, u.deg))
        righ.append(float(c.ra.degree))
        dela.append(float(c.dec.degree))
        fluxJy.append(float(x[6]))
        freq = float(x[3][-6:-2])
        fGhz.append(float(freq))
        wmicron.append(299792458/freq)
        efluxJy.append('--')
        filter.append('SPLASH2')
        tab.append('J/ApJS/239/15/table2')

    url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=J/ApJS/227/26/table1&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    p = open(xmlfile, "wb")
    p.write(req.content)
    votable = parse(xmlfile)
    table = votable.get_first_table()
    data = table.array
    p.close()

    if len(data)>0:
        sp1flag = 1
        x = str(data[0]).split(", ")

        print(url)
        hmsra = x[3][2:-1]
        hmsdec = x[4][2:-1]
        print(hmsra, hmsdec)
        if hmsra != "" and hmsdec != "":
            c = SkyCoord(str(hmsra)+' '+str(hmsdec), unit=(u.hourangle, u.deg))
            righ.append(float(c.ra.degree))
            dela.append(float(c.dec.degree))
        else:
            righ.append(0)
            dela.append(0)
        fluxJy.append(float(x[6]))
        freq = float(x[3][-6:-2])
        fGhz.append(float(freq))
        wmicron.append(299792458/freq)
        efluxJy.append('--')
        filter.append('SPLASH1')
        tab.append('J/ApJS/227/26/table1')

    #IRAS ii/126/sources
    url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=II/126/sources&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    p = open(xmlfile, "wb")
    p.write(req.content)
    exfg = 0
    with open(xmlfile,'r') as x:
        lines = x.readlines()
        for line in lines:
            if 'fatalError' in line:
                exfg =1
    x.close()

    if exfg ==0:
        votable = parse(xmlfile, invalid = 'mask')
        table = votable.get_first_table()
        data = table.array
        p.close()
    else:
        data = []
        p.close()

    if len(data)>0:
        for j in range(len(data)):
            x = str(data[j]).split(", ")

            hmsra = x[3][2:-1]
            hmsdec = x[4][2:-1]
            c = SkyCoord(str(hmsra)+' '+str(hmsdec), unit=(u.hourangle, u.deg))
            righ.append(float(c.ra.degree))
            dela.append(float(c.dec.degree))
            fluxJy.append(float(x[6]))
            efluxJy.append('--')
            if str(x[10]) == '1':
                filter.append('IRAS:12-uplim')
            else:
                filter.append('IRAS:12')
            wmicron.append('12')
            tab.append('II/126/sources')


            righ.append(float(c.ra.degree))
            dela.append(float(c.dec.degree))
            fluxJy.append(float(x[7]))
            efluxJy.append('--')
            if str(x[11]) == '1':
                filter.append('IRAS:25-uplim')
            else:
                filter.append('IRAS:25')
            wmicron.append('25')
            tab.append('II/126/sources')


            righ.append(float(c.ra.degree))
            dela.append(float(c.dec.degree))
            fluxJy.append(float(x[8]))
            efluxJy.append('--')
            if str(x[12]) == '1':
                filter.append('IRAS:60-uplim')
            else:
                filter.append('IRAS:60')
            wmicron.append('60')
            tab.append('II/126/sources')


            righ.append(float(c.ra.degree))
            dela.append(float(c.dec.degree))
            fluxJy.append(float(x[9]))
            efluxJy.append('--')
            if str(x[13]) == '1':
                filter.append('IRAS:100-uplim')
            else:
                filter.append('IRAS:100')
            wmicron.append('100')
            tab.append('II/126/sources')

    #unWISE catalogue
    url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=II/363/unwise&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    p = open(xmlfile, "wb")
    p.write(req.content)
    exfg = 0
    with open(xmlfile,'r') as x:
        lines = x.readlines()
        for line in lines:
            if 'fatalError' in line:
                exfg =1
    x.close()

    if exfg ==0:
        votable = parse(xmlfile, invalid = 'mask')
        table = votable.get_first_table()
        data = table.array
        p.close()
    else:
        data = []
        p.close()

    if len(data)>0:
        for j in range(len(data)):
            x = str(data[j]).split(", ")

            righ.append(float(x[2]))
            dela.append(float(x[3]))
            mag = float(x[8])/309.05
            fluxJy.append(float(mag/1e9))
            emag = float(x[14])/309.05
            efluxJy.append(float(emag/1e9))
            wmicron.append(3.4)
            fGhz.append(88174.252352941)
            filter.append('unWISE:W1')
            tab.append('II/363/unwise')

            righ.append(float(x[2]))
            dela.append(float(x[3]))
            mag = float(x[9])/309.05
            fluxJy.append(float(mag/1e9))
            emag = float(x[15])/309.05
            efluxJy.append(float(emag/1e9))
            wmicron.append(4.6)
            fGhz.append(65172.273478261)
            filter.append('unWISE:W2')
            tab.append('II/363/unwise')


    #J/AcA/62/247/stars catalogue
    url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=J/AcA/62/247/stars&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    p = open(xmlfile, "wb")
    p.write(req.content)
    votable = parse(xmlfile)
    table = votable.get_first_table()
    data = table.array
    p.close()

    if len(data)>0:
        for j in range(len(data)):
            x = str(data[j]).split(", ")

            c = SkyCoord(str(x[3][2:-1])+' '+str(x[4][2:-1]), unit=(u.hourangle, u.deg))
            righ.append(float(c.ra.degree))
            dela.append(float(c.dec.degree))
            mag = float(x[7])
            fluxJy.append(float(2416.0*10**(-0.4*mag)))
            emag = float(x[11])
            efluxJy.append(float(2416.0*10**(-0.4*emag)))
            wmicron.append(0.8)
            fGhz.append(374740.5725)
            filter.append('OGLE:I')
            tab.append('J/AcA/62/247/stars')

            righ.append(float(c.ra.degree))
            dela.append(float(c.dec.degree))
            mag = float(x[8])
            fluxJy.append(float(3636.0*10**(-0.4*mag)))
            emag = float(x[12])
            efluxJy.append(float(3636.0*10**(-0.4*emag)))
            wmicron.append(0.55)
            fGhz.append(545077.19636)
            filter.append('OGLE:V')
            tab.append('J/AcA/62/247/stars')

    #SDSS
    url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=V/147/sdss12&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    p = open(xmlfile, "wb")
    p.write(req.content)
    votable = parse(xmlfile)
    table = votable.get_first_table()
    data = table.array
    p.close()

    if len(data)>0:
        for j in range(len(data)):
            x = str(data[j]).split(", ")

            righ.append(float(x[1]))
            dela.append(float(x[2]))
            mag = float(x[20])
            fluxJy.append(float(3631.0*10**(-0.4*mag)))
            emag = float(x[21])
            efluxJy.append(float(3631.0*10**(-0.4*emag)))
            wmicron.append(0.3543)
            fGhz.append(846154.27039232)
            filter.append('SDSS:u')
            tab.append('V/147/sdss12')

            righ.append(float(x[1]))
            dela.append(float(x[2]))
            mag = float(x[22])
            fluxJy.append(float(3631.0*10**(-0.4*mag)))
            emag = float(x[23])
            efluxJy.append(float(3631.0*10**(-0.4*emag)))
            wmicron.append(0.4770)
            fGhz.append(628495.71907757)
            filter.append('SDSS:g')
            tab.append('V/147/sdss12')

            righ.append(float(x[1]))
            dela.append(float(x[2]))
            mag = float(x[24])
            fluxJy.append(float(3631.0*10**(-0.4*mag)))
            emag = float(x[25])
            efluxJy.append(float(3631.0*10**(-0.4*emag)))
            wmicron.append(0.6231)
            fGhz.append(481130.56973199)
            filter.append('SDSS:r')
            tab.append('V/147/sdss12')

            righ.append(float(x[1]))
            dela.append(float(x[2]))
            mag = float(x[26])
            fluxJy.append(float(3631.0*10**(-0.4*mag)))
            emag = float(x[27])
            efluxJy.append(float(3631.0*10**(-0.4*emag)))
            wmicron.append(0.7625)
            fGhz.append(393170.43672131)
            filter.append('SDSS:i')
            tab.append('V/147/sdss12')

            righ.append(float(x[1]))
            dela.append(float(x[2]))
            mag = float(x[28])
            fluxJy.append(float(3631.0*10**(-0.4*mag)))
            emag = float(x[29])
            efluxJy.append(float(3631.0*10**(-0.4*emag)))
            wmicron.append(0.9134)
            fGhz.append(328215.96014889)
            filter.append('SDSS:z')
            tab.append('V/147/sdss12')
    rejectlist = obj[i] + ":\n"

    offset = []
    #calculate diff in coords in arcsec
    for m in range(len(righ)):
        c1 = SkyCoord(ra[i]*u.deg, dec[i]*u.deg)
        c2 = SkyCoord(righ[m]*u.deg, dela[m]*u.deg)
        sep = c1.separation(c2)
        offset.append(float(sep.arcsecond))


    #fixing labels
    for n in range(len(filter)):
        if "Spitzer" in filter[n]:
            filter[n] = "ii/293/glimpse/" + filter[n][10:]
        if "24um" in filter[n]:
            filter[n] = "mipsgal:24um"
        if "allwise" in tab[n] and "WISE" in filter[n]:
            filter[n] = 'ALLWISE' + ":" + filter[n][-1:]
        if "J/ApJS/188/123/table1" in tab[n] or "J/ApJ/799/29/table5" in tab[n] or "J/ApJS/208/14/bgpsv2_1" in tab[n] or "J/ApJS/195/14/table3" in tab[n] or 'J/ApJS/218/1/table3' in tab[n] or 'J/ApJ/805/157/table1' in tab[n] or "J/ApJ/741/110/table1" in tab[n]:
            filter[n] = 'BOLOCAM'
        if "J/A+A/599/A139/table1" in tab[n] or "J/A+A/549/A45/atl-csc" in tab[n] or "J/A+A/568/A41/atl-csc" in tab[n] or "J/MNRAS/443/1555/table3" in tab[n] or 'J/MNRAS/435/400/table3' in tab[n]:
            filter[n] = 'ATLASGAL'
        if "J/ApJS/221/26/table2" in tab[n]:
            filter[n] = "WISE"
        if "VIII/88/erc" in tab[n] or 'J/A+A/619/A94/pcnt' in tab[n] or "VIII/91/pccs1" in tab[n] or 'J/A+A/594/A26/pccs217e' in tab[n]:
            filter[n] = 'PLANCK'
        if "VIII/85A/waste" in tab[n] or 'VIII/85A/spectra' in tab[n]:
            filter[n] = 'SPECFIND'
        if "J/A+A/501/539/table4" in tab[n] or 'J/A+A/461/11/radio' in tab[n]:
            filter[n] = 'RMS-Survey'
        if "J/MNRAS/471/100/hcatalog" in tab[n]:
            filter[n] = 'HI-GAL'
        if "J/MNRAS/469/2163/table2" in tab[n] or 'J/MNRAS/453/4264/table1' in tab[n]:
            filter[n] = 'JPS'
        if "J/A+A/627/A175/table7" in tab[n]:
            filter[n] = 'GLOSTAR'
        if "J/A+A/588/A97/catalog" in tab[n] or "J/A+A/619/A124" in tab[n]:
            filter[n] = 'THOR'
        if 'II/161/catalog' in tab[n]:
            filter[n] = 'EIC'
        if 'J/ApJ/737/45/table1' in tab[n]:
            filter[n] = 'FIRST-NVSS(VAR)'
        if 'J/ApJS/91/347/table2' in tab[n]:
            filter[n] = 'VLA-INNERGAL(BECKER94)'
        if 'J/AJ/92/787/table1' in tab[n]:
            filter[n] = 'VLA-OUTERGAL(FICH86)'
        if 'J/A+A/579/A71/tablea1' in tab[n]:
            filter[n] = 'HI-GAL-CORNISH(HIIREGIONS)'
        if 'J/A+A/537/A1/IRfluxes' in tab[n]:
            filter[n] = 'PN-HIIREGIONS(ANDERSON12)'
        if 'J/ApJ/838/139/outside' in tab[n]:
            filter[n] = '3FGL-ATCA-VLA(SCHINZEL17)'
        if 'J/ApJ/727/114/table1' in tab[n]:
            filter[n] = 'BLAST'
        if 'J/ApJS/205/1/catalog' in tab[n] or 'J/A+A/615/A103/uchiicat' in tab[n]:
            filter[n] = 'CORNISH'
        if ('II/305/archive' in tab[n] or 'II/305/catalog' in tab[n]) and ":=" in filter[n]:
            filter[n] = 'SAGE'
        if "II/181/sources" in tab[n] or "VII/185/catalog" in tab[n] or "II/274/iras_r" in tab[n] or "VII/221/pscz" in tab[n]:
            filter[n] = 'IRAS'
        if "J/A+A/531/A26/table1" in tab[n]:
            filter[n] = 'W31(BEUTHER13)'
        if 'J/MNRAS/363/405/table5' in tab[n]:
            filter[n] = "SEST"
        if 'VIII/81B/sumss212' in tab[n]:
            filter[n] = 'SUMSS'
        if "II/336/apass9" in tab[n]:
            filter[n] = "APASS(" + filter[n] + ")"
        if "J/MNRAS/402/2403/at20gcat" in tab[n] or 'J/MNRAS/434/956/table2' in tab[n] or 'J/MNRAS/417/2651/catalog' in tab[n]:
            filter[n] = 'AT20G'
        if 'J/ApJS/171/61/table5' in tab[n]:
            filter[n] = 'CRATES'
        if 'IX/10A/cor_iras' in tab[n]:
            filter[n] = '1RXS'
        if 'J/ApJS/217/17/table2' in tab[n]:
            filter[n] = 'SAFIRES'
        if "J/other/RAA/10.67/table1" in tab[n] or 'J/A+A/550/A21/obs' in tab[n]:
            filter[n] = 'Lin:' + filter[n]
            lins.append([])
            lins[-1].append(righ[n])
            lins[-1].append(dela[n])
            lins[-1].append(offset[n])
            lins[-1].append(filter[n])
            lins[-1].append(wmicron[n])
            lins[-1].append(fluxJy[n])
            lins[-1].append(efluxJy[n])
            lins[-1].append(tab[n])
        if "J/A+A/547/A49/tablec" in tab[n]:
            filter[n] = 'Herschel'
        if "J/AJ/131/2525/table2" in tab[n]:
            filter[n] = 'MAGPIS'
        if 'J/ApJ/749/47/table3' in tab[n] or 'J/A+A/532/A127/table1' in tab[n] or 'VIII/42/txs' in tab[n] or "J/ApJS/194/32/hrds" in tab[n]:
            filter[n] = 'Lin:' + filter[n]
        if "J/MNRAS/451/3089/photo" in tab[n]:
            filter[n] = 'IRDC'
        if 'J/AJ/130/156/match' in tab[n]:
            filter[n] = 'HII-REGIONS(GIVEON05)'
        if 'II/53/catalog' in tab[n]:
            filter[n] = 'CELOBJ(HALL74)'




    print('Querying IRSA')
#IRSA
    pacs = ['ppsc_70', 'ppsc_100', 'ppsc_160']
    pacs_w = [70, 100, 160]
    spire = ['spsc250', 'spsc350', 'spsc500']
    spire_w = [250, 350, 500]

    #PACS
    snrnoise = []
    for g in range(len(pacs)):

        #get pacs data from irsa
        url = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog='+ pacs[g] + '&spatial=cone&radius=' + str(rad) + '&radunits='+'arcsec&outfmt=3&objstr='+ str(ra[i]) + '+' + str(dec[i])
        req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
        #req = requests.get(url,allow_redirects=True, verify=False, headers=headers)

        #put photometry in a votable and extract info
        if flag == '1':
            xmlfile = subdir + "/" + str(obj[i]) + "_pacs_phot.xml"
        else:
            xmlfile = dir + "/" + "phot.xml"
        p = open(xmlfile, "wb")
        p.write(req.content)
        votable = parse(xmlfile)
        table = votable.get_first_table()
        data = table.array
        data = str(data).split(", ")
        p.close()
        if len(data) < 2:
            continue
        else:
            #extract relevant information from data list and append
            righ.append(float(data[2]))
            dela.append(float(data[3]))
            offset.append(float(data[-2]))
            indx = pacs[g].index('_')
            filter.append(str(pacs[g][:indx])+str(pacs[g][indx+1:]))
            wmicron.append(float(pacs_w[g]))
            fluxJy.append(float(data[8])/(10**3))
            efluxJy.append(float(data[11])/(10**3))
            tab.append("-")
            #snrnoise.append(data[10])

    #SPIRE
    for g in range(len(spire)):

        #get spire data from irsa
        url = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog='+ spire[g] + '&spatial=cone&radius=' + str(rad) + '&radunits='+'arcsec&outfmt=3&objstr='+ str(ra[i]) + '+' + str(dec[i])
        req = requests.get(url,allow_redirects=True, verify=False, headers=headers)

        #put photometry in a votable and extract info
        if flag == '1':
            xmlfile = subdir + "/" + str(obj[i]) + "_spire_phot.xml"
        else:
            xmlfile = dir + "/" + "phot.xml"
        p = open(xmlfile, "wb")
        p.write(req.content)
        votable = parse(xmlfile)
        table = votable.get_first_table()
        data = table.array
        data = str(data).split(", ")
        p.close()
        if len(data) < 2:
            continue
        else:
            #extract relevant information from data list and append
            righ.append(float(data[2]))
            dela.append(float(data[3]))
            offset.append(float(data[-2]))
            filter.append(str(spire[g]))
            wmicron.append(float(spire_w[g]))
            fluxJy.append(float(data[13])/(10**3))
            efluxJy.append(float(data[14])/(10**3))
            tab.append("-")

    #CATWISE
    #begin query irsa cat
    #get catwise data from irsa
    url = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=cwcat2&spatial=cone&radius=' + str(rad) + '&radunits='+'arcsec&outfmt=3&objstr='+ str(ra[i]) + '+' + str(dec[i])
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    #put photometry in a votable and extract info
    print(url)
    p = open(xmlfile, "wb")
    if "The catalog is not on the list" not in str(req.content):

        p.write(req.content)
        votable = parse(xmlfile)
        table = votable.get_first_table()
        data = table.array
        p.close()
        #end query irsa cat
        if len(data) >0:
            #extract relevant information from data list and append
            data = str(data).split(", ")

            #w1
            if data[34] != '--':
                righ.append(float(data[2]))
                dela.append(float(data[3]))
                offset.append(float(data[-2]))
                filter.append('CATWISE:W1')
                wmicron.append(3.35001071)
                w1mag = float(data[34])
                flux = 306.68*10**(-0.4*w1mag)
                fluxJy.append(flux)
                if data[35] == '--':
                    eflux = data[35]
                else:
                    ew1mag = float(data[35])
                    eflux = 306.68*10**(-0.4*ew1mag)
                efluxJy.append(eflux)
                tab.append("-")

            #w2
            if data[38] != '--':
                righ.append(float(data[2]))
                dela.append(float(data[3]))
                offset.append(float(data[-2]))
                filter.append('CATWISE:W2')
                wmicron.append(4.6000193)
                w2mag = float(data[38])
                flux2 = 170.66*10**(-0.4*w2mag)
                fluxJy.append(flux2)
                if data[39] == '--':
                    eflux2 = data[39]
                else:
                    ew2mag = float(data[39])
                    eflux2 = 170.66*10**(-0.4*ew2mag)
                efluxJy.append(eflux2)
                tab.append("-")

    #throwing out points with flux/error < 3
    efl = []
    flu = []
    fl = []
    r = []
    d = []
    t = []
    w = []
    o = []
    for l in range(len(efluxJy)):
        if efluxJy[l] != '--' and float(efluxJy[l]) != 0:
            rat = float(float(fluxJy[l])/float(efluxJy[l]))
            if rat >= 3:
                efl.append(efluxJy[l])
                flu.append(fluxJy[l])
                fl.append(filter[l])
                r.append(righ[l])
                d.append(dela[l])
                t.append(tab[l])
                w.append(wmicron[l])
                o.append(offset[l])
            else:
                why = filter[l] + ' rejected because flux/error is ' + str(rat) + '\n'
                rejectlist += why
        elif efluxJy[l] =='--' or float(efluxJy[l]) == 0:
            efl.append(efluxJy[l])
            flu.append(fluxJy[l])
            fl.append(filter[l])
            r.append(righ[l])
            d.append(dela[l])
            t.append(tab[l])
            w.append(wmicron[l])
            o.append(offset[l])
    efluxJy = efl
    fluxJy = flu
    filter = fl
    righ = r
    dela = d
    tab = t
    wmicron = w
    offset = o

    #throwing out points with offset > astrometric accuracy
    efl = []
    flu = []
    fl = []
    r = []
    d = []
    t = []
    w = []
    o = []


    for l in range(len(filter)):
        fg2 = 0
        efl2 = efl
        if "MSX" in filter[l] and "V/114/msx6_gp" in tab[l]:
            url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=' + urllib.parse.quote('V/114/msx6_gp')+ '&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
            req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
            p = open(xmlfile, "wb")
            p.write(req.content)
            votable = parse(xmlfile)
            table = votable.get_first_table()
            data = table.array
            p.close()

            x = str(data[0]).split(", ")
            am = (float(x[5])**2 + float(x[6])**2)**(1/2)
            tot = (am**2+inp_cerr**2)**(1/2)
            if offset[l] != "":
                fg2 = 1
            if offset[l] != "" and offset[l]!= "." and tot > float(offset[l]):
                efl.append(efluxJy[l])
                flu.append(fluxJy[l])
                fl.append(filter[l])
                r.append(righ[l])
                d.append(dela[l])
                t.append(tab[l])
                w.append(wmicron[l])
                o.append(offset[l])
            elif fg2 ==1:
                rejectlist+= filter[l] + ' rejected because offset is ' + str(offset[l]) + ', astrometric uncertainty is ' + str(tot) + '\n'
        elif "IRAS" in filter[l] and "II/125/main" in tab[l]:
            url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=' + urllib.parse.quote('II/125/main')+ '&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
            req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
            p = open(xmlfile, "wb")
            p.write(req.content)
            votable = parse(xmlfile)
            table = votable.get_first_table()
            data = table.array
            p.close()

            x = str(data[0]).split(", ")
            if filter[l][-2:]=='12':
                if str(x[11]) == '1':
                    filter[l] = filter[l] + '-uplim'
            if filter[l][-2:]=='25':
                if str(x[14]) == '1':
                    filter[l] = filter[l] + '-uplim'
            if filter[l][-2:]=='60':
                if str(x[17]) == '1':
                    filter[l] = filter[l] + '-uplim'
            if filter[l][-2:]=='100':
                if str(x[20]) == '1':
                    filter[l] = filter[l] + '-uplim'
            am = (float(x[5])**2 + float(x[6])**2)**(1/2)
            tot = (am**2+inp_cerr**2)**(1/2)
            if offset[l] != "":
                fg2 = 1
            if offset[l] != "" and offset[l]!= "." and tot > float(offset[l]):
                efl.append(efluxJy[l])
                flu.append(fluxJy[l])
                fl.append(filter[l])
                r.append(righ[l])
                d.append(dela[l])
                t.append(tab[l])
                w.append(wmicron[l])
                o.append(offset[l])
            elif fg2 ==1:
                rejectlist+= filter[l] + ' rejected because offset is ' + str(offset[l]) + ', astrometric uncertainty is ' + str(tot) + '\n'
        elif "EIC" in filter[l] and "II/161/catalog" in tab[l]:
            url = 'https://vizier.u-strasbg.fr/viz-bin/votable?-source=' + urllib.parse.quote('II/161/catalog')+ '&-c='+ urllib.parse.quote(points)+'&-c.rs=' + str(rad) + '&-out.all'
            req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
            p = open(xmlfile, "wb")
            p.write(req.content)
            votable = parse(xmlfile)
            table = votable.get_first_table()
            data = table.array
            p.close()

            x = str(data[0]).split(", ")
            am = ((float(x[6])*15)**2 + float(x[8])**2)**(1/2)
            tot = (am**2+inp_cerr**2)**(1/2)
            if offset[l] != "":
                fg2 = 1
            if offset[l] != "" and offset[l]!= "." and tot > float(offset[l]):
                efl.append(efluxJy[l])
                flu.append(fluxJy[l])
                fl.append(filter[l])
                r.append(righ[l])
                d.append(dela[l])
                t.append(tab[l])
                w.append(wmicron[l])
                o.append(offset[l])
            elif fg2 ==1:
                rejectlist+= filter[l] + ' rejected because offset is ' + str(offset[l]) + ', astrometric uncertainty is ' + str(tot) + '\n'
        else:
            for m in range(len(amac)):
                if amacmission[m] in filter[l]:
                    tot = (amac[m]**2+inp_cerr**2)**(1/2)
                    if offset[l] != "":
                        fg2 = 1
                    if offset[l] != "" and offset[l]!= "." and tot > float(offset[l]):
                        efl.append(efluxJy[l])
                        flu.append(fluxJy[l])
                        fl.append(filter[l])
                        r.append(righ[l])
                        d.append(dela[l])
                        t.append(tab[l])
                        w.append(wmicron[l])
                        o.append(offset[l])
                    elif fg2 ==1:
                        rejectlist+= filter[l] + ' rejected because offset is ' + str(offset[l]) + ' , astrometric uncertainty is ' + str(tot) + '\n'

        if len(efl2) == len(efl) and fg2 == 0:
            efl.append(efluxJy[l])
            flu.append(fluxJy[l])
            fl.append(filter[l])
            r.append(righ[l])
            d.append(dela[l])
            t.append(tab[l])
            w.append(wmicron[l])
            o.append(offset[l])

    efluxJy = efl
    fluxJy = flu
    filter = fl
    righ = r
    dela = d
    tab = t
    wmicron = w
    offset = o

    #getting rid of certain catalogues
    efl = []
    flu = []
    fl = []
    r = []
    d = []
    t = []
    w = []
    o = []

    for l in range(len(tab)):
        if "J/AJ/122/1844/table2" not in tab[l]:
            efl.append(efluxJy[l])
            flu.append(fluxJy[l])
            fl.append(filter[l])
            r.append(righ[l])
            d.append(dela[l])
            t.append(tab[l])
            w.append(wmicron[l])
            o.append(offset[l])

    efluxJy = efl
    fluxJy = flu
    filter = fl
    righ = r
    dela = d
    tab = t
    wmicron = w
    offset = o


    #fixing names for repeated filters
    filterch =filter
    for j in range(len(filter)):
        fil = filter[j]
        count = 2
        for k in range(len(filter)):
            if fil == filter[k] and j!=k and ":=" not in fil:
                filter[j] = fil + "(1)"
                filter[k] = filter[k] + "(" + str(count) + ")"
                count = count+1
    efl = []
    flu = []
    fl = []
    r = []
    d = []
    t = []
    w = []
    o = []

    for l in range(len(filter)):
        if 'Lin:' not in filter[l]:
            efl.append(efluxJy[l])
            flu.append(fluxJy[l])
            fl.append(filter[l])
            r.append(righ[l])
            d.append(dela[l])
            t.append(tab[l])
            w.append(wmicron[l])
            o.append(offset[l])
    for j in range(len(lins)):
        r.append(lins[j][0])
        d.append(lins[j][1])
        o.append(lins[j][2])
        fl.append(lins[j][3])
        w.append(lins[j][4])
        flu.append(lins[j][5])
        efl.append(lins[j][6])
        t.append(lins[j][7])

    efluxJy = efl
    fluxJy = flu
    filter = fl
    righ = r
    dela = d
    tab = t
    wmicron = w
    offset = o
    #determine how long header is
    if len(filter) > num:
            num = len(filter)

    #update mission list
    for e in range(len(filter)):
        if filter[e] not in mission:
            mission.append(filter[e])
            tab2.append(tab[e])

    #find type of object
    url = 'http://simbad.u-strasbg.fr/simbad/sim-coo?Coord='+urllib.parse.quote(points)+'&Radius=5&Radius.unit=arcsec&output.format=ASCII&list.otypesel=on&otypedisp=S'
    req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
    nam = dir + "/" + "simbad-objtype.out"
    p = open(nam, "wb")
    p.write(req.content)
    p.close()

    lines = 'str'
    with open(nam, "r") as f:
        lines = f.readlines()
        if len(lines)<6:
            type = '--'
            lines = ""
        else:
            lines = lines[5]
    if '---' in lines:
        ind = lines.index('---')
        ind2 = lines.index('---', ind+1)

        type = lines[ind+3:ind2].strip()
    else:
        type = '--'

    #create list that has all the info together
    row = []

    row.append(str(obj[i]))
    row.append(str(ra[i]))
    row.append(str(dec[i]))
    row.append(str(type))

    for k in range(len(wmicron)):
        row.append(str(filter[k]))
        row.append(str(tab[k]))
        row.append(str(righ[k]))
        row.append(str(dela[k]))
        row.append(str(offset[k]))
        row.append(str(wmicron[k]))
        row.append(str(fluxJy[k]))
        row.append(str(efluxJy[k]))

#create .csv file that has photometry of each object in one row
    with open(noheader, "a") as photometry:
        pm_writer = csv.writer(photometry, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        pm_writer.writerow(row)
    with open(rejects, "a") as r:
        r.write(rejectlist)
photometry.close()

print('creating header')
#creating header
header = ['object', 'ra(source)', 'dec(source)', 'type']

for w in range(num):
    header.append('filter')
    header.append('catalogue')
    header.append('ra')
    header.append('dec')
    header.append('offset(arcsec)')
    header.append('wavelength(micron)')
    header.append('flux(Jy)')
    header.append('error_flux(Jy)')

#writing header
h = dir + "/" + root + "_phot_raw.csv"
with open(h, "w") as ph, open(noheader, "r") as photometry:
    p_writer = csv.writer(ph, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    p_writer.writerow(header)
    csv_reader = csv.reader(photometry, delimiter=',')
    for row in csv_reader:
        p_writer.writerow(row)
photometry.close()
ph.close()

os.remove(noheader)

#writing mission list
missionfile = dir + "/" + root +  "_missions.csv"
mission2d = []
for i in range(len(mission)):
    new = []
    new.append(mission[i])
    new.append(' '+ tab2[i])
    mission2d.append(new)


with open(missionfile, "w") as missions:
    m_writer = csv.writer(missions, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for a in range(len(mission2d)):
        m_writer.writerows([mission2d[a]])
missions.close()

with open(dir+ '/' + coordfile, 'w') as c, open(coordfile, 'r') as corig:
    c_writer = csv.writer(c, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    csv_reader = csv.reader(corig, delimiter=',')
    for row in csv_reader:
        c_writer.writerow(row)
c.close()
corig.close()

copy(sys.argv[1], dir)


if flag != '1':
    os.remove(dir + "/" + "phot.xml")



#DIRECTORY PART
print('creating directory')
data = dir + "/" + root + "_phot_raw.csv"
mis = dir + "/" + root + "_missions.csv"

missions = []
#cata = []
cats = []
cata = []

with open(mis, mode = 'r') as m:
    m_reader = csv.reader(m, delimiter=',')
    for row in m_reader:
        #row = str(row)
        missions.append((row)[0])
        cata.append((row)[1])
        #if len(row) > 1:
        #    cata.append(row[1])
        #else:
        #    cata.append('0')
m.close()
#info = []
mis_nml = []
mis_ml = []
#molline #1
for i in range(len(missions)):
    if "J/ApJ/702/1615/masers" in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/ApJ/702/1615/masers')
    elif ":=1612" in missions[i] and 'J/A+AS/90/327/catalog' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/A+AS/90/327/catalog')
    elif ":=1.1" in missions[i] and 'J/ApJ/799/29/table5' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/ApJ/799/29/table5')
    elif ":=6668" in missions[i] and 'VIII/96/catalog' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('VIII/96/catalog')
    elif ":=6.7" in missions[i] and 'J/AJ/152/92/table1' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/AJ/152/92/table1')
    elif ":=36.169" in missions[i] and 'J/ApJS/227/10/table2' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/ApJS/227/10/table2')
    elif 'J/AN/333/634/tablea1' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/AN/333/634/tablea1')
    elif 'J/A+A/569/A125/tablea1' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/A+A/569/A125/tablea1')
    elif 'J/MNRAS/439/2584/table1' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/MNRAS/439/2584/table1')
    elif 'J/A+A/325/255/table1' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/A+A/325/255/table1')
    elif 'J/A+A/368/845/table2' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/A+A/368/845/table2')
    elif 'J/A+A/399/1083/table1' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/A+A/399/1083/table1')
    elif 'J/AZh/79/328/table2' in cata[i]:
        mis_ml.append(missions[i])
        cats.append('J/AZh/79/328/table2')
    elif missions[i] == 'SPLASH3':
        mis_ml.append(missions[i])
        cats.append('J/ApJS/247/5/table1')
    elif missions[i] == 'SPLASH2':
        mis_ml.append(missions[i])
        cats.append('J/ApJS/239/15/table2')
    elif missions[i] == 'SPLASH1':
        mis_ml.append(missions[i])
        cats.append('J/ApJS/227/26/table1')
    else:
        mis_nml.append(missions[i])

header = ['object', 'ra(source)', 'dec(source)', 'type']
for i in range(len(mis_nml)):
    header.append(mis_nml[i]+'_catalogue')
    header.append(mis_nml[i]+'_ra')
    header.append(mis_nml[i]+'_dec')
    header.append(mis_nml[i]+'_offset(arcsec)')
    header.append(mis_nml[i]+'_wavelength(micron)')
    header.append(mis_nml[i]+'_flux(Jy)')
    header.append(mis_nml[i]+'_error_flux(Jy)')

#molline #2
for i in range(len(mis_ml)):

    if ":=1612" in mis_ml[i] and 'J/A+AS/90/327/catalog' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_peaks')
        header.append('Lin'+mis_ml[i]+'_VL(km/s)')
        header.append('Lin'+mis_ml[i]+'_VH(km/s)')
    elif 'J/A+A/569/A125/tablea1' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VLSR(km/s)')
    elif 'J/A+A/399/1083/table1' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VLSR(km/s)')
    elif 'J/AZh/79/328/table2' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_1VLSR(km/s)')
        header.append('Lin'+mis_ml[i]+'_2VLSR(km/s)')
    elif 'J/MNRAS/439/2584/table1' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_1VL(km/s)')
        header.append('Lin'+mis_ml[i]+'_1VH(km/s)')
        header.append('Lin'+mis_ml[i]+'_2VL(km/s)')
        header.append('Lin'+mis_ml[i]+'_2VH(km/s)')
    elif 'J/A+A/325/255/table1' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VL(km/s)')
        header.append('Lin'+mis_ml[i]+'_VH(km/s)')
    elif 'J/ApJ/702/1615/masers' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VLSR(km/s)')
    elif 'J/AN/333/634/tablea1' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VLSR(km/s)')
    elif ":=1.1" in mis_ml[i] and 'J/ApJ/799/29/table5' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VLSR(km/s)')
    elif ":=6668" in mis_ml[i] and 'VIII/96/catalog' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VL(km/s)')
        header.append('Lin'+mis_ml[i]+'_VH(km/s)')
    elif 'J/A+A/368/845/table2' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VL(km/s)')
        header.append('Lin'+mis_ml[i]+'_VH(km/s)')
        header.append('Lin'+mis_ml[i]+'_VPEAK(km/s)')
    elif ':=6.7' in mis_ml[i] and 'J/AJ/152/92/table1' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_VPEAK(km/s)')
        header.append('Lin'+mis_ml[i]+'_dist(kpc)')
        header.append('Lin'+mis_ml[i]+'_lum(Lsun)')
    elif ':=36.169' in mis_ml[i] and 'J/ApJS/227/10/table2' in cats[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_vel(km/s)')
    elif 'SPLASH' in mis_ml[i]:
        header.append('Lin'+mis_ml[i]+'_catalogue')
        header.append('Lin'+mis_ml[i]+'_ra')
        header.append('Lin'+mis_ml[i]+'_dec')
        header.append('Lin'+mis_ml[i]+'_offset(arcsec)')
        header.append('Lin'+mis_ml[i]+'_wavelength(micron)')
        header.append('Lin'+mis_ml[i]+'_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_error_flux(Jy)')
        header.append('Lin'+mis_ml[i]+'_iflux(Jy*km/s)')
        header.append('Lin'+mis_ml[i]+'_pkvel(km/s)')
        header.append('Lin'+mis_ml[i]+'_minvel(km/s)')
        header.append('Lin'+mis_ml[i]+'_maxvel(km/s)')
        header.append('Lin'+mis_ml[i]+'_assoc')


#MOLECULAR LINE header

#h20
molin = path + '/' + str(root) + "_h2o_data.csv"
count = 0
if os.path.exists(molin):
    with open(molin, 'r') as m:
        lines = m.readlines()[2:]
        for line in lines:
            x = line.split(",")
            if "Yes" in x[4]:
                count = count+1

    m.close()
    if count == 1:
        header.append('Lin:H2O_ra')
        header.append('Lin:H2O_dec')
        header.append('Lin:H2O_angdist(arcsec)')
        header.append('Lin:H2O_fluxpeak(Jy)')
        header.append('Lin:H2O_fluxpeak_error(Jy)')
        header.append('Lin:H2O_vpeak(km/s)')
        header.append('Lin:H2O_vpeak_error(km/s)')
    else:
        for k in range(count):
            header.append('Lin:H2O(' + str(k+1) + ')'+ '_ra')
            header.append('Lin:H2O(' + str(k+1) + ')'+ '_dec')
            header.append('Lin:H2O(' + str(k+1) + ')'+'_angdist(arcsec)')
            header.append('Lin:H2O(' + str(k+1) + ')'+ '_fluxpeak(Jy)')
            header.append('Lin:H2O(' + str(k+1) + ')'+ '_fluxpeak_error(Jy)')
            header.append('Lin:H2O(' + str(k+1) + ')'+ '_vpeak(km/s)')
            header.append('Lin:H2O(' + str(k+1) + ')'+'_vpeak_error(km/s)')

#oh
molin = path + '/' + str(root) + "_oh_data.csv"
count = 0
if os.path.exists(molin):
    with open(molin, 'r') as m:
        lines = m.readlines()[2:]
        for line in lines:
            x = line.split(",")
            if "Yes" in x[4]:
                count = count+1

    m.close()
    if count == 1:
        header.append('Lin:OH_ra')
        header.append('Lin:OH_dec')
        header.append('Lin:OH_angdist(arcsec)')
        header.append('Lin:OH_freq(MHz)')
        header.append('Lin:OH_fluxpeak(Jy)')
        header.append('Lin:OH_fluxpeak_error(Jy)')
        header.append('Lin:OH_vpeak(km/s)')
        header.append('Lin:OH_vpeak_error(km/s)')
        header.append('Lin:OH_vpeak2(km/s)')
        header.append('Lin:OH_vpeak2_error(km/s)')
    else:
        for k in range(count):
            header.append('Lin:OH(' + str(k+1) + ')'+'_ra')
            header.append('Lin:OH(' + str(k+1) + ')'+'_dec')
            header.append('Lin:OH(' + str(k+1) + ')'+'_angdist(arcsec)')
            header.append('Lin:OH(' + str(k+1) + ')'+'_freq(MHz)')
            header.append('Lin:OH(' + str(k+1) + ')'+'_fluxpeak(Jy)')
            header.append('Lin:OH(' + str(k+1) + ')'+'_fluxpeak_error(Jy)')
            header.append('Lin:OH(' + str(k+1) + ')'+'_vpeak(km/s)')
            header.append('Lin:OH(' + str(k+1) + ')'+'_vpeak_error(km/s)')
            header.append('Lin:OH(' + str(k+1) + ')'+'_vpeak2(km/s)')
            header.append('Lin:OH(' + str(k+1) + ')'+'_vpeak2_error(km/s)')

#sio
molin = path + '/' + str(root) + "_sio_data.csv"
count = 0
if os.path.exists(molin):
    with open(molin, 'r') as m:
        lines = m.readlines()[2:]
        for line in lines:
            x = line.split(",")
            if "Yes" in x[4]:
                count = count+1

    m.close()
    if count ==1:
        header.append('Lin:SIO_ra')
        header.append('Lin:SIO_dec')
        header.append('Lin:SIO_angdist(arcsec)')
        header.append('Lin:SIO_fluxpeak(Jy)')
        header.append('Lin:SIO_fluxpeak_error(Jy)')
        header.append('Lin:SIO_vpeak(km/s)')
        header.append('Lin:SIO_vpeak_error(km/s)')
    else:
        for k in range(count):
            header.append('Lin:SIO(' + str(k+1) + ')'+'_ra')
            header.append('Lin:SIO(' + str(k+1) + ')'+'_dec')
            header.append('Lin:SIO(' + str(k+1) + ')'+'_angdist(arcsec)')
            header.append('Lin:SIO(' + str(k+1) + ')'+'_fluxpeak(Jy)')
            header.append('Lin:SIO(' + str(k+1) + ')'+'_fluxpeak_error(Jy)')
            header.append('Lin:SIO(' + str(k+1) + ')'+'_vpeak(km/s)')
            header.append('Lin:SIO(' + str(k+1) + ')'+'_vpeak_error(km/s)')

#ch30h class i
molin = path + '/' + str(root) + "_met1_data.csv"
count = 0
if os.path.exists(molin):
    with open(molin, 'r') as m:
        lines = m.readlines()[2:]
        for line in lines:
            x = line.split(",")
            if "Yes" in x[4]:
                count = count+1

    m.close()
    if count ==1:
        header.append('Lin:CH3OHclass1_ra')
        header.append('Lin:CH3OHclass1_dec')
        header.append('Lin:CH3OHclass1_angdist(arcsec)')
        header.append('Lin:CH3OHclass1_fluxpeak(Jy)')
        header.append('Lin:CH3OHclass1_fluxpeak_error(Jy)')
        header.append('Lin:CH3OHclass1_vpeak(km/s)')
        header.append('Lin:CH3OHclass1_vpeak_error(km/s)')
    else:
        for k in range(count):
            header.append('Lin:CH3OHclass1(' + str(k+1) + ')'+'_ra')
            header.append('Lin:CH3OHclass1(' + str(k+1) + ')'+'_dec')
            header.append('Lin:CH3OHclass1(' + str(k+1) + ')'+'_angdist(arcsec)')
            header.append('Lin:CH3OHclass1(' + str(k+1) + ')'+'_fluxpeak(Jy)')
            header.append('Lin:CH3OHclass1(' + str(k+1) + ')'+'_fluxpeak_error(Jy)')
            header.append('Lin:CH3OHclass1(' + str(k+1) + ')'+'_vpeak(km/s)')
            header.append('Lin:CH3OHclass1(' + str(k+1) + ')'+'_vpeak_error(km/s)')

#ch30h class ii
molin = path + '/' + str(root) + "_met2_data.csv"
count = 0
if os.path.exists(molin):
    with open(molin, 'r') as m:
        lines = m.readlines()[2:]
        for line in lines:
            x = line.split(",")
            if "Yes" in x[4]:
                count = count+1

    m.close()
    if count ==1:
        header.append('Lin:CH3OHclass2_ra')
        header.append('Lin:CH3OHclass2_dec')
        header.append('Lin:CH3OHclass2_angdist(arcsec)')
        header.append('Lin:CH3OHclass2_fluxpeak(Jy)')
        header.append('Lin:CH3OHclass2_fluxpeak_error(Jy)')
        header.append('Lin:CH3OHclass2_vpeak(km/s)')
        header.append('Lin:CH3OHclass2_vpeak_error(km/s)')
    else:
        for k in range(count):
            header.append('Lin:CH3OHclass2(' + str(k+1) + ')'+'_ra')
            header.append('Lin:CH3OHclass2(' + str(k+1) + ')'+'_dec')
            header.append('Lin:CH3OHclass2(' + str(k+1) + ')'+'_angdist(arcsec)')
            header.append('Lin:CH3OHclass2(' + str(k+1) + ')'+'_fluxpeak(Jy)')
            header.append('Lin:CH3OHclass2(' + str(k+1) + ')'+'_fluxpeak_error(Jy)')
            header.append('Lin:CH3OHclass2(' + str(k+1) + ')'+'_vpeak(km/s)')
            header.append('Lin:CH3OHclass2(' + str(k+1) + ')'+'_vpeak_error(km/s)')

#iras lrs spectra header
header.append('Sp:IRAS')
#SDSS optical spectra header
op = path + '/' + str(root) + "_optical.csv"
if os.path.exists(op):
    header.append('Sp:SDSS_optical')

dictphot = dir + "/" + root + "_phot.csv"
with open(dictphot, mode='w') as dict:
    dict_writer = csv.writer(dict, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    dict_writer.writerow(header)

rthd = dir + "/" + root + "_hdr_column.csv"
with open(rthd, 'w') as h:
    h_writer = csv.writer(h, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for i in range(len(header)):
        h_writer.writerow([header[i], " " + str(i+1)])

info = []
with open(data, mode='r') as phot:
    phot_reader = csv.reader(phot, delimiter=',')
    next(phot_reader)
    for row in phot_reader:
        info.append(row)

print('filling in molecular line data')
for i in range(len(info)):
    row = []
    row.append(info[i][0]) #obj
    row.append(info[i][1]) #ra
    row.append(info[i][2]) #dec
    row.append(info[i][3]) #type

    for k in range(len(mis_nml)):
        if mis_nml[k] not in info[i]:
            row.append(-999)
            row.append(-999)
            row.append(-999)
            row.append(-999)
            row.append(-999)
            row.append(-999)
            row.append(-999)
        else:
            index = info[i].index(mis_nml[k])
            row.append(info[i][index+1]) #catalogue
            row.append(info[i][index+2]) #ra
            row.append(info[i][index+3]) #dec
            row.append(info[i][index+4]) #offset
            row.append(info[i][index+5]) #wav
            row.append(info[i][index+6]) #f
            row.append(info[i][index+7]) #ef


    for k in range(len(mis_ml)):


        if "SPLASH" in mis_ml[k]:
            if sp3flag ==1 or sp2flag ==1 or sp1flag ==1:
                points = str(info[i][1]) + ' ' + str(info[i][2])

                url = 'http://vizier.u-strasbg.fr/viz-bin/votable?-source='+ urllib.parse.quote(cats[k]) +'&-c='+ urllib.parse.quote(points) + '&-c.rs=' + str(rad) + '&-out.all'
                req = requests.get(url,allow_redirects=True, verify=False, headers=headers)


                #put photometry in a votable and extract info
                if flag == '1':
                    xmlfile = subdir + "/" + str(obj[i]) + "_vizier_" + str(cats[k])+"_phot.xml"
                else:
                    xmlfile = dir + "/" + "phot.xml"
                p = open(xmlfile, "wb")
                p.write(req.content)
                with open(xmlfile,'r') as x:
                    lines = x.readlines()
                    for line in lines:
                        if 'COOSYS' in line:
                            exfg =1
                if exfg ==1:
                    votable = parse(xmlfile)
                    table = votable.get_first_table()
                    data = table.array
                    p.close()
                else:
                    data = []

        else:
            points = str(info[i][1]) + ' ' + str(info[i][2])
            print(cats[k])

            url = 'http://vizier.u-strasbg.fr/viz-bin/votable?-source='+ urllib.parse.quote(cats[k]) +'&-c='+ urllib.parse.quote(points) + '&-c.rs=' + str(rad) + '&-out.all'
            print(url)
            req = requests.get(url,allow_redirects=True, verify=False, headers=headers)


            #put photometry in a votable and extract info
            if flag == '1':
                xmlfile = subdir + "/" + str(obj[i]) + "_vizier_" + str(cats[k])+"_phot.xml"
            else:
                xmlfile = dir + "/" + "phot.xml"
            p = open(xmlfile, "wb")
            p.write(req.content)
            if url == 'http://vizier.u-strasbg.fr/viz-bin/votable?-source=J/A%2BA/569/A125/tablea1&-c=342.372854%2059.915645&-c.rs=5&-out.all' or url=='http://vizier.u-strasbg.fr/viz-bin/votable?-source=J/A%2BA/569/A125/tablea1&-c=342.3708333333333%2059.91555555555556&-c.rs=5&-out.all':
                data = []
            else:
                votable = parse(xmlfile)
                table = votable.get_first_table()
                data = table.array
                p.close()
        #mollin #3
        if ":=1612" in mis_ml[k] and 'J/A+AS/90/327/catalog' in cats[k]:
            if len(data) ==0 or mis_ml[k] not in info[i]:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
            elif mis_ml[k] in info[i]:
                peaks = []
                vl = []
                vh = []

                for j in range(len(data)):
                    x = str(data[j]).split(", ")
                    peaks.append(x[5])
                    vl.append(x[8])
                    vh.append(x[9])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(peaks[0])
                row.append(vl[0])
                row.append(vh[0])
        elif 'J/AZh/79/328/table2' in cats[k]:
            if len(data) ==0:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:


                x = str(data[0]).split(", ")
                vlsr1 = (x[3])
                vlsr2 = (x[6])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vlsr1)
                row.append(vlsr2)
        elif 'J/A+A/399/1083/table1' in cats[k]:
            if len(data) ==0:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:


                x = str(data[0]).split(", ")
                vlsr = (x[-10])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vlsr)
        elif 'J/A+A/569/A125/tablea1' in cats[k]:
            print(str(info[i][0]))
            #print(data)
            if len(data) ==0 or mis_ml[k] not in info[i]:
                #print('y')
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            else:#if mis_ml[k] in info[i]:
                print('y')

                x = str(data[0]).split(", ")
                vlsr = (x[7])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vlsr)
        elif 'J/AN/333/634/tablea1' in cats[k]:
            if len(data) ==0:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:


                x = str(data[0]).split(", ")
                vlsr = (x[11])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vlsr)

        elif "J/ApJ/702/1615/masers" in cats[k]:
            if len(data) ==0:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:


                x = str(data[0]).split(", ")
                vlsr = (x[-1][:-1])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vlsr)

        elif ":=1.1" in mis_ml[k] and 'J/ApJ/799/29/table5' in cats[k]:
            if len(data) ==0 or mis_ml[k] not in info[i]:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:
                vlsr = []

                for j in range(len(data)):
                    x = str(data[j]).split(", ")
                    vlsr.append(x[7])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vlsr[0])
        elif 'J/A+A/325/255/table1' in cats[k]:
            if len(data) ==0:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:


                x = str(data[0]).split(", ")
                vl = (x[6])
                vh = x[7]

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vl)
                row.append(vh)
        elif 'J/MNRAS/439/2584/table1' in cats[k]:
            if len(data) ==0:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:


                x = str(data[0]).split(", ")
                vl = (x[9])
                vh = x[10]
                vl2 =x[13]
                vh2 =x[14]

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vl)
                row.append(vh)
                row.append(vl2)
                row.append(vh2)
        elif ":=6668" in mis_ml[k] and 'VIII/96/catalog' in cats[k]:
            if len(data) ==0 or mis_ml[k] not in info[i]:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:
                vl = []
                vh = []

                for j in range(len(data)):
                    x = str(data[j]).split(", ")
                    vl.append(x[7])
                    vh.append(x[8])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vl[0])
                row.append(vh[0])
        elif 'J/A+A/368/845/table2' in cats[k]:
            if len(data) ==0 or mis_ml[k] not in info[i]:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:
                vl = []
                vh = []
                vpeak = []

                for j in range(len(data)):
                    x = str(data[j]).split(", ")
                    vl.append(x[12])
                    vh.append(x[13])
                    vpeak.append(x[14])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vl[0])
                row.append(vh[0])
                row.append(vpeak[0])
        elif ":=6.7" in mis_ml[k] and 'J/AJ/152/92/table1' in cats[k]:
            if len(data) ==0 or mis_ml[k] not in info[i]:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:
                vpeak = []
                dist = []
                lum = []

                for j in range(len(data)):
                    x = str(data[j]).split(", ")
                    vpeak.append(x[8])
                    dist.append(x[10])
                    lum.append(x[11])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vpeak[0])
                row.append(dist[0])
                row.append(lum[0])
        elif ":=36.169" in mis_ml[k] and 'J/ApJS/227/10/table2' in cats[k]:
            if len(data) ==0 or mis_ml[k] not in info[i]:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

            elif mis_ml[k] in info[i]:
                vel = []

                for j in range(len(data)):
                    x = str(data[j]).split(", ")
                    vel.append(x[11])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(vel[0])

        #SPLASH
        elif 'SPLASH' in mis_ml[k] and sp3flag ==1 or sp2flag==1 or sp1flag ==1:
            if len(data) ==0:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
            else:
                ifl = []
                pkvel = []
                minvel = []
                maxvel = []
                assoc = []


                x = str(data[0]).split(", ")
                ifl.append(x[7])
                pkvel.append(x[8])
                minvel.append(x[9])
                maxvel.append(x[10])
                assoc.append(x[11][1:])

                index = info[i].index(mis_ml[k])
                row.append(info[i][index+1]) #catalogue
                row.append(info[i][index+2]) #ra
                row.append(info[i][index+3]) #dec
                row.append(info[i][index+4]) #offset
                row.append(info[i][index+5]) #wav
                row.append(info[i][index+6]) #f
                row.append(info[i][index+7]) #ef
                row.append(ifl[0])
                row.append(pkvel[0])
                row.append(minvel[0])
                row.append(maxvel[0])
                row.append(assoc[0])


    #MOLECULAR LINE DATA
    #h20
    molin = path + '/' + str(root) + "_h2o_data.csv"

    if os.path.exists(molin):
        object = []
        riha = []
        deca = []
        angd = []
        fpeak = []
        fpe = []
        vpeak = []
        vpe = []
        with open(molin, 'r') as m:
            lines = m.readlines()[2:]
            for line in lines:
                x = line.split(",")
                if "Yes" in x[4]:
                    object.append(str(x[1][1:-1]))
                    if x[5][1] == "+":
                        r = x[5][2:-1]
                    else:
                        r = x[5][2:-1]
                    d = x[6][1:-1]
                    coord = r +" "+ d
                    ra, dec = pyasl.coordsSexaToDeg(coord)
                    riha.append(ra)
                    deca.append(dec)
                    angd.append(x[7][1:-1])
                    if x[8] != "":
                        if "(" not in x[8]:
                            fpeak.append(x[8][1:-1])
                            fpe.append(-999)
                        else:
                            s = x[8][1:-1].find("(")
                            e = x[8][1:-1].find(")")
                            fpeak.append(x[8][1:s+1] + x[8][e+2:-1])
                            fpe.append(float(x[8][s+2:e+1]))
                    else:
                        fpeak.append(-999)
                        fpe.append(-999)
                    if x[9] != "":
                        if "(" not in x[9]:
                            vpeak.append(x[9][1:-1])
                            vpe.append(-999)
                        else:
                            s = x[9][1:-1].find("(")
                            e = x[9][1:-1].find(")")
                            vpeak.append(float(x[9][1:s+1]))
                            vpe.append(float(x[9][s+2:e+1]))
                    else:
                        vpeak.append(-999)
                        vpe.append(-999)
        m.close()

        for k in range(len(object)):
            if object[k] == info[i][0]:

                row.append(riha[k])
                row.append(deca[k])
                row.append(angd[k])
                row.append(fpeak[k])
                row.append(fpe[k])
                row.append(vpeak[k])
                row.append(vpe[k])
            else:

                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

    #oh
    molin = path + '/' + str(root) + "_oh_data.csv"

    if os.path.exists(molin):
        object = []
        riha = []
        deca = []
        angd = []
        freq = []
        fpeak = []
        fpe = []
        vpeak = []
        vpe = []
        vpeak2 = []
        vpe2 = []
        with open(molin, 'r') as m:
            lines = m.readlines()[2:]
            for line in lines:
                x = line.split(",")
                if "Yes" in x[4]:
                    object.append(str(x[1][1:-1]))
                    if x[5][1] == "+":
                        r = x[5][2:-1]
                    else:
                        r = x[5][2:-1]
                    d = x[6][1:-1]
                    coord = r +" "+ d
                    ra, dec = pyasl.coordsSexaToDeg(coord)
                    riha.append(ra)
                    deca.append(dec)
                    angd.append(float(x[7][1:-1]))
                    freq.append(float(x[8][1:-1]))
                    if x[10] != "":
                        if "(" not in x[10]:
                            fpeak.append(x[10][1:-1])
                            fpe.append(-999)
                        else:
                            s = x[10][1:-1].find("(")
                            e = x[10][1:-1].find(")")
                            fpeak.append(x[10][1:s+1] + x[10][e+2:-1])
                            fpe.append(float(x[10][s+2:e+1]))
                    else:
                        fpeak.append(-999)
                        fpe.append(-999)
                    if x[11] != "":
                        if "(" not in x[11]:
                            vpeak.append(x[11][1:-1])
                            vpe.append(-999)
                        else:
                            s = x[11][1:-1].find("(")
                            e = x[11][1:-1].find(")")
                            vpeak.append(float(x[11][1:s+1]))
                            vpe.append(float(x[11][s+2:e+1]))
                    else:
                        vpeak.append(-999)
                        vpe.append(-999)
                    if x[12] != "":
                        if "(" not in x[12]:
                            vpeak2.append(x[12][1:-1])
                            vpe2.append(-999)
                        else:
                            s = x[12][1:-1].find("(")
                            e = x[12][1:-1].find(")")
                            vpeak2.append(float(x[12][1:s+1]))
                            vpe2.append(float(x[12][s+2:e+1]))
                    else:
                        vpeak2.append(-999)
                        vpe2.append(-999)
        m.close()
        for k in range(len(object)):

            if object[k] == info[i][0]:
                row.append(riha[k])
                row.append(deca[k])
                row.append(angd[k])
                row.append(freq[k])
                row.append(fpeak[k])
                row.append(fpe[k])
                row.append(vpeak[k])
                row.append(vpe[k])
                row.append(vpeak2[k])
                row.append(vpe2[k])
            else:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

    #sio
    molin = path + '/' + str(root) + "_sio_data.csv"

    if os.path.exists(molin):
        object = []
        riha = []
        deca = []
        angd = []
        fpeak = []
        fpe = []
        vpeak = []
        vpe = []
        with open(molin, 'r') as m:
            lines = m.readlines()[2:]
            for line in lines:
                x = line.split(",")
                if "Yes" in x[4]:
                    object.append(str(x[1][1:-1]))
                    if x[5][1] == "+":
                        r = x[5][2:-1]
                    else:
                        r = x[5][2:-1]
                    d = x[6][1:-1]
                    coord = r +" "+ d
                    ra, dec = pyasl.coordsSexaToDeg(coord)
                    riha.append(ra)
                    deca.append(dec)
                    angd.append(x[7][1:-1])
                    if x[8] != "":
                        if "(" not in x[8]:
                            fpeak.append(x[8][1:-1])
                            fpe.append(-999)
                        else:
                            s = x[8][1:-1].find("(")
                            e = x[8][1:-1].find(")")
                            fpeak.append(x[8][1:s+1] + x[8][e+2:-1])
                            fpe.append(float(x[8][s+2:e+1]))
                    else:
                        fpeak.append(-999)
                        fpe.append(-999)
                    if x[9] != "":
                        if "(" not in x[9]:
                            vpeak.append(x[9][1:-1])
                            vpe.append(-999)
                        else:
                            s = x[9][1:-1].find("(")
                            e = x[9][1:-1].find(")")
                            vpeak.append(float(x[9][1:s+1]))
                            vpe.append(float(x[9][s+2:e+1]))
                    else:
                        vpeak.append(-999)
                        vpe.append(-999)
        m.close()
        for k in range(len(object)):
            if object[k] == info[i][0]:
                row.append(riha[k])
                row.append(deca[k])
                row.append(angd[k])
                row.append(fpeak[k])
                row.append(fpe[k])
                row.append(vpeak[k])
                row.append(vpe[k])
            else:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

    #ch30h class i
    molin = path + '/' + str(root) + "_met1_data.csv"

    if os.path.exists(molin):
        object = []
        riha = []
        deca = []
        angd = []
        fpeak = []
        fpe = []
        vpeak = []
        vpe = []
        with open(molin, 'r') as m:
            lines = m.readlines()[2:]
            for line in lines:
                x = line.split(",")
                if "Yes" in x[4]:
                    object.append(str(x[1][1:-1]))
                    if x[5][1] == "+":
                        r = x[5][2:-1]
                    else:
                        r = x[5][2:-1]
                    d = x[6][1:-1]
                    coord = r +" "+ d
                    ra, dec = pyasl.coordsSexaToDeg(coord)
                    riha.append(ra)
                    deca.append(dec)
                    angd.append(x[7][1:-1])
                    if x[9] != "":
                        if "(" not in x[9]:
                            fpeak.append(x[9][1:-1])
                            fpe.append(-999)
                        else:
                            s = x[9][1:-1].find("(")
                            e = x[9][1:-1].find(")")
                            fpeak.append(x[9][1:s+1] + x[9][e+2:-1])
                            fpe.append(float(x[9][s+2:e+1]))
                    else:
                        fpeak.append(-999)
                        fpe.append(-999)
                    if x[10] != "":
                        if "(" not in x[10]:
                            vpeak.append(x[10][1:-1])
                            vpe.append(-999)
                        else:
                            s = x[10][1:-1].find("(")
                            e = x[10][1:-1].find(")")
                            vpeak.append(float(x[10][1:s+1]))
                            vpe.append(float(x[10][s+2:e+1]))
                    else:
                        vpeak.append(-999)
                        vpe.append(-999)
        m.close()
        for k in range(len(object)):
            if object[k] == info[i][0]:
                row.append(riha[k])
                row.append(deca[k])
                row.append(angd[k])
                row.append(fpeak[k])
                row.append(fpe[k])
                row.append(vpeak[k])
                row.append(vpe[k])
            else:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)

    #ch30h class ii
    molin = path + '/' + str(root) + "_met2_data.csv"
    count = 0
    if os.path.exists(molin):
        object = []
        riha = []
        deca = []
        angd = []
        fpeak = []
        fpe = []
        vpeak = []
        vpe = []
        with open(molin, 'r') as m:
            lines = m.readlines()[2:]
            for line in lines:
                x = line.split(",")
                if "Yes" in x[4]:
                    object.append(str(x[1][1:-1]))
                    if x[5][1] == "+":
                        r = x[5][2:-1]
                    else:
                        r = x[5][2:-1]
                    d = x[6][1:-1]
                    coord = r +" "+ d
                    ra, dec = pyasl.coordsSexaToDeg(coord)
                    riha.append(ra)
                    deca.append(dec)
                    angd.append(x[7][1:-1])
                    if x[9] != "":
                        if "(" not in x[9]:
                            fpeak.append(x[9][1:-1])
                            fpe.append(-999)
                        else:
                            s = x[9][1:-1].find("(")
                            e = x[9][1:-1].find(")")
                            fpeak.append(x[9][1:s+1] + x[9][e+2:-1])
                            fpe.append(float(x[9][s+2:e+1]))
                    else:
                        fpeak.append(-999)
                        fpe.append(-999)
                    if x[10] != "":
                        if "(" not in x[10]:
                            vpeak.append(x[10][1:-1])
                            vpe.append(-999)
                        else:
                            s = x[10][1:-1].find("(")
                            e = x[10][1:-1].find(")")
                            vpeak.append(float(x[10][1:s+1]))
                            vpe.append(float(x[10][s+2:e+1]))
                    else:
                        vpeak.append(-999)
                        vpe.append(-999)
        m.close()
        for k in range(len(object)):
            if object[k] == info[i][0]:
                row.append(riha[k])
                row.append(deca[k])
                row.append(angd[k])
                row.append(fpeak[k])
                row.append(fpe[k])
                row.append(vpeak[k])
                row.append(vpe[k])
            else:
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)
                row.append(-999)


    #iras lrs spectra
    objc = info[i][0]
    if 'IRAS' in objc:
        print('IRAS LRS Spectra')
        name = urllib.parse.quote(objc[5:])
        url = 'http://cdsarc.u-strasbg.fr/viz-bin/vizExec/getsp?III/197&' + str(name)

        req = requests.get(url,allow_redirects=True, verify=False, headers=headers)

        nam = dir + "/" + objc.replace(" ", "") +"_spec.csv"
        p = open(nam, "wb")
        p.write(req.content)

        fg = 0
        with open(nam, "r") as f:
            lines = f.readlines()
            if len(lines) < 10:
                fg = 1

        if fg == 0:
            lam = []
            flux = []
            hea = ['#wavelength(micron),flux(W/m^2/micron)']
            with open(nam, "r") as f:
                lines = f.readlines()

                for line in lines:
                    line = line.strip()
                    if len(line) > 0 and 'ALIGN=RIGHT' in line:
                        ind1st = line.find('NOWRAP>', 1)
                        ind1en = line.find('</TD><TD', 1)
                        ind2st =line.find('NOWRAP>', ind1st+1)
                        ind2en =line.find('</TD></TR>', ind1en+1)
                        lam.append(float(line[ind1st+7:ind1en]))
                        flux.append(float(line[ind2st+7:ind2en]))
            f.close()
            with open(nam,"w") as f:
                f_writer = csv.writer(f)
                f_writer.writerow(hea)
                for g in range(len(lam)):
                    r = []
                    r.append(lam[g])
                    r.append(flux[g])
                    f_writer.writerow(r)

            f.close()
            row.append(nam)
        else:
            row.append(-999)
            os.remove(nam)

    else:
        row.append(-999)


    #sdss optical spectra
    if os.path.exists(op):
        print('SDSS optical spectra')
        ra = float(info[i][1])
        dec = float(info[i][2])
        fg = 0
        fg2 = 0
        plate = 0
        mjd = 0
        fiberid = 0
        with open(op, "r") as f:
            lines = f.readlines()[2:]
            if len(lines) < 1:
                fg = 1
            else:
                for line in lines:
                    x = line.split(",")
                    if np.isclose(ra, float(x[2]), atol=0.01) and np.isclose(dec, float(x[3]), atol=0.01):
                        plate = x[4]
                        mjd = x[5]
                        fiberid = x[6]
                        fg2 = 1
        f.close()
        if fg == 0 and fg2 ==1:
            lam = []
            flux = []
            hea = ['#wavelength(micron)','flux(10**(-17) erg/cm^s/s/angstrom)']
            url = 'https://dr14.sdss.org/optical/spectrum/view/data/format=csv/spec=lite?plateid=' + plate+ '&mjd='+ mjd+'&fiberid=' + fiberid
            req = requests.get(url,allow_redirects=True, verify=False, headers=headers)
            nam = dir + "/" + objc.replace(" ", "") +"_opticalspec.csv"
            p = open(nam, "wb")
            p.write(req.content)
            with open(nam,"r") as f:
                sp_reader = csv.reader(f, delimiter=',')
                next(sp_reader)
                for row2 in sp_reader:
                    x = row2
                    lam.append(float(x[0])*0.0001)
                    flux.append(float(x[1]))
            f.close()
            with open(nam,"w") as f:
                f_writer = csv.writer(f)
                f_writer.writerow(hea)
                for g in range(len(lam)):
                    r = []
                    r.append(lam[g])
                    r.append(flux[g])
                    f_writer.writerow(r)
            row.append(nam)
        else:
            row.append(-999)


    with open(dictphot, mode='a') as dict:
        dict_writer = csv.writer(dict, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        dict_writer.writerow(row)
dict.close()

reader = csv.DictReader(open(dictphot))

result = {}
for row in reader:
    for column, value in row.items():
        result.setdefault(column, []).append(value)

wav = []
fils = []
ob = []
for x, y in result.items():
    if 'wavelength' in x:
        if 'Lin' != x[:3]:
            for i in range(len(y)):
                if str(y[i]) != '-999':

                    w = y[i]
            wav.append(w)
            inde = x.find('_')
            fils.append(x[:inde])
    if 'object' in x:
        for j in range(len(y)):
            ob.append(y[j])

wav, fils = (list(t) for t in zip(*sorted(zip(wav, fils))))

h = []
h.append('object')
h.append('ra(source)')
h.append('dec(source)')
h.append('type')
for i in range(len(fils)):
    h.append(fils[i]+'_catalogue')
    h.append(fils[i]+'_offset(arcsec)')
    h.append(fils[i]+'_wavelength(micron)')
    h.append(fils[i]+'_flux(Jy)')
    h.append(fils[i]+'_error_flux(Jy)')

ord = dir + "/" + root + "_phot_order.csv"

with open(ord, mode='w') as dict:
    dict_writer = csv.writer(dict, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    dict_writer.writerow(h)
dict.close()

for i in range(len(ob)):
    row = []
    row.append(result['object'][i])
    row.append(result['ra(source)'][i])
    row.append(result['dec(source)'][i])
    row.append(result['type'][i])
    for j in range(len(fils)):
        row.append(result[fils[j]+'_catalogue'][i])
        row.append(result[fils[j]+'_offset(arcsec)'][i])
        row.append(result[fils[j]+'_wavelength(micron)'][i])
        row.append(result[fils[j]+'_flux(Jy)'][i])
        row.append(result[fils[j]+'_error_flux(Jy)'][i])

    with open(ord, mode='a') as dict:
        dict_writer = csv.writer(dict, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        dict_writer.writerow(row)
    dict.close()


if flag != '1' and os.path.exists(dir + "/" + "phot.xml"):
    os.remove(dir + "/" + "phot.xml")
    
#SED-RAW    
sed_raw = dir + "/sed-for-MoD"

os.makedirs(sed_raw, exist_ok=True)

#Section 1: This section of SED-RAW code outputs the first raw. It goes to the phot_order.csv and extracts the filter name, flux, flux error, wavelength, VizieR reference catalog, and offset. This part also converts the fluxes and flux errors to milliJansky from Jansky. The name of the object, coordinates (deg), and headers are also placed at the top of the file. This file (and all of the other sed outputs) will be sorted in terms of ascending wavelength.
file = np.loadtxt(ord, delimiter=",", dtype=str)


objects=[]
ra = []
dec = []

for i in range(len(file[1:])):
    objects.append(file[i+1][0])
    ra.append(file[i+1][1])
    dec.append(file[i+1][2])
    

data = file.T[4:]

flux = []
flux_error = []
wavelength = []
filters = []
catalog = []
offset = []

i=0
while i<len(data):
    catalog_name = data[i][0]
    catalog_name = catalog_name.replace("_catalogue", "")
    catalog_name = catalog_name.replace("'", "")
    filters.append(catalog_name)
    wavelength.append(data[i + 2][1:])
    flux.append(data[i + 3][1:])
    flux_error.append(data[i + 4][1:])
    catalog.append(data[i][1:])
    offset.append(data[i + 1][1:])
    i=i+5


flux = np.asarray(flux, dtype=float)

new_flux_error = [[]]

for l in range(len(flux_error)):
    for m in range(len(flux_error[l])):
        if (flux_error[l][m] != '--'):
            new_flux_error[l].append(float(flux_error[l][m]))
        else:
            new_flux_error[l].append(flux_error[l][m])
    new_flux_error.append([])
name = []
for c in range(len(objects)):
    name.append(objects[c].replace(" ", "_"))
    
obj_duplicates = [o for o,t in Counter(name).items() if t>1]
for q in range(len(obj_duplicates)):
    the_duplicates = (np.array(name) == obj_duplicates[q])
    for r in range(len(np.array(name)[the_duplicates])):
        first_occurrence = name.index(obj_duplicates[q])
        name[first_occurrence] = name[first_occurrence] + '_' + str(r+1)
for x in range(len(objects)):
    master_flux = []
    master_flux_error = []
    master_wavelength = []
    master_filter = []
    master_catalog = []
    master_offset = []
    for y in range(len(flux)):
        if (flux[y][x] != -999):
            master_filter.append(filters[y])
            master_flux.append(flux[y][x]*1000)
            if (flux_error[y][x] != '--'):
                master_flux_error.append(new_flux_error[y][x]*1000)
            else:
                master_flux_error.append(new_flux_error[y][x])
            master_wavelength.append(wavelength[y][x])
            master_catalog.append(catalog[y][x])
            master_offset.append(offset[y][x])
            
    phot = list(zip(master_filter, master_flux, master_flux_error, master_wavelength, master_offset, master_catalog))
    #top = objects[x] + ','+ra[x]+','+dec[x]
    top= 'Filter, Flux(mJy), Err(mJy), Wavelength(micron), Offset(arcsec), VizierCatalog'
    f_out = (sed_raw + '/sed2_' + name[x] + '_raw1.csv')
    np.savetxt(f_out, phot, fmt='%s', delimiter=", ", header=top)
    call('sort -nk4 -t"," -o' + f_out + ' ' + f_out, shell=True)
    call('sed -i \'\' -e \'1 i\\' + '\n' +
         '#' + objects[x] + ', ' + ra[x] + ', ' + dec[x] + '\' ' + f_out, shell=True)
#Section 1 finished.

#Section 2: This section of SED-RAW outputs the second raw. This part takes the filter name and replaces it with the appropriate name given by the MoD manual. Some filters have names that are not in the manual, in which case the line is preceded by a # to comment it out, and other filters may have issues, in which case they are also commented out. For Gaia, we want to use Gaia3 as much as possible. If the filter is not present phot_order, then we use Gaia2 and GaiaG.
    filter_name, fluxes, flux_errors, wavelengths, offsets, catalogs = np.loadtxt(f_out, delimiter=",", dtype=str, unpack=True)
    for b in range(len(filter_name)):
        if 'VISTA' in filter_name[b]:
            filter_name[b] = '#' + filter_name[b]
#             if 'Y' in filter_name[b]:
#                 filter_name[b] = '" vistaY"'
#             elif 'J' in filter_name[b]:
#                 filter_name[b] = '" vistaJ"'
#             elif 'Z' in filter_name[b]:
#                 filter_name[b] = '" vistaZ"'
#             elif 'K' in filter_name[b]:
#                 filter_name[b] = '" vistaK"'
#             elif 'H' in filter_name[b]:
#                 filter_name[b] = '" vistaH"'
#             else:
#                 filter_name[b] = '#' + filter_name[b]
        elif 'Johnson' in filter_name[b]:
            if 'B' in filter_name[b]:
                filter_name[b] = '"  BesB"'
            elif 'V' in filter_name[b]:
                filter_name[b] = '"  BesV"'
            elif 'U' in filter_name[b]:
                filter_name[b] = '"  BesU"'
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'IRAS' in filter_name[b]:
            if '100' in filter_name[b]:
                filter_name[b] = '"ira100"'
                if '-' in fluxes[b]:
                    fluxes[b] = fluxes[b].replace("-", "")
                    filter_name[b] = '#' + filter_name[b]
            elif '25' in filter_name[b]:
                filter_name[b] = '"iras25"'
                if '-' in fluxes[b]:
                    fluxes[b] = fluxes[b].replace("-", "")
                    filter_name[b] = '#' + filter_name[b]
            elif '60' in filter_name[b]:
                filter_name[b] = '"iras60"'
                if '-' in fluxes[b]:
                    fluxes[b] = fluxes[b].replace("-", "")
                    filter_name[b] = '#' + filter_name[b]
            elif '12' in filter_name[b]:
                filter_name[b] = '"iras12"'
                if '-' in fluxes[b]:
                    fluxes[b] = fluxes[b].replace("-", "")
                    filter_name[b] = '#' + filter_name[b]
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'WISE' in filter_name[b]:
            if 'unWISE' not in filter_name[b]: 
                if'ALLWISE' in filter_name[b]:
                    if ':1' in filter_name[b]:
                        filter_name[b] = '" WISE1"'
                    elif ':2' in filter_name[b]:
                        filter_name[b] = '" WISE2"'
                    elif ':3' in filter_name[b]:
                        filter_name[b] = '" WISE3"'
                    elif ':4' in filter_name[b]:
                        filter_name[b] = '" WISE4"'
                    else:
                        filter_name[b] = '#' + filter_name[b]
                else:
                    if 'W1' in filter_name[b]:
                        filter_name[b] = '" WISE1"'
                    elif 'W2' in filter_name[b]:
                        filter_name[b] = '" WISE2"'
                    elif 'W3' in filter_name[b]:
                        filter_name[b] = '" WISE3"'
                    elif 'W4' in filter_name[b]:
                        filter_name[b] = '" WISE4"'
                    else:
                        filter_name[b] = '#' + filter_name[b]
            else:
                filter_name[b] = '#' + filter_name[b]
        elif '2MASS' in filter_name[b]:
            if 'J' in filter_name[b]:
                filter_name[b] = '"2massJ"'
            elif 'H' in filter_name[b]:
                filter_name[b] = '"2massH"'
            elif 'K' in filter_name[b]:
                filter_name[b] = '"2massK"'
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'AKARI' in filter_name[b]:
            if '9W' in filter_name[b]:
                filter_name[b] = '"AkaS9W"'
            elif 'L18W' in filter_name[b]:
                filter_name[b] = '"AkL18W"'
            elif 'N60' in filter_name[b]:
                filter_name[b] = '"AkaN60"'
            elif 'WIDE-S' in filter_name[b]:
                filter_name[b] = '" AkaWS"'
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'GALEX' in filter_name[b]:
            if 'NUV' in filter_name[b]:
                filter_name[b] = '"GALNUV"'
            elif 'FUV' in filter_name[b]:
                filter_name[b] = '"GALFUV"'
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'HIP' in filter_name[b]:
            if 'V' in filter_name[b]:
                filter_name[b] = '"Hippac"'
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'GAIA' in filter_name[b] or 'Gaia' in filter_name[b]:
            #filter_name[b] = '#' + filter_name[b]
            if ':G' in filter_name[b]:
                if 'Grp' in filter_name[b]:
                    if 'GAIA/GAIA3:Grp' in filter_name or '"GAIARp"' in filter_name:
                        if "GAIA/GAIA3:Grp" in filter_name[b]:
                            filter_name[b] = '"GAIARp"'
                        else:
                            filter_name[b] = '#' + filter_name[b]
                    elif 'GAIA/GAIA2:Grp' in filter_name or '"GAIARp"' in filter_name:
                        if 'GAIA/GAIA2:Grp' in filter_name[b]:
                            filter_name[b] = '"GAIARp"'
                        else:
                            filter_name[b] = '#' + filter_name[b]
                    else:
                        filter_name[b] = '#' + filter_name[b]
                elif 'Gbp' in filter_name[b]:
                    if 'GAIA/GAIA3:Gbp' in filter_name or '"GAIABp"' in filter_name:
                        if "GAIA/GAIA3:Gbp" in filter_name[b]:
                            filter_name[b] = '"GAIABp"'
                        else:
                            filter_name[b] = '#' + filter_name[b]  
                    elif 'GAIA/GAIA2:Gbp' in filter_name or '"GAIAGBp"' in filter_name:
                        if 'GAIA/GAIA2:Gbp' in filter_name[b]:
                            filter_name[b] = '"GAIABp"'
                        else:
                            filter_name[b] = '#' + filter_name[b]
                    else:
                        filter_name[b] = '#' + filter_name[b]
                elif ':G' in filter_name[b]:
                    if 'GAIA/GAIA3:G' in filter_name or '" GAIAG"' in filter_name:
                        if 'GAIA/GAIA3:G' in filter_name[b]:
                            filter_name[b] = '" GAIAG"'
                        else:
                            filter_name[b] = '#' + filter_name[b]
                    elif 'GAIA/GAIA2:G' in filter_name or '" GAIAG"' in filter_name:
                        if 'GAIA/GAIA2:G' in filter_name[b]:
                            filter_name[b] = '" GAIAG"'
                        else:
                            filter_name[b] = '#' + filter_name[b]
                    else:
                        filter_name[b] = '#' + filter_name[b]
                else:
                    filter_name[b] = '#' + filter_name[b]
            else:
                filter_name[b] = '#' + filter_name[b]
#         elif 'Gaia' in filter_name[b]:
#             filter_name[b] = '#' + filter_name[b]
#             if ':G' in filter_name[b]:
#                 filter_name[b] = '" GAIAG"'
#             else:
#                 filter_name[b] = '#' + filter_name[b]
        elif 'DIRBE' in filter_name[b]:
            if '1' in filter_name[b]:
                filter_name[b] = '" COBE1"'
            elif '2' in filter_name[b]:
                filter_name[b] = '" COBE2"'
            elif '3' in filter_name[b]:
                filter_name[b] = '" COBE3"'
            elif '4' in filter_name[b]:
                filter_name[b] = '" COBE4"'
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'Cousins' in filter_name[b]:
            if 'I' in filter_name[b]:
                filter_name[b] = '" CousI"'
            elif 'R' in filter_name[b]:
                filter_name[b] = '" CousR"'
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'glimpse' in filter_name[b]:
            if '3.6' in filter_name[b]:
                filter_name[b] = '"irac36"'
            elif '4.5' in filter_name[b]:
                filter_name[b] = '"irac45"'
            elif '5.8' in filter_name[b]:
                filter_name[b] = '"irac58"'
            elif '8.0' in filter_name[b]:
                filter_name[b] = '"irac80"'
            else:
                filter_name[b] = '#' + filter_name[b]
        elif 'MSX' in filter_name[b]:
            if 'A' in filter_name[b]:
                filter_name[b] = '"  MSXA"'
                if '-' in fluxes[b]:
                    fluxes[b] = fluxes[b].replace("-", "")
                    filter_name[b] = '#' + filter_name[b]
            else:
                filter_name[b] = '#' + filter_name[b]
        else:
            filter_name[b] = '#' + filter_name[b]
        
    phot = list(zip(filter_name, fluxes, flux_errors, wavelengths, offsets, catalogs))
    f_out = (sed_raw + '/sed2_' + name[x] + '_raw2.csv')
    np.savetxt(f_out, phot, fmt='%s', delimiter=", ", header=top)

    call('sed -i \'\' -e \'1 i\\' + '\n' +
         '#' + objects[x] + ', ' + ra[x] + ', ' + dec[x] + '\' ' + f_out, shell=True)
#Section 2 finished.

#Section 3: This section outputs the third raw. Duplicate lines that are identical in filter name, flux, flux error, and wavelength will be deleted.
    call('sort -nk5 -t"," ' + f_out + ' | awk \'BEGIN{FS=", ";OFS=", "} !seen[$1,$2,$3,$4]++\' | sort -nk4 -t"," > ' + sed_raw + '/sed2_' + name[x] + '_raw3.csv', shell=True)
    
    #call('awk \'BEGIN{FS=", ";OFS=", "} !seen[$1,$2,$3,$4]++\' ' + f_out + ' > ' + sed_raw + '/sed2_' + objects[x] + '_raw3.csv', shell=True)
    
    f_out = (sed_raw + '/sed2_' + name[x] + '_raw3.csv')
    
    call('sed -i \'\' -e \'2d\' ' + f_out, shell=True)
    
    call('sed -i \'\' -e \'1 i\\' + '\n' +
         '#' + objects[x] + ', ' + ra[x] + ', ' + dec[x] + '\' ' + f_out, shell=True)
#Section 3 finished.

#Section 4: This section outputs the fourth raw. For lines that share the same filter, we will take the median average (weighted median average) of the flux. We will be using the flux errors as the weights. First, we need to exclude any lines that have flux that is more or less than 3 standard deviations from the center value, iterated over 10 times. If some lines don't have flux errors (denoted by --) we need to assign it a flux error. If this line is a part of a group of lines that share filter names, for each line that DOES have a flux error, we do flux error/flux for each one and then take the median of them (fractional error). We then take this fractional error and multiply by the flux to get the guess flux error. If this is just a single line that doesn't have flux error, we assume the fractional error to be a value (iras has certain guess fractional errors while others have 0.3) and multiply it by the flux to get a flux error. Lines that were commented out are not affected by this step.
    #taking median averages
    filter_name2, flux2, flux_error2, wavelength2, offset2, catalog2 = np.loadtxt(f_out, delimiter = ", ", dtype=str, unpack=True)

    res = []
    for i in filter_name2:
        if i not in res:
            res.append(i)


    flux2 = np.asarray(flux2, dtype=np.float)
    wavelength2 = np.asarray(wavelength2, dtype=np.float)

    master = []

    def return_indices_of_a(a, b):
        b_set = set(b)
        return [i for i, v in enumerate(a) if v in b_set]


    original_error = []
    for i in range(len(res)):
        frac_errors = []
        name_check = (filter_name2 == res[i])
        flux_mini = flux2[name_check]
        flux_error_mini = flux_error2[name_check]
        wavelength_mini = wavelength2[name_check]
        offset_mini = offset2[name_check]
        catalog_mini = catalog2[name_check]
        clip_flux = sigma_clip(flux_mini, sigma=3, maxiters=10, masked=False)
        indices = return_indices_of_a(flux_mini, clip_flux)

        new_flux_mini = []
        new_flux_error_mini = []
        new_wavelength_mini = []
        new_offset_mini = []
        new_catalog_mini = []

        for s in indices:
            new_flux_mini.append(flux_mini[s])
            new_flux_error_mini.append(flux_error_mini[s])
            new_wavelength_mini.append(wavelength_mini[s])
            new_offset_mini.append(offset_mini[s])
            new_catalog_mini.append(catalog_mini[s])

        flux_mini = np.array(new_flux_mini)
        flux_error_mini = np.array(new_flux_error_mini)
        wavelength_mini = np.array(new_wavelength_mini)
        offset_mini = np.array(new_offset_mini)
        catalog_mini = np.array(new_catalog_mini)

        if (len(flux_mini) == 1):
            if flux_error_mini[0] == ' --':
                flux_error_mini = list(flux_error_mini)
                if res[i] == '"iras12"':
                    frac_error = 0.15
                elif res[i] == '"iras25"':
                    frac_error = 0.15
                elif res[i] == '"iras60"':
                    frac_error = 0.20
                else:
                    frac_error = 0.3
                flux_error_mini[0] = frac_error * flux_mini[0] 
                original_error.append(' --')
            else:
                original_error.append(flux_error_mini[0])

            master.append([res[i], flux_mini[0], flux_error_mini[0], wavelength_mini[0], offset_mini[0], catalog_mini[0]])
        else:
            if len(set(flux_mini)) != len(flux_mini) and len(set(wavelength_mini)) != len(wavelength_mini):
                duplicates = [k for k,v in Counter(flux_mini).items() if v>1]
                trim_flux_mini = []
                trim_flux_error_mini = []
                trim_wavelength_mini = []
                trim_offset_mini = []
                trim_catalog_mini = []
                for d in range(len(duplicates)):
                    duplicate_check = (flux_mini == duplicates[d])
                    if ' --' in flux_error_mini[duplicate_check]:
                        dash = (flux_error_mini[duplicate_check] != ' --')
                        delete_flux = flux_mini[duplicate_check][dash]
                        delete_flux_error = flux_error_mini[duplicate_check][dash]
                        delete_wavelength = wavelength_mini[duplicate_check][dash]
                        delete_offset = offset_mini[duplicate_check][dash]
                        delete_catalog = catalog_mini[duplicate_check][dash]
                        #print(duplicates)
                        trim_flux_error_mini.append(delete_flux_error)
                        trim_flux_mini.append(delete_flux)
                        trim_wavelength_mini.append(delete_wavelength)
                        trim_offset_mini.append(delete_offset)
                        trim_catalog_mini.append(delete_catalog)
                    
                    else:
                        trim_flux_error_mini.append(flux_error_mini[duplicate_check])
                        trim_flux_mini.append(flux_mini[duplicate_check])
                        trim_wavelength_mini.append(wavelength_mini[duplicate_check])
                        trim_offset_mini.append(offset_mini[duplicate_check])
                        trim_catalog_mini.append(catalog_mini[duplicate_check])

#                     flux_temp = np.array(trim_flux_mini).flatten()
#                     flux_error_temp = np.array(trim_flux_error_mini).flatten()
#                     wavelength_temp = np.array(trim_wavelength_mini).flatten()
#                     offset_temp = np.array(trim_offset_mini).flatten()
#                     catalog_temp = np.array(trim_catalog_mini).flatten()

                    flux_temp = np.concatenate(trim_flux_mini).ravel()
                    flux_error_temp = np.concatenate(trim_flux_error_mini).ravel()
                    wavelength_temp = np.concatenate(trim_wavelength_mini).ravel()
                    offset_temp = np.concatenate(trim_offset_mini).ravel()
                    catalog_temp = np.concatenate(trim_catalog_mini).ravel()
                for y in range(len(flux_mini)):
                    if clip_flux[y] not in flux_temp:
                        flux_temp = np.append(flux_temp, flux_mini[y])
                        flux_error_temp = np.append(flux_error_temp, flux_error_mini[y])
                        wavelength_temp = np.append(wavelength_temp, wavelength_mini[y])
                        offset_temp = np.append(offset_temp, offset_mini[y])
                        catalog_temp = np.append(catalog_temp, catalog_mini[y])

                flux_mini = flux_temp
                flux_error_mini = flux_error_temp
                wavelength_mini = wavelength_temp
                offset_mini = offset_temp
                catalog_mini = catalog_temp

            if (np.all(flux_error_mini == ' --') == False):
                if ' --' in flux_error_mini:
                    for j in range(len(flux_error_mini)):
                        if flux_error_mini[j] != ' --':
                            frac_errors.append(float(flux_error_mini[j]) / float(flux_mini[j]))

                    median_errors2 = ws.median(frac_errors)
                    
                    flux_error_mini = list(flux_error_mini)

                    for j in range(len(flux_error_mini)):
                        if flux_error_mini[j] == ' --':
                            flux_error_mini[j] = median_errors2 * flux_mini[j]

                clip_flux_error = np.asarray(flux_error_mini, dtype=float)

                clip_flux_error = list(clip_flux_error)
                clip_flux = list(flux_mini)


                clip_flux_median = ws.weighted_median(clip_flux, weights=clip_flux_error)


                if clip_flux_median in clip_flux:

                    median_check = (flux_mini == clip_flux_median)
                    
                    original_error.append(np.array(clip_flux_error)[median_check][0])

                    master.append([res[i], clip_flux_median, np.array(clip_flux_error)[median_check][0], np.array(wavelength_mini)[median_check][0], np.array(offset_mini)[median_check][0], np.array(catalog_mini)[median_check][0]])
                else:
                    offset_final = ';'.join(str(p) for p in offset_mini)
                    catalog_final = ';'.join(str(p) for p in catalog_mini)
                    
                    if np.sum(clip_flux_error) == 0:
                        flux_aver = np.average(clip_flux)
                        err = 0
                    else:
                        flux_aver = np.average(clip_flux,weights=clip_flux_error)
                        err = flux_aver * np.sqrt(((clip_flux_error[0] / clip_flux[0])**2 + (clip_flux_error[1] / clip_flux[1])**2)/2)
                    
                    original_error.append(err)

                    master.append([res[i], flux_aver, err, wavelength_mini[0], offset_final, catalog_final])
            else:
#                 clip_flux_median = ws.median(flux_mini)
#                 offset_final = ';'.join(str(p) for p in offset_mini)
#                 catalog_final = ';'.join(str(p) for p in catalog_mini)

#                 master.append([res[i], clip_flux_median, ' --', wavelength_mini[0], offset_final, catalog_final])
                if ' --' in flux_error_mini:
                    flux_error_mini = list(flux_error_mini)
                    for j in range(len(flux_error_mini)):
                        if flux_error_mini[j] == ' --':
                            if res[i] == '"iras12"':
                                frac_errors.append(0.15)
                            elif res[i] == '"iras25"':
                                frac_errors.append(0.15)
                            elif res[i] == '"iras60"':
                                frac_errors.append(0.20)
                            else:
                                frac_errors.append(0.3)

                    median_errors2 = ws.median(frac_errors)

                    for j in range(len(flux_error_mini)):
                        if flux_error_mini[j] == ' --':
                            flux_error_mini[j] = median_errors2 * flux_mini[j]

                clip_flux_error = np.asarray(flux_error_mini, dtype=float)

                clip_flux_error = list(clip_flux_error)
                clip_flux = list(flux_mini)


                clip_flux_median = ws.weighted_median(clip_flux, weights=clip_flux_error)
                
                original_error.append(' --')


                if clip_flux_median in clip_flux:

                    median_check = (flux_mini == clip_flux_median)

                    master.append([res[i], clip_flux_median, np.array(clip_flux_error)[median_check][0], np.array(wavelength_mini)[median_check][0], np.array(offset_mini)[median_check][0], np.array(catalog_mini)[median_check][0]])
                else:
                    offset_final = ';'.join(str(p) for p in offset_mini)
                    catalog_final = ';'.join(str(p) for p in catalog_mini)

                    flux_aver = np.average(clip_flux,weights=clip_flux_error)

                    err = flux_aver * np.sqrt(((clip_flux_error[0] / clip_flux[0])**2 + (clip_flux_error[1] / clip_flux[1])**2)/2)

                    master.append([res[i], flux_aver, err, wavelength_mini[0], offset_final, catalog_final])

    rows=[]
    with open(f_out, "r") as fp:
        line = fp.readlines()
        for row in line:
            row = row.split(", ")
            rows.append(row)

    hashtag = []

    for b in range(len(rows)):
        if b !=0 and b != 1:
            if '#' in rows[b][0]:
                master.append(rows[b])
    
    f_out = (sed_raw + '/sed2_' + name[x] + '_raw4.csv')
    np.savetxt(f_out, master, fmt='%s', delimiter=", ", header=top)
    call('sort -nk4 -t"," -o' + f_out + ' ' + f_out, shell=True)
    call('sed -i \'\' -e \'/^$/d\' ' + f_out, shell=True)
    call('sed -i \'\' -e \'1 i\\' + '\n' +
         '#' + objects[x] + ', ' + ra[x] + ', ' + dec[x] + '\' ' + f_out, shell=True)
#Section 4 finished.

#Section 5: This section outputs the fifth raw. Some MoD filters require the flux and flux error be given in magnitude (you can view this in the manual). This section converts the fluxes for the filters that require it from mJy to mag. Each filter has a different zero-point, and this can be viewed in MyDusty_v16a.f (or MyDusty_v16a_noPGP.f) of MoD. This also outputs the original flux and flux errors if the flux and flux error was converted to mag and the nuFnu and nuFnu error.
    filter_name3, flux3, flux_error3, wavelength3, offset3, catalog3 = np.loadtxt(f_out, delimiter = ", ", dtype=str, unpack=True)
    
    total = []
    jy = ['" pacsb"', '" pacsg"', '" pacsr"', '"spire1"', '"spire2"', '"spire3"', '"apx870"', '"ap1200"', '"ukt350"', '"ukt450"', '"ukt800"', '"ukt850"', '"uk1100"', '"uk1300"', '"uk2000"', '"se1300"', '" COBE1"', '" COBE2"', '" COBE3"', '" COBE4"', '" COBE5"']
    
    for i in range(len(filter_name3)):
        if filter_name3[i] not in jy:
            if filter_name3[i] == '"  BesU"': #if the name of the filter equals "  BesU"
                zp = 4.175e-8 #zero-point magnitude for "  BesU"

            elif filter_name3[i] == '"  BesB"':
                zp = 6.255e-8

            elif filter_name3[i] == '"  BesV"':
                zp = 3.627e-8

            elif filter_name3[i] == '"2massJ"':
                zp = 3.137e-9

            elif filter_name3[i] == '"2massH"':
                zp = 1.130e-9

            elif filter_name3[i] == '"2massK"':
                zp = 4.312e-10

            elif filter_name3[i] == '"iras12"':
                l0= 12.
                zp= 40.5 *2.998e-12/l0/l0

            elif filter_name3[i] == '"iras25"':
                l0= 25.
                zp= 8.97 *2.998e-12/l0/l0

            elif filter_name3[i] == '"iras60"':
                l0= 60.
                zp= 1.45 *2.998e-12/l0/l0

            elif filter_name3[i] == '"ira100"':
                l0= 100.
                zp= 0.423 *2.998e-12/l0/l0

            elif filter_name3[i] == '"  MSXA"':
                zp= 2.364e-12

            elif filter_name3[i] == '" WISE1"':
                zp= 8.12e-11

            elif filter_name3[i] == '" WISE2"':
                zp= 2.39e-11

            elif filter_name3[i] == '" WISE3"':
                zp= 6.100e-13

            elif filter_name3[i] == '" WISE4"':
                zp= 4.93e-14

            elif filter_name3[i] == '"Hippac"':
                zp= 3.916e-8

            elif filter_name3[i] == '"GALFUV"':
                zp= 6.74e-8

            elif filter_name3[i] == '"GALNUV"':
                zp= 4.36e-8

            elif filter_name3[i] == '" CousR"':
                zp= 2.137e-8

            elif filter_name3[i] == '" CousI"':
                zp= 1.200e-8

            elif filter_name3[i] == '" GAIAG"':
                zp= 2.170e-8

            elif filter_name3[i] == '"GAIABp"':
                zp= 3.816e-8

            elif filter_name3[i] == '"GAIARp"':
                zp= 1.237e-8
                
            elif filter_name3[i] == '"GAIRVS"':
                zp= 0.9088e-8

            elif filter_name3[i] == '"AkaS9W"':
                zp= 4.49e-12

            elif filter_name3[i] == '"AkL18W"':
                zp= 1.93e-12

            elif filter_name3[i] == '"AkaN60"':
                zp=6.72e-16

            elif filter_name3[i] == '" AkaWS"':
                zp=2.24e-16

            elif filter_name3[i] == '"irac36"':
                zp= 6.63e-11

            elif filter_name3[i] == '"irac45"':
                zp= 2.67e-11

            elif filter_name3[i] == '"irac58"':
                zp= 10.81e-12

            elif filter_name3[i] == '"irac80"':
                zp= 3.08e-12
                
            #converting from mJy to mag and replacing the current flux
            jansky_flux = float(flux3[i])
            
            nufnu_flux = float(flux3[i]) * 2.9979e-12 / (float(wavelength3[i]))
            nufnu_error = float(flux_error3[i]) * 2.9979e-12 / (float(wavelength3[i]))

            if original_error[i] == ' --':
                jansky_error = ' --'
            else:
                jansky_error = float(flux_error3[i])

            f_lambda = float(flux3[i]) * 2.9979e-15/(float(wavelength3[i]) ** 2)
            mag = 2.5 * np.log10(zp/f_lambda)
            flux3[i] = mag

            if flux_error3[i] != " --":
                mag_error = 2.5 * np.log10(1 + float(flux_error3[i])/jansky_flux)
                flux_error3[i] = mag_error
        else:
            jansky_flux = float(flux3[i])
            jansky_error = float(flux_error3[i])
            
            nufnu_flux = float(flux3[i]) * 2.9979e-12 / (float(wavelength3[i]))
            nufnu_error = float(flux_error3[i]) * 2.9979e-12 / (float(wavelength3[i]))
            
            flux3[i] = ' --'
            flux_error3[i] = ' --'
    
        total.append([filter_name3[i], flux3[i], flux_error3[i], jansky_flux, jansky_error, nufnu_flux, nufnu_error, wavelength3[i], offset3[i], catalog3[i]])
    
    rows=[]
    with open(f_out, "r") as fp:
        line = fp.readlines()
        for row in line:
            row = row.split(", ")
            rows.append(row)

    hashtag = []
    hashtag_flux = []
    hashtag_error = []
    hashtag_wavelength = []
    hashtag_offset = []
    hashtag_catalog = []

    for b in range(len(rows)):
        if b !=0 and b != 1:
            if '#' in rows[b][0]:
                hashtag.append(rows[b][0])
                hashtag_flux.append(rows[b][1])
                hashtag_error.append(rows[b][2])
                hashtag_wavelength.append(rows[b][3])
                hashtag_offset.append(rows[b][4])
                hashtag_catalog.append(rows[b][5])
    
    for u in range(len(hashtag)):
        total.append([hashtag[u], ' --', ' --', hashtag_flux[u], hashtag_error[u], ' --', ' --', hashtag_wavelength[u], hashtag_offset[u], hashtag_catalog[u]])
    
    f_out = (sed_raw + '/sed2_' + name[x] + '_raw5.csv')
    #top = objects[x] + ','+ra[x]+','+dec[x]
    top= 'Filter, Flux(mag), Err(mag), Flux(mJy), Err(mJy), nuFnu(erg/cm^2/s), nuFnu(error), Wavelength(micron), Offset(arcsec), VizierCatalog'
    np.savetxt(f_out, total, fmt='%s', delimiter=", ", header=top)
    call('sort -nk8 -t"," -o' + f_out + ' ' + f_out, shell=True)
    call('sed -i \'\' -e \'/^$/d\' ' + f_out, shell=True)
    call('sed -i \'\' -e \'1 i\\' + '\n' +
         '#' + objects[x] + ', ' + ra[x] + ', ' + dec[x] + '\' ' + f_out, shell=True)
#Section 5 finished.

#Section 6. This section outputs the .dat file that will be fed into MoD. We need to take the raw5 and format it in a way that MoD recognizes. We only preserve the filter name, flux, and flux error (mag or mJy whatever is appropriate). We also put #, #, # at the bottom.
    filter_name4, flux_mag, flux_mag_error, flux_jansky, flux_jansky_error, nufnu, nufnu_err, wavelength4, offset4, catalog4 = np.loadtxt(f_out, delimiter = ", ", dtype=str, unpack=True)
    
    final_hopefully = []
    for z in range(len(filter_name4)):
        if flux_mag[z] == ' --':
            final_hopefully.append([filter_name4[z], flux_jansky[z], flux_jansky_error[z]])
        else:
            final_hopefully.append([filter_name4[z], flux_mag[z], flux_mag_error[z]])
    
    f_out = (sed_raw + '/sed2_' + name[x] + '.dat')
    np.savetxt(f_out, final_hopefully, fmt='%s', delimiter=", ")
    call('sed -i \'\' -e \'s/"//g\' ' + f_out, shell=True)
    call('echo "#, #, #," >> ' + f_out, shell=True)

#Section 6 finished.

print(str(len(obj)) + ' objects done in ' "--- %s seconds ---" % (time.time() - start_time))
