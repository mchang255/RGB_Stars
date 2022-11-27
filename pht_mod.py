#This program takes the sed2_name_raw5.csv that was created when get-dict-phot.py was ran and makes it into a LaTeX table (name_pht_MoD.tex). Specifically, the table outputs the following columns in the following order: filter name, wavelength, flux, flux error, offset, and catalog of name
#Run this program in the same directory as your outputs from get-dict-phot.py. In the terminal, type python pht_mod.py name. What I mean by name is that if you look at the sed2 file names for your star, you will notice that there are underscores that replace the spaces. Make sure you write down the exact way the star is written in the file names. Ex: If I want to create a table for the star IRAS 23361+6437, I would type python pht_mod.py IRAS_23361+6437 because that is how the sed2 file name is written in the sed-for-MoD folder in ohir20-2.5asec.

import sys
import glob
import numpy as np
from pathlib import Path
from subprocess import call

name = sys.argv[1]
path = str(Path(__file__).parent.absolute())
i = 1

while True:
    if i < 10:
        file_name = glob.glob(path + '/ohir0' + str(i) +'-2.5asec/sed-for-MoD/sed2_' + name + '_raw5.csv')
    else:
        file_name = glob.glob(path + '/ohir' + str(i) +'-2.5asec/sed-for-MoD/sed2_' + name + '_raw5.csv')
    
    print(file_name)
    if file_name != []:
        break
        
    i = i + 1

if name[-2] == "_":
    name = name[:-2]
    
call('sed \'s/#/%/g\' ' + file_name[0] + ' | awk \'BEGIN{FS=", ";OFS=", "} {print $1,$4,$5,$8,$9,$10}\' > ' + name + '_latex.csv', shell=True)
call('sed -i \'\' -e \'1,2d\' ' + name + '_latex.csv', shell=True)

filter_name, flux, error, wvl, offset, catalog = np.loadtxt(name + '_latex.csv', delimiter = ", ", dtype = str, unpack = True)

sig_fig = [flux, error, wvl, offset]

for f in sig_fig:
    for x in range(len(f)):
        if f[x] != ' --':
            if ';' not in f[x]:
                f[x] = str('%.4f' % float(f[x]))
            else:
                temp_list = f[x].split('; ')
                for a in range(len(temp_list)):
                    temp_list[a] = float(temp_list[a])
                    temp_list[a] = str('%.4f' % temp_list[a])
                f[x] = '; '.join(temp_list)

catalog_new = []
for g in range(len(catalog)):
    new = str(catalog[g]) + ' \\\\'
    catalog_new.append(new)
    
top = "\\begin{deluxetable*}{cccccc}\n"
top += "\\tablecaption{" + name + " SED Input Data\label{tab:" + name.lower() + "}}"
top += "\n\\tabletypesize{\small}"
top += "\n\\tablewidth{0pt}"
top += "\n\\tablehead{\colhead{Filter Name} & \colhead{Wavelength} & \colhead{Flux} & \colhead{Flux Error} & \colhead{Offset} & \colhead{Catalog}\\\\"
top += "\n\colhead{} & \colhead{$\mu$m} & \colhead{mJy} & \colhead{mJy} & \colhead{arcsec} & \colhead{}}"
top += "\n\\startdata"

data = list(zip(filter_name, wvl, flux, error, offset, catalog_new))
np.savetxt(name + '_pht_MoD.tex', data, fmt = "%s", delimiter = " & ", header=top)

call('sed -i \'\' -e \'s/#//g;s/^ //;1,/_/s// /;s/"//g\' ' + name + '_pht_MoD.tex', shell=True)
call('echo "\\enddata" >> ' + name + '_pht_MoD.tex', shell = True )
call('echo "\\end{deluxetable*}" >> ' + name + '_pht_MoD.tex', shell = True )