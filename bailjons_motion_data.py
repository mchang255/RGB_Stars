#This program goes into the file ohir-all-final_Thres_eDR3-bailjons-cln.csv and extracts the name of source, input RA, input dec, output ra, output dec, offset, parallax, parallax error, proper motion, RA proper motion, RA proper motion error, dec proper motion, dec proper motion error, distance, lower limit distance, and upper limit distance of select stars (bailjons_motion_stars.inp) and outputs it into a .tex file formatted as a landscape table (bailjons_table.tex).
#Run this program wherever you want but I did it in the photometry folder as that's where the other .py programs are. Make sure you have the correct file paths for your files. Have your list of stars ready (bailjons_motion_stars.inp). To run this program, type python bailjons_motion_data.py

import sys
import glob
import numpy as np
from pathlib import Path
from subprocess import call

file_path = '/Users/miranda/mchang-JPL2021/ohir-all-final_Thres_eDR3-bailjons-cln.csv'
star_file = 'bailjons_motion_stars.inp'

star_lines = []
with open(star_file, 'r') as star:
    lines = star.readlines()
    for line in lines:
        if line[:1] == '#':
            continue
        else:
            line = line.strip()
            star_lines.append(line)

start = 'grep'
for x in range(len(star_lines)):
    start = start + ' -e "' + star_lines[x] + '"'
    
print(start)
call('awk \'BEGIN{FS=",";OFS=","} {print $1,$2,$3,$5,$6,$4,$16,$17,$19,$20,$21,$22,$23,$128,$129,$130}\' ' + file_path + ' | ' + start + '> new_file.csv', shell=True )

name, ra_input, dec_input, ra_output, dec_output, offset, plx, e_plx, pm, pmra, e_pmra, pmdec, e_pmdec, dist, low_dist, upp_dist = np.loadtxt('new_file.csv', dtype=str, delimiter = ",", unpack=True)

coords = [ra_input, dec_input, ra_output, dec_output]

for f in coords:
    for x in range(len(f)):
        temp_list = f[x].split(' ')
        temp_list[2] = str('%.3f' % float(temp_list[2]))
        f[x] = ' '.join(temp_list)
        
dsn = [dist, low_dist, upp_dist]

for a in dsn:
    for c in range(len(a)):
        a[c] = str('%.3f' % float(a[c]))
        
plxs = [plx, e_plx]

for b in plxs:
    for e in range(len(b)):
        b[e] = str('%.3f' % float(b[e]))
        

for n in range(len(upp_dist)):
    upp_dist[n] = upp_dist[n] + ' \\\\'

data = list(zip(name, ra_input, dec_input, ra_output, dec_output, offset, plx, e_plx, pm, pmra, e_pmra, pmdec, e_pmdec, dist, low_dist, upp_dist))

top = "\\begin{longrotatetable}\n"
top += "\\begin{deluxetable*}{ccccccccccccccccc}\n"
top += "\\tablecaption{Bailer-Jones \& Gaia Data}"
top += "\n\\tablewidth{0pt}"
top += "\n\\tabletypesize{\\tiny}"
top += "\n\\tablehead{\colhead{Name} & \colhead{Input RA} & \colhead{Input Dec} & \colhead{Output RA} & \colhead{Output Dec} & \colhead{Offset} & \colhead{Plx} & \colhead{Plx Err} & \colhead{PM} & \colhead{RA PM} & \colhead{RA PM Err} & \colhead{Dec PM} & \colhead{Dec PM Err} & \colhead{Avg Dist} & \colhead{Low Dist} & \colhead{Upp Dist} \\\\"
top += "\n\colhead{} & \colhead{sexagesimal} & \colhead{sexagesimal} & \colhead{sexagesimal} & \colhead{sexagesimal} & \colhead{arcsec} & \colhead{mas} & \colhead{mas/yr} & \colhead{mas/yr} & \colhead{mas/yr} & \colhead{mas/yr} & \colhead{mas/yr} & \colhead{mas/yr} & \colhead{pc} & \colhead{pc} & \colhead{pc}}"
top += "\n\startdata"

np.savetxt('bailjons_table.tex', data, fmt = "%s", delimiter = " & ", header=top)

call('sed -i \'\' -e \'s/#//g;s/^ //\' bailjons_table.tex', shell=True)
call('echo "\\enddata" >> bailjons_table.tex', shell = True )
call('echo "\\end{deluxetable*}" >> bailjons_table.tex', shell = True )
call('echo "\\end{longrotatetable}" >> bailjons_table.tex', shell = True )
call('rm new_file.csv', shell = True)