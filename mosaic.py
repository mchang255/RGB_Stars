#Program Purpose: This program takes the outputs from the postage_stamp.py program and creates a mosaic of a user-specified size (specify in root-mosaic.inp). We are using GAIA EDR3 coordinates.
#To run this program, type: python mosaic.py root-mosaic.inp

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.wcs.utils import skycoord_to_pixel
from PyAstronomy import pyasl
from regions import CircleSkyRegion

root = sys.argv[1][:-11]
coordfile = root + "-prgb-final-candidates.csv"
inpfile = root + "-mosaic.inp"

#reading in root-mosaic.inp
with open(inpfile, 'r') as inp:
    lines = inp.readlines()
    for line in lines:
        if line[:1] == '#':
            continue
        else:
            index = line.find('=')
            if 'length' in line[:index]:
                length = float(line[index+1:].strip())
            elif 'width' in line[:index]:
                width = float(line[index+1:].strip())
            elif 'circle' in line[:index]:
                circle = float(line[index+1:].strip())
            elif 'c_format' in line[:index]:
                c_format = str(line[index+1:].strip())
            elif 'sep' in line[:index]:
                sep = str(line[index+1:].strip())
            elif 'usr_flg' in line[:index]:
                usr_flg = int(line[index+1:].strip())
            elif 'noplot' in line[:index]:
                noplot = str(line[index+1:].strip())
            elif 'subsetsize' in line[:index]:
                subsetsize = float(line[index+1:].strip())
                
ra = []
dec = []
obj = []
new_names = []

#reading in coordinates from root_coords.csv
with open(coordfile, 'r') as c:
    lines = c.readlines()[1:]
    for line in lines:
        if line[:1] == "#":
            continue
        else:
            x = line.split(sep)
            obj.append(str(x[0].strip()))
            if c_format == 'd':
                ra.append(float(x[2]))
                dec.append(float(x[3]))
            elif c_format =='s':
                c = SkyCoord(str(x[2])+' '+str(x[3]), unit=(u.hourangle, u.deg))
                ra.append(float(c.ra.degree))
                dec.append(float(c.dec.degree))

avg_lum = []
upp_lum = []
low_lum = []

avg_dist = []
upp_dist = []
low_dist = []

blue_vel = []
red_vel = []

#reading luminosities, distances, and OH maser line velocities of each star
with open(coordfile, 'r') as c:
    lines = c.readlines()[1:]
    for line in lines:
        if line[:1] == "#":
            continue
        else:
            x = line.split(',')
            avg_lum.append(float(x[7]))
            low_lum.append(float(x[8]))
            upp_lum.append(float(x[9]))
            
            avg_dist.append(float(x[10]))
            low_dist.append(float(x[11]))
            upp_dist.append(float(x[12]))
            
            if x[24] == '-999':
                blue_vel.append('Unavailable')
            else:
                blue_vel.append(float(x[24]))
                
            if x[25] == '-999\n':
                red_vel.append('Unavailable')
            else:
                red_vel.append(float(x[25]))
                
#replacing spaces in the names of the stars
for b in range(len(obj)):
    new_name = obj[b].replace(" ", "_")
    new_names.append(new_name)

#creating the mosaic
for a in range(len(ra)):
    print('Creating mosaic for ' + new_names[a])
    files = os.listdir(root + '_img/' + new_names[a] + '_imgs') #listing the files in each of the stars' directories
    fig = plt.figure(figsize=(30,30))

    new_filenames = []
    
    #Creating the list of fit files we are going to use in the directory. sdss and iras surveys currently give errors when opening them, so we will skip them.
    for c in range(len(files)):

        if 'sdss' in files[c] or 'iras' in files[c]:
            continue
        elif '.fits' not in files[c]:
            continue
        else:
            new_filenames.append(files[c])
            
    #arranging images by increasing wavelength (left to right, top to bottom)        
    wvl_order = ['dss1blue', 'dss2blue', 'dss1red', 'dss2red', 'dss2ir', '2mass_j', '2mass_h', '2mass_k', 'wise_1', 'irac1', 'irac2', 'wise_2', 'irac3', 'irac4', 'wise_3', 'wise_4', 'mips', 'akari_fisn60', 'akari_fiswides', 'akari_fiswidel', 'akari_fisn160(160)']
    new_filenames_order = []
    for p in range(len(wvl_order)):
        for e in range(len(new_filenames)):
            if wvl_order[p] in new_filenames[e]:
                new_filenames_order.append(new_filenames[e])

    #if we specified that there are images we want to skip, we will skip them            
    if noplot != 'none':
        li = list(noplot.split(","))
        for ab in li:
            if ab == 'AKARI':
                new_filenames_order = [ u for u in new_filenames_order if 'akari' not in u ]
            elif ab == 'SEIP':
                new_filenames_order = [ u for u in new_filenames_order if 'spitzer' not in u ]
            elif ab == 'DSS':
                new_filenames_order = [ u for u in new_filenames_order if 'dss' not in u ]
            elif ab == '2MASS':
                new_filenames_order = [ u for u in new_filenames_order if '2mass' not in u ]
            elif ab == 'WISE':
                new_filenames_order = [ u for u in new_filenames_order if 'wise' not in u ]
    
    sqr_root = math.ceil(np.sqrt(len(new_filenames_order))) #try to make a square shape of the mosaic as much as possible
    print('How many images used in the mosaic:' + str(len(new_filenames_order)))

    mins = []
    maxs = []
    stretch = []
    
    #if we want to set our own min, max intensities and stretch
    if usr_flg == 1:
        min_max_inten_stretch = root + '_img/' + new_names[a] + '_imgs/' + root + '_' + str(ra[a]) + '_' + str(dec[a]) + '-min-max.csv'
        name, mins, maxs, stretch = np.loadtxt(min_max_inten_stretch, delimiter = ", ", dtype = str, unpack=True)
    
    for d in range(len(new_filenames_order)):
        filename = get_pkg_data_filename(root + '_img/' + new_names[a] + '_imgs/' + new_filenames_order[d])

        print(filename)

        hdu = fits.open(filename)[0]
        wcs = WCS(hdu.header)
        
        sky_coord = SkyCoord(ra[a], dec[a], frame='fk5', unit='deg')
        pixel_x,pixel_y = skycoord_to_pixel(sky_coord, wcs)
        
        size1 = len(hdu.data)
        size2 = len(hdu.data[0])
        if size1 > size2:
            img_size = (length * size1) / (subsetsize * 60)
        else:
            img_size = (length * size2) / (subsetsize * 60)
        
        #cropping the image based on the sizes we specified in the .inp
        x_lower_limit = math.ceil(pixel_x - (img_size/2))
        x_upper_limit = math.ceil(pixel_x + (img_size/2))

        y_lower_limit = math.ceil(pixel_y - (img_size/2))
        y_upper_limit = math.ceil(pixel_y + (img_size/2))
        
        if (x_lower_limit < 0):
            difference_x = 0 - x_lower_limit
            x_upper_limit = x_upper_limit + difference_x
            if (y_lower_limit < 0):
                difference_y = 0 - y_lower_limit
                y_upper_limit = y_upper_limit + difference_y
                wcs_cut = wcs[0:y_upper_limit, 0:x_upper_limit]
                hdu_cut = hdu.data[0:y_upper_limit, 0:x_upper_limit]
            else:
                wcs_cut = wcs[y_lower_limit:y_upper_limit, x_lower_limit:x_upper_limit]
                hdu_cut = hdu.data[y_lower_limit:y_upper_limit, x_lower_limit:x_upper_limit]
        elif (y_lower_limit < 0):
            difference_y = 0 - y_lower_limit
            y_upper_limit = y_upper_limit + difference_y
            if (x_lower_limit < 0):
                difference_x = 0 - x_lower_limit
                x_upper_limit = x_upper_limit + difference_x
                wcs_cut = wcs[0:y_upper_limit, 0:x_upper_limit]
                hdu_cut = hdu.data[0:y_upper_limit, 0:x_upper_limit]
            else:
                wcs_cut = wcs[y_lower_limit:y_upper_limit, x_lower_limit:x_upper_limit]
                hdu_cut = hdu.data[y_lower_limit:y_upper_limit, x_lower_limit:x_upper_limit]
        else:
            wcs_cut = wcs[y_lower_limit:y_upper_limit, x_lower_limit:x_upper_limit]
            hdu_cut = hdu.data[y_lower_limit:y_upper_limit, x_lower_limit:x_upper_limit]
            
        ax = fig.add_subplot(sqr_root,sqr_root,d+1, projection=wcs_cut) 
        ax.set_title(new_filenames_order[d], fontsize = 10)
        ax.axis('off')
        
        if usr_flg == 0:
            ax.imshow(hdu_cut)
            #getting intensities from: https://stackoverflow.com/questions/57930402/is-it-possible-to-get-vmin-vmax-from-axes-after-using-imshow?rq=1
            imgs = ax.get_images()
            if len(imgs) > 0:
                vmin, vmax = imgs[0].get_clim()

            mins.append(vmin)
            maxs.append(vmax)
            stretch.append('linear')
        elif usr_flg == 1:
            norm = simple_norm(hdu_cut, stretch[d], min_cut = float(mins[d]), max_cut = float(maxs[d]))
            ax.imshow(hdu_cut, norm=norm)

        #placing circle around GAIA EDR3 coordinates. Akari survey images shrink down greatly when a circle is placed on top, so we will not place a circle on them.
        if 'akari' not in filename:
            sky_radius = Angle((0,0,circle), 'deg')
            sky_region = CircleSkyRegion(sky_coord, sky_radius)
            pixel_region = sky_region.to_pixel(wcs_cut)
            pixel_region.plot(lw = 2, edgecolor='red', facecolor='none')
          
        #adding labels
        if usr_flg == 0:
            caption = 'STRETCH: linear\n MIN: ' + str(vmin) + ' pixels\n MAX: ' + str(vmax) + ' pixels\n'
        elif usr_flg == 1:
            caption = 'STRETCH: ' + stretch[d] + '\n MIN: ' + mins[d] + ' pixels\n MAX: ' + maxs[d] + ' pixels\n'
    
        ax.text(0.5, -0.1, caption, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
    if usr_flg == 0:
        top = obj[a] + ", Min Intensity (pixels), Max Intensity (pixels), Stretch"
        info = list(zip(new_filenames_order, mins, maxs, stretch))
        f_out = root + '_img/' + new_names[a] + '_imgs/' + root + '_' + str(ra[a]) + '_' + str(dec[a]) + '-min-max.csv'
        np.savetxt(f_out, info, fmt="%s", delimiter = ", ", header=top)
        
    diameter = 2 * circle
    sexag = pyasl.coordsDegToSexa(ra[a], dec[a])
    
    if blue_vel[a] == 'Unavailable':
        source = 'Name: ' + obj[a] + '\n Average Luminosity: ' + str('%.3f' % avg_lum[a]) + ' $L_{\odot}$ \n Lower Limit Luminosity: ' + str('%.3f' % low_lum[a]) + ' $L_{\odot}$ \n Upper Limit Luminosity: ' + str('%.3f' % upp_lum[a]) + ' $L_{\odot}$ \n Average Distance: ' + str('%.3f' % avg_dist[a]) + ' kpc \n Lower Limit Distance: ' + str('%.3f' % low_dist[a]) + ' kpc \n Upper Limit Distance: ' + str('%.3f' % upp_dist[a]) + ' kpc \n Lower Velocity (Blue-Shifted Velocity): ' + blue_vel[a] + ' \n Higher Velocity (Red-Shifted Velocity): ' + red_vel[a] + ' '
    else:
        source = 'Name: ' + obj[a] + '\n Average Luminosity: ' + str('%.3f' % avg_lum[a]) + ' $L_{\odot}$ \n Lower Limit Luminosity: ' + str('%.3f' % low_lum[a]) + ' $L_{\odot}$ \n Upper Limit Luminosity: ' + str('%.3f' % upp_lum[a]) + ' $L_{\odot}$ \n Average Distance: ' + str('%.3f' % avg_dist[a]) + ' kpc \n Lower Limit Distance: ' + str(low_dist[a]) + ' kpc \n Upper Limit Distance: ' + str('%.3f' % upp_dist[a]) + ' kpc \n Lower Velocity (Blue-Shifted Velocity): ' + str('%.3f' % blue_vel[a]) + ' km/s \n Higher Velocity (Red-Shifted Velocity): ' + str('%.3f' % red_vel[a]) + ' km/s '
    
    mosaic='Size of image: ' + str(length) + ' x ' + str(width) + ' arcsec\n Diameter of circle: ' + str(diameter) + ' arcsec\n Coordinates of center of circle: ' + str(sexag)
    
    txt = source + '\n \n' + mosaic
    
    plt.figtext(0.5, 0.9, txt, wrap=True, horizontalalignment='center', fontsize=12)
    
    plt.savefig(root + '_img/' + new_names[a] + '_imgs/' + new_names[a] + '_' + str(ra[a]) + '_' + str(dec[a]) + '_mosaic.png', dpi = 150)
    print('Done creating mosaic for ' + new_names[a])
