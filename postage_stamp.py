#Program Purpose: This program acquires images from various missions specified in the root-img.inp file from IRSA. The coordinates used to search the IRSA finder chart are from GAIA EDR3.
#To run this program, type: python postage_stamp.py root-img.inp

import sys
import requests
import os
import shutil
import zipfile
import io
from astropy.coordinates import SkyCoord
from astropy import units as u
from pathlib import Path
from subprocess import call

root = sys.argv[1][:-8]
coordfile = root + "-prgb-final-candidates.csv"
inpfile = root + "-img.inp"

#reading in root_img.inp parameters
with open(inpfile, 'r') as inp:
    lines = inp.readlines()
    for line in lines:
        if line[:1] == '#':
            continue
        else:
            index = line.find('=')
            if 'survey' in line[:index]:
                survey = str(line[index+1:].strip())
            elif 'subsetsize' in line[:index]:
                subsetsize = float(line[index+1:].strip())
            elif 'c_format' in line[:index]:
                c_format = str(line[index+1:].strip())
            elif 'sep' in line[:index]:
                sep = str(line[index+1:].strip())
    
ra = []
dec = []
obj = []
new_names = []

#creating root_img directory to store images for each star
path = str(Path(__file__).parent.absolute())
dir = path + '/' + str(root) + '_img'
os.makedirs(dir, exist_ok=True)

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
#             elif c_format =='b':
#                 if x[3].strip() =='d':
#                     ra.append(float(x[1]))
#                     dec.append(float(x[2]))
#                 elif x[3].strip()=='h':
#                     c = SkyCoord(str(x[1])+' '+str(x[2]), unit=(u.hourangle, u.deg))
#                     ra.append(float(c.ra.degree))
#                     dec.append(float(c.dec.degree))


#creating image directories for every star
for b in range(len(obj)):
    new_name = obj[b].replace(" ", "_")
    new_names.append(new_name)
    os.makedirs(dir + '/' + new_name + '_imgs', exist_ok=True)
                    
                    
#UPDATE (11/14/21): not sure if this code is needed anymore. function to be used in the code block below. takes all images and moves them into the parent directory and removes any empty folders. this code was taken from: https://stackoverflow.com/questions/8428954/move-child-folder-contents-to-parent-folder-in-python (by George)                    
def move_to_root_folder(root_path, cur_path):
    for filename in os.listdir(cur_path):
        if os.path.isfile(os.path.join(cur_path, filename)):
            shutil.move(os.path.join(cur_path, filename), os.path.join(root_path, filename))
        elif os.path.isdir(os.path.join(cur_path, filename)):
            move_to_root_folder(root_path, os.path.join(cur_path, filename))
        else:
            sys.exit("Should never reach here.")
    # remove empty folders
    if cur_path != root_path:
        os.rmdir(cur_path)
                    
#for each star, we search for cutouts on the IRSA website by entering a link containing our parameters and retrieving a ZIP file, which contains the fit files for each star. Unzipping code taken from: https://stackoverflow.com/questions/9419162/download-returned-zip-file-from-url                    
for n in range(len(ra)):
    link_image = "https://irsa.ipac.caltech.edu/applications/finderchart/servlet/api?mode=getImage&locstr=" + str(ra[n]) + "+" + str(dec[n]) + "&survey=" + survey + "&subsetsize=" + str(subsetsize) + "&marker=true;reproject=true"
    
    print(link_image)
    req = requests.get(link_image,allow_redirects=True, verify=False)
    z = zipfile.ZipFile(io.BytesIO(req.content))
    z.extractall(dir + '/' + new_names[n] + '_imgs')
    

    

    