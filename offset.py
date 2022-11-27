#This program calculates the offset between two user-specified coordinates in arcseconds
#To run this program: python offset.py RA1 Dec1 RA2 Dec2 (make sure the coordinates are in degrees)
import astropy.units as u
from astropy.coordinates import SkyCoord
import sys

c1 = SkyCoord(sys.argv[1], sys.argv[2], unit="deg", frame = 'icrs')
c2 = SkyCoord(sys.argv[3], sys.argv[4], unit="deg", frame = 'icrs')

separation=c1.separation(c2).to(u.arcsec)
print(separation)