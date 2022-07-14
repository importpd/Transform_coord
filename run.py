from transform import GridToGeodetic, GeodeticToGeocentric

ellipsoid = 'grs1980'
utm = 'test'

########################################################

x = 1135809.413803
y = 555304.016555

t = GridToGeodetic(x, y, ellipsoid, utm)
latlon = t.grid_to_geodetic()
#print(latlon)

########################################################

psi = 58.0
lam = 17.0
h = 30.0

t = GeodeticToGeocentric(psi, lam, h, ellipsoid)
XYZ = t.geodetic_to_geocentric()
print(XYZ)