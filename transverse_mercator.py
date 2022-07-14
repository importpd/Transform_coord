import numpy as np


class TransverseMercator:
    def __init__(self, utm):
        self.utm = utm

    def ddmmss_to_decimal(self, ddmmss):
        return ddmmss[0] + ddmmss[1] / 60.0 + ddmmss[2] / 3600.0

    def transverse_mercator_params(self):
        """
        Transverse mercator parameters
            central_meridian    central meridian
            lam0                longitude of the central meridian
            k0                  scale factor of the central meridian
            FN                  false northing [m]
            FE                  false easting [m]
        """

        deg_to_rad = np.pi / 180.0

        if self.utm == 'test':
            central_meridian = [13, 35, 7.692000]
            lam0 = self.ddmmss_to_decimal(central_meridian) * deg_to_rad
            k0 = 1.000002540000
            FN = -6226307.8640
            FE = 84182.8790

        return lam0, k0, FN, FE
