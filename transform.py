import numpy as np
from ellipsoid import Ellipsoid
from transverse_mercator import TransverseMercator


class GridToGeodetic:
    def __init__(self, x, y, ellipsoid, transverse_mercator):
        self.x = x
        self.y = y

        e = Ellipsoid(ellipsoid)
        e.get_ellipsoid()
        self.e2, self.n, self.a_hat = e.ellipsoid_variables()

        tm = TransverseMercator(transverse_mercator)
        self.lam0, self.k0, self.FN, self.FE = tm.transverse_mercator_params()

    def grid_to_geodetic(self):
        epsilon = (self.x - self.FN) / (self.k0 * self.a_hat)
        new = (self.y - self.FE) / (self.k0 * self.a_hat)

        delta1 = self.n / 2.0 - 2.0 * self.n * self.n / 3.0 + 37.0 * self.n * self.n * self.n / 96.0 - self.n * self.n * self.n * self.n / 360.0
        delta2 = self.n * self.n / 48.0 + self.n * self.n * self.n / 15.0 - 437.0 * self.n * self.n * self.n * self.n / 1440.0
        delta3 = 17.0 * self.n * self.n * self.n / 480.0 - 37.0 * self.n * self.n * self.n * self.n / 840.0
        delta4 = 4397.0 * self.n * self.n * self.n * self.n / 161280.0

        epsilon_prime = (epsilon - delta1 * np.sin(2 * epsilon) * np.cosh(2 * new) -
                         delta2 * np.sin(4 * epsilon) * np.cosh(4 * new) -
                         delta3 * np.sin(6 * epsilon) * np.cosh(6 * new) -
                         delta4 * np.sin(8 * epsilon) * np.cosh(8 * new))

        new_prime = (new - delta1 * np.cos(2 * epsilon) * np.sinh(2 * new) -
                     delta2 * np.cos(4 * epsilon) * np.sinh(4 * new) -
                     delta3 * np.cos(6 * epsilon) * np.sinh(6 * new) -
                     delta4 * np.cos(8 * epsilon) * np.sinh(8 * new))

        """
        The conformal latitude φ* and the difference in longitude dλ are obtained by the formulas
        """
        psi_star = np.arcsin(np.sin(epsilon_prime) / np.cosh(new_prime))

        delta_lam = np.arctan(np.sinh(new_prime) / np.cos(epsilon_prime))

        """
        Finally, the latitude φ and the longtidue λ are obtained by the formulas
        """
        e = np.sqrt(self.e2)
        A_star = (e ** 2 + e ** 4 + e ** 6 + e ** 8)
        B_star = -(1 / 6) * (7 * e ** 4 + 17 * e ** 6 + 30 * e ** 8)
        C_star = (1 / 120) * (224 * e ** 6 + 889 * e ** 8)
        D_star = -(1 / 1260) * (4279 * e ** 8)

        lam = self.lam0 + delta_lam
        psi = (psi_star + np.sin(psi_star) * np.cos(psi_star) *
               (A_star +
                B_star * (np.sin(psi_star)) ** 2 +
                C_star * (np.sin(psi_star)) ** 4 +
                D_star * (np.sin(psi_star)) ** 6))

        latlon = (lam * 180 / np.pi, psi * 180 / np.pi)

        return latlon


class GeodeticToGeocentric:
    def __init__(self, psi, lam, h, ellipsoid):
        deg_to_rad = np.pi / 180.0

        self.psi = psi * deg_to_rad
        self.lam = lam * deg_to_rad
        self.h = h

        e = Ellipsoid(ellipsoid)
        self.a, self.f = e.get_ellipsoid()

    def geodetic_to_geocentric(self):
        """
        X,Y,Z           geocentric cartesian coordinates
        psi,lambda,h    geodetic coordinates (h = height above ellipsoid)
        N_prime         tvärkrökningsradien
        e               first eccentricity
        f               flattening
        a               half major-axis
        """
        e2 = self.f * (2.0 - self.f)
        N_prime = self.a / np.sqrt(1.0 - e2 * (np.sin(self.psi))**2.0)

        X = (N_prime + self.h) * np.cos(self.psi) * np.cos(self.lam)
        Y = (N_prime + self.h) * np.cos(self.psi) * np.sin(self.lam)
        Z = (N_prime * (1 - e2) + self.h) * np.sin(self.psi)

        XYZ = (X, Y, Z)

        return XYZ
