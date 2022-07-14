class Ellipsoid:
    def __init__(self, ellipsoid):
        self.ellipsoid = ellipsoid

    def ellipsoid_variables(self):
        """
        The following variables are defined out of the ellipsoidal parameters a and f
        """

        e2 = self.f * (2.0 - self.f)
        n = self.f / (2.0 - self.f)
        a_hat = self.a / (1.0 + n) * (1.0 + n * n / 4.0 + n * n * n * n / 64.0)

        return e2, n, a_hat

    def get_ellipsoid(self):
        """
        Get ellipsoid parameters a and f

            a   semi-major axis [m]
            f   flattening
        """

        if self.ellipsoid == 'grs1980':
            self.a = 6378137.0
            self.f = 1.0 / 298.257222101

        return self.a, self.f
