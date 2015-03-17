from timeit import default_timer as clock
from sympy import QQ
from sympy.polys.ring_series import rs_mul, rs_square, rs_series_from_list, rs_pow
from sympy.polys.rings import sring, PolyElement


class FormalPowerSeries(object):

    """
    Formal power series 'alpha':

    alpha = [a0, a1, a2, ....] = a0 + a1*X + a2*X^2 + ...

    It is called formal, because X is just a formal symbol, you cannot
    substitute for it. All questions of convergence are absent. The series is
    fully specified by the numbers [a0, a1, a2, ...].

    The series always starts at i=0. Formal Laurent series then allows to start
    from a negative exponent i = -n.

    The formal power series is infinite. But we only store finite number of
    terms. As such, the order term is determined as follows:

    _data = [a0, a1, a2, ...., an]
        = a0 + a1*X + a2*X^2 + ... + an*X^n + O(X^{n+1})

    The representation is sparse.

    References:

    [1] http://en.wikipedia.org/wiki/Formal_power_series

    [2] https://www.maa.org/sites/default/files/pdf/upload_library/22/Ford/IvanNiven.pdf

    [3] http://planetmath.org/invertibleformalpowerseries

    [4] http://www.math.ucsd.edu/~jverstra/264A-LECTUREB.pdf
    """

    def __init__(self, data):
        if isinstance(data, PolyElement):
            self.series = data
            self.R = data.ring
        else:
             self.R, self.series = sring(data, domain=QQ)

    def __repr__(self):
        from sympy.printing import sstr
        return sstr(self.series, order=None)

    def __str__(self):
        from sympy.printing import sstr
        return sstr(self.series, order=None)

    def __add__(self, other):
        return FormalPowerSeries(self.series + other.series)

    def __mul__(self, other):
        if(self.R.ngens == 1 and other.R.ngens == 1):
            x = ((self.R).gens)[0]
            prec = (self.series).degree() + (other.series).degree() + 1
            return FormalPowerSeries(rs_mul(self.series, other.series, x, prec))

    def __pow__(self, n):
        if(self.R.ngens == 1):
            x = ((self.R).gens)[0]
            prec = (self.series).degree() + 1
            return FormalPowerSeries(rs_pow(self.series, n, x, prec))

    def __div__(self, other):
        if(self.R.ngens == 1 and other.R.ngens == 1):
            x = ((self.R).gens)[0]
            return self * (other)**1
