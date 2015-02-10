from sympy import S

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

    The representation is dense.

    References:

    [1] http://en.wikipedia.org/wiki/Formal_power_series

    [2] https://www.maa.org/sites/default/files/pdf/upload_library/22/Ford/IvanNiven.pdf

    [3] http://planetmath.org/invertibleformalpowerseries

    [4] http://www.math.ucsd.edu/~jverstra/264A-LECTUREB.pdf
    """

    def __init__(self, data):
        # _data = [a0, a1, a2, ...]
        self._data = [S(x) for x in data]
        pass

    def __str__(self):
        s = ""
        for i in range(len(self._data)):
            if self._data[i] != 0:
                s += "(%s)*X^%d + " % (self._data[i], i)
        return s[:-3]

def add(a, b):
    """
    Sums two series (a+b). The order is determined by the shorter one.
    """
    an = len(a._data)
    bn = len(b._data)
    cn = min(an, bn)
    cdata = [0]*cn
    for n in range(cn):
        cdata[n] = a._data[n]+b._data[n]
    return FormalPowerSeries(cdata)

def mul(a, b):
    """
    Multiplies two series (a*b). The order is determined by the shorter one.
    """
    an = len(a._data)
    bn = len(b._data)
    cn = min(an, bn)
    cdata = [0]*cn
    for n in range(cn):
        s = 0
        for j in range(n+1):
            s = s + a._data[j]*b._data[n-j]
        cdata[n] = s
    return FormalPowerSeries(cdata)

def pow_m1(a):
    """
    Calculates 1/a.  The order is determined by 'a'.
    """
    if a._data[0] == 0:
        raise ValueError("It must hold: a0 != 0")
    an = len(a._data)
    bdata = [0]*an
    bdata[0] = 1/a._data[0]
    for n in range(1, an):
        s = 0
        for i in range(1, n+1):
            s = s + a._data[i]*bdata[n-i]
        bdata[n] = - s / a._data[0]
    return FormalPowerSeries(bdata)

def div(a, b):
    """
    Calculates a/b.
    """
    return mul(a, pow_m1(b))

# sin(x)
sin = FormalPowerSeries([0, 1, 0, -S(1)/6, 0, S(1)/120, 0, -S(1)/5040, 0,
    S(1)/362880, 0])
# cos(x)
cos = FormalPowerSeries([1, 0, -S(1)/2, 0, S(1)/24, 0, -S(1)/720, 0,
    S(1)/40320, 0])
x = FormalPowerSeries([0, 1, 0, 0, 0, 0, 0, 0, 0])
onemx = FormalPowerSeries([1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
print sin
print cos
print mul(sin, cos)
print add(sin, cos)
print mul(sin, x)
print div(sin, cos)
