print "import"
from timeit import default_timer as clock
from sympy import S
print "done"

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

    def __str__(self):
        s = ""
        for i in range(len(self._data)):
            if self._data[i] != 0:
                s += "(%s)*X^%d + " % (self._data[i], i)
        return s[:-3]

class FormalPowerSeriesSparse(object):

    def __init__(self, data):
        # Only non-zero elements are stored
        # _data = {0: a0, 1: a1, 2: a2, ...}
        self._data = {}
        for n, x in enumerate(data):
            y = S(x)
            if y != 0: self._data[n] = y

    def __str__(self):
        items = list(self._data.items())
        items.sort(key=lambda e: e[0])
        s = ""
        for e, c in items:
            s += "(%s)*X^%d + " % (c, e)
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

def mul_sparse(p1, p2, prec):
    """
    Multiplies two series (a*b). The order is determined by the shorter one.
    """
    items2 = list(p2._data.items())
    items2.sort(key=lambda e: e[0])
    p = {}
    get = p.get
    for exp1, v1 in p1._data.items():
        for exp2, v2 in items2:
            exp = exp1 + exp2
            if exp < prec:
                p[exp] = get(exp, 0) + v1*v2
            else:
                break
    return p

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

n = 200
# sin(x)
data = [0]
t = S(1)
print "0"
for i in range(1, n):
    t = t/i
    if i % 2 == 0:
        data.append(0)
        t = -t
    else:
        data.append(t)
sin = FormalPowerSeries(data)
# cos(x)
data = [1]
t = S(1)
for i in range(1, n):
    t = t/i
    if i % 2 == 1:
        data.append(0)
        t = -t
    else:
        data.append(t)
cos = FormalPowerSeries(data)
# exp(x)
data = [1]
t = S(1)
for i in range(1, n):
    t = t/i
    data.append(t)
exp = FormalPowerSeriesSparse(data)
# log(1+x)
data = [0]
t = S(1)
for i in range(1, n):
    data.append(t/i)
    t = -t
log = FormalPowerSeriesSparse(data)
#x = FormalPowerSeries([0, 1, 0, 0, 0, 0, 0, 0, 0])
#onemx = FormalPowerSeries([1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
print "4"
t1 = clock()
s = mul_sparse(log, exp, n)
t2 = clock()
#print s
print t2-t1
