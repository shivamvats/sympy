print "import"
from timeit import default_timer as clock
from sympy.polys.domains import ZZ, QQ
from sympy.polys.rings import ring
from sympy.polys.ring_series import rs_exp, rs_log
print "done"

def rs_mul(p1, p2, x, prec):
    ring = p1.ring
    p = {}
    if ring.__class__ != p2.ring.__class__ or ring != p2.ring:
        raise ValueError('p1 and p2 must have the same ring')
    iv = ring.gens.index(x)
    if ring == p2.ring:
        get = p.get
        items1 = list(p1.items())
        items1.sort(key=lambda e: e[0][iv])
        items2 = list(p2.items())
        items2.sort(key=lambda e: e[0][iv])
        items1 = [(a[0][0], a[1]) for a in items1]
        items2 = [(a[0][0], a[1]) for a in items2]
        assert ring.ngens == 1
        t1 = clock()
        for exp1, v1 in items1:
            for exp2, v2 in items2:
                exp = exp1 + exp2
                if exp < prec:
                    p[exp] = get(exp, 0) + v1*v2
                else:
                    break
        t2 = clock()
        print t2-t1

R, x = ring('x', QQ)

n = 200
a = rs_log(1+x, x, n)
b = rs_exp(x, x, n)
print "1"
t1 = clock()
s = rs_mul(a, b, x, n)
t2 = clock()
print "2"
#print s
print t2-t1
