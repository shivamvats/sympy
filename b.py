print "import"
from timeit import default_timer as clock
from sympy.polys.domains import ZZ, QQ
from sympy.polys.rings import ring
from sympy.polys.ring_series import rs_mul, rs_exp, rs_log
print "done"

R, x = ring('x', QQ)

n = 400
print "initialize series"
a = rs_log(1+x, x, n)
b = rs_exp(x, x, n)
print "done"
print "start"
t1 = clock()
s = rs_mul(a, b, x, n)
t2 = clock()
print "stop"
#print s
print t2-t1
