import math
n = 90333
of = open("plot.3D", "wt")
vf = open("value.data", "r")
yf = open("yy.data", "r")
xf = open("xx.data", "r")
of.write("x\ty\tz\tvalue\n");
for i in range(0,3):
    x= xf.readline().strip()
    y= yf.readline().strip()
    v= vf.readline().strip()
for i in range(4, n):
    x= xf.readline().strip()
    y= yf.readline().strip()
    v= vf.readline().strip()
    of.write("%s\t%s\t%s\t%s\n" %(x, y, 0.0, v))
    #t = float(i) / float(n-1)
    #angle = t * (math.pi * 2.) * 50.
    #r = t * 10.
    #x = r * math.cos(angle)
    #y = r * math.sin(angle)
    #z = t * 10.
    #value = math.sqrt(x*x + y*y + z*z)
    #f.write("%g %g %g %g\n" % (x,y,z,value))
of.close()
vf.close()
xf.close()
yf.close()
