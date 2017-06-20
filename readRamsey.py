import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import numpy as np
import math
from scipy import optimize

numPhases = 12

delBlist = np.array([45,90,135])

delB = float (delBlist[1])
print delB

#Read in data
data = np.loadtxt('RamseyEcho_%dPulse.txt' % (delB))
file = open('RamseyResults_%dPulse.txt' % (delB), 'w')

#Unpack data
for i in range(len(data)/numPhases):
    startval = (i)*numPhases
    endval = (i+1)*numPhases
    block = data[startval:endval,:]
    phase = block[:,1]
    frac = -block[:,2]
    time = block[0,0]

    #Fit sinusoid
    p = np.array([0,0.7])
    def peval(phase,p):
        return p[1]*np.cos((p[0]+phase)*math.pi/180)
    def residuals(p,phase,frac):
        return frac - peval(phase,p)
    pf, C, info, message, success = optimize.leastsq(residuals,p,args=(phase,frac),full_output=1)
    print pf
    file.write('%f \t %f \t %f \n' % (time,pf[0],math.fabs(pf[1])))


    #Plot
    plt.plot(phase,frac,'o')
    plt.plot(phase,peval(phase,pf),'-')
    plt.hold(True)

plt.xlim(-190,190)
plt.ylim(-1.0,1.0)
plt.xlabel('Phase')
plt.ylabel('Fraction')
plt.show()
file.close()
