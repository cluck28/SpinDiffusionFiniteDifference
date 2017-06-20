import numpy
import math
import time

def writeData(outputFile,array,t,gridSize,spinDimension):
    for i in range(spinDimension):
        outputFile.write('%f \t' % (t))
        for j in range(gridSize):
            outputFile.write(str(array[i][j])+' ')
        outputFile.write('\n')

def writeMagData(outputFile,array,t,gridSize):
    outputFile.write('%f ' % (t))
    for i in range(gridSize):
        outputFile.write(str(array[i])+' ')
    outputFile.write('\n')

def writeMathOutput(outputFile,gridAvgMag,errMag,t):
    outputFile.write('%f \t %f \t %f\n' % (t,gridAvgMag,errMag))

def writeRamsey(outputFile,time,phase,fraction):
    outputFile.write('%f \t %f \t %f\n' % (time,phase,fraction))

#Start Time
print time.strftime("%H:%M:%S")
#Initialize variables
spinDimension = 3
gridLength1 = 5
gridLength = gridLength1*0.2*10**-5
gridSize = 300+1
deltaZ = gridLength/gridSize
print deltaZ
deltaT = (10**-7)/1.0
ZStep = gridLength/gridSize
#Make times list
finalT1 = 1.5/2.0 #1.5
initialT1 = 0.06
timeStep1 = 0.06 #0.06
lengthtimelist1 = int((finalT1-initialT1)/timeStep1)
timelist = numpy.zeros(lengthtimelist1)
for i in range(lengthtimelist1):
    timelist[i] = (i*timeStep1+initialT1)*10**-3
finalT2 = 3.0/2.0 #3.0
initialT2 = finalT1
timeStep2 = 0.15 #0.15
lengthtimelist2 = int((finalT2-initialT2)/timeStep2)
timelist2 = numpy.zeros(lengthtimelist2)
for i in range(lengthtimelist2):
    timelist2[i] = (i*timeStep2+initialT2)*10**-3
finalT3 = 5.5/2.0 #5.5
initialT3 = finalT2
timeStep3 = 0.35 #0.35
lengthtimelist3 = int((finalT3-initialT3)/timeStep3)
timelist3 = numpy.zeros(lengthtimelist3)
for i in range(lengthtimelist3):
    timelist3[i] = (i*timeStep3+initialT3)*10**-3
timelist = numpy.concatenate((timelist,timelist2))
timelist = numpy.concatenate((timelist,timelist3))

#Relevant constants
uB = 9.274*10**-24
h = 6.626*10**-34
hbar = h/(2*math.pi)
grad = 20 
delB = grad*10**-2
gamma = 0.112*uB/hbar
u = 1.66*10**-27
m = 39.0983*u
D=2.0*hbar/m
pitime1 = 120
pitime = pitime1*10**-6
print timelist*2+pitime
rabi =  1./pitime/2.
rabi = 2*math.pi*rabi
tipAngle = math.pi/2.0 
tipAngleDeg = tipAngle*180/math.pi
M0 = math.sin(tipAngle)
G = delB*gamma
uu = -1.0
posScale = (G/D)**(1./3.)
timeScale = ((G**2.)*D)**(1./3.)
timelist = timelist*timeScale
deltaT = deltaT*timeScale
deltaZ = deltaZ*posScale
print (deltaT)/(deltaZ*deltaZ)
tAdjust = 2
deltaT1 = deltaT/tAdjust
deltaZ1 = deltaZ/2
uu = uu*M0
pitime = pitime*timeScale
rabi = rabi/timeScale

#Open file to save data
outputFile = open('SpinComponentsEcho_%dPulse.txt' % (tipAngleDeg), 'w')
outputFile2 = open('MagnitudeMagnetizationEcho_%dPulse.txt' % (tipAngleDeg),'w')
outputFile3 = open('ZProjectionEcho_%dPulse.txt' % (tipAngleDeg),'w')
outputFile4 = open('TrapAverageMag_%dPulse.txt' % (tipAngleDeg), 'w')
outputFile5 = open('XProjectionEcho_%dPulse.txt' % (tipAngleDeg), 'w')
outputFile6 = open('YProjectionEcho_%dPulse.txt' % (tipAngleDeg), 'w')
outputFile7 = open('RamseyEcho_%dPulse.txt' % (tipAngleDeg), 'w')

#Write header
outputFile.write('#Spin dimension '+str(spinDimension)+', grid size '+str(gridSize)+', delta t '+str(deltaT)+'\n')
outputFile2.write('#Spin dimension '+str(spinDimension)+', grid size '+str(gridSize)+', delta t '+str(deltaT)+'\n')
outputFile3.write('#Spin dimension '+str(spinDimension)+', grid size '+str(gridSize)+', delta t '+str(deltaT)+'\n')
outputFile4.write('#time \t Mag \t e_Mag\n')
outputFile5.write('#Spin dimension '+str(spinDimension)+', grid size '+str(gridSize)+', delta t '+str(deltaT)+'\n')
outputFile6.write('#Spin dimension '+str(spinDimension)+', grid size '+str(gridSize)+', delta t '+str(deltaT)+'\n')
outputFile7.write('#time \t phase \t fraction\n')

#Make spatially dependent lamour freq
omega = numpy.zeros(gridSize)
omegagen = numpy.zeros(gridSize)
detuning = numpy.zeros(gridSize)
#Rotations during pi pulse
Rot11 = numpy.zeros(gridSize)
Rot12 = numpy.zeros(gridSize)
Rot13 = numpy.zeros(gridSize)
Rot21 = numpy.zeros(gridSize)
Rot22 = numpy.zeros(gridSize)
Rot23 = numpy.zeros(gridSize)
Rot31 = numpy.zeros(gridSize)
Rot32 = numpy.zeros(gridSize)
Rot33 = numpy.zeros(gridSize)
for i in range(gridSize):
    scaledpos = i-gridSize/2
    omega[i] = posScale*(scaledpos*ZStep)
    omegagen[i] = -numpy.sqrt(omega[i]*omega[i]+rabi*rabi) #why minus sign
    detuning[i] = numpy.arctan(omega[i]/rabi)
    Rot11[i] = numpy.cos(detuning[i])*numpy.cos(detuning[i])+numpy.sin(detuning[i])*numpy.cos(omegagen[i]*deltaT1)*numpy.sin(detuning[i])
    Rot12[i] = numpy.sin(detuning[i])*numpy.sin(omegagen[i]*deltaT1)
    Rot13[i] = -numpy.cos(detuning[i])*numpy.sin(detuning[i])+numpy.cos(detuning[i])*numpy.cos(omegagen[i]*deltaT1)*numpy.sin(detuning[i])
    
    Rot21[i] = -numpy.sin(detuning[i])*numpy.sin(omegagen[i]*deltaT1)
    Rot22[i] = numpy.cos(omegagen[i]*deltaT1)
    Rot23[i] = -numpy.cos(detuning[i])*numpy.sin(omegagen[i]*deltaT1)
    
    Rot31[i] = -numpy.cos(detuning[i])*numpy.sin(detuning[i])+numpy.cos(detuning[i])*numpy.cos(omegagen[i]*deltaT1)*numpy.sin(detuning[i])
    Rot32[i] = numpy.cos(detuning[i])*numpy.sin(omegagen[i]*deltaT1)
    Rot33[i] = numpy.cos(detuning[i])*numpy.cos(detuning[i])*numpy.cos(omegagen[i]*deltaT1)+numpy.sin(detuning[i])*numpy.sin(detuning[i])
print omega[gridSize-1]/posScale

phaselist = numpy.arange(-math.pi,math.pi,math.pi/6)
#Make rotation for each phase for Ramsey
A = numpy.zeros(len(phaselist))
B = numpy.zeros(len(phaselist))
C = numpy.zeros(len(phaselist))
Dnew = numpy.zeros(len(phaselist))
E = numpy.zeros(len(phaselist))
F = numpy.zeros(len(phaselist))
Gnew = numpy.zeros(len(phaselist))
H = numpy.zeros(len(phaselist))
I = numpy.zeros(len(phaselist))
for i in range(len(phaselist)):
    A[i] = numpy.cos(phaselist[i])*numpy.cos(phaselist[i])
    B[i] = -numpy.cos(phaselist[i])*numpy.sin(phaselist[i])
    C[i] = -numpy.sin(phaselist[i])
    Dnew[i] = -numpy.cos(phaselist[i])*numpy.sin(phaselist[i])
    E[i] = numpy.sin(phaselist[i])*numpy.sin(phaselist[i])
    F[i] = -numpy.cos(phaselist[i])
    Gnew[i] = numpy.sin(phaselist[i])
    H[i] = numpy.cos(phaselist[i])
    I[i] = 0

#Make grid
array = numpy.zeros((spinDimension,gridSize))
arrayOld = numpy.zeros((spinDimension,gridSize))
arrayMag = numpy.zeros((gridSize))
zProj = numpy.zeros((gridSize))
xProj = numpy.zeros((gridSize))
yProj = numpy.zeros((gridSize))
d2m = numpy.zeros((spinDimension))
m2 = 0

#Initial conditions
arrayOld[0,:] = 0
arrayOld[1,:] = math.sin(tipAngle)/M0
arrayOld[2,:] = math.cos(tipAngle)/M0

#An array that acts as a placeholder before the pi pulse
placeHoldArray = numpy.zeros((spinDimension,gridSize))
placeHoldArray = numpy.copy(arrayOld)
placeHoldTime = 0.

#An array that is a placeholder before the Ramsey
ramseyArray = numpy.zeros((spinDimension,gridSize))

#Loop for different pi times
for tpi in timelist:
    t0 = time.time()

    #Get starting point
    t = placeHoldTime
    arrayOld = numpy.copy(placeHoldArray)

    #Carry out numerical integration
    #Precess
    while t<tpi:
    #End points
        for i in range(spinDimension):
            array[i,0] = 0
            array[i,gridSize-1] = 0
        #Calculate relevant arrays
        for i in range(1,gridSize-1):
            #Calculate M^2
            m2 = arrayOld[0,i]*arrayOld[0,i]+arrayOld[1,i]*arrayOld[1,i]+arrayOld[2,i]*arrayOld[2,i]
            #Calculate second derivatives
            d2m[0] = (arrayOld[0,i-1]-2.0*arrayOld[0,i]+arrayOld[0,i+1])/(deltaZ*deltaZ)
            d2m[1] = (arrayOld[1,i-1]-2.0*arrayOld[1,i]+arrayOld[1,i+1])/(deltaZ*deltaZ)
            d2m[2] = (arrayOld[2,i-1]-2.0*arrayOld[2,i]+arrayOld[2,i+1])/(deltaZ*deltaZ)
            array[0,i] = arrayOld[0,i]-deltaT*arrayOld[1,i]*omega[i]+(deltaT/(1+uu*uu*m2))*(d2m[0]+uu*(arrayOld[1,i]*d2m[2]-arrayOld[2,i]*d2m[1]))
            array[1,i] = arrayOld[1,i]+deltaT*arrayOld[0,i]*omega[i]+(deltaT/(1+uu*uu*m2))*(d2m[1]+uu*(arrayOld[2,i]*d2m[0]-arrayOld[0,i]*d2m[2]))
            array[2,i] = arrayOld[2,i]+(deltaT/(1+uu*uu*m2))*(d2m[2]+uu*(arrayOld[0,i]*d2m[1]-arrayOld[1,i]*d2m[0]))
        arrayOld = numpy.copy(array)
        t = t+deltaT

    #Save placeholder
    placeHoldArray = numpy.copy(array)
    placeHoldTime = t

    #Echo pi-x
    while t<tpi+2*pitime:
        #Endpoints
        for i in range(spinDimension):
            array[i,0] = 0
            array[i,gridSize-1]=0
        for i in range(1,gridSize-1):
            #Do rotation
            array[0,i] = Rot11[i]*arrayOld[0,i]+Rot12[i]*arrayOld[1,i]+Rot13[i]*arrayOld[2,i]
            array[1,i] = Rot21[i]*arrayOld[0,i]+Rot22[i]*arrayOld[1,i]+Rot23[i]*arrayOld[2,i]
            array[2,i] = Rot31[i]*arrayOld[0,i]+Rot32[i]*arrayOld[1,i]+Rot33[i]*arrayOld[2,i]
        arrayOld = numpy.copy(array)
        for i in range(1,gridSize-1):
            #Calculate M^2
            m2 = arrayOld[0,i]*arrayOld[0,i]+arrayOld[1,i]*arrayOld[1,i]+arrayOld[2,i]*arrayOld[2,i]
            #Calculate second derivatives
            d2m[0] = (arrayOld[0,i-1]-2.0*arrayOld[0,i]+arrayOld[0,i+1])/(deltaZ*deltaZ)
            d2m[1] = (arrayOld[1,i-1]-2.0*arrayOld[1,i]+arrayOld[1,i+1])/(deltaZ*deltaZ)
            d2m[2] = (arrayOld[2,i-1]-2.0*arrayOld[2,i]+arrayOld[2,i+1])/(deltaZ*deltaZ)
            array[0,i] = arrayOld[0,i]+(deltaT/(1+uu*uu*m2))*(d2m[0]+uu*(arrayOld[1,i]*d2m[2]-arrayOld[2,i]*d2m[1]))
            array[1,i] = arrayOld[1,i]+(deltaT/(1+uu*uu*m2))*(d2m[1]+uu*(arrayOld[2,i]*d2m[0]-arrayOld[0,i]*d2m[2]))
            array[2,i] = arrayOld[2,i]+(deltaT/(1+uu*uu*m2))*(d2m[2]+uu*(arrayOld[0,i]*d2m[1]-arrayOld[1,i]*d2m[0]))
        arrayOld = numpy.copy(array)
        t = t+deltaT
#    array[0,:] = -arrayOld[0,:]
#    array[1,:] = arrayOld[1,:]
#    array[2,:] = arrayOld[2,:]
#    arrayOld = numpy.copy(array)
#    t = t+pitime

    #Precess
    while t<2*tpi+2*pitime:
        #End points
        for i in range(spinDimension):
            array[i,0] = 0
            array[i,gridSize-1] = 0
        for i in range(1,gridSize-1):
            #Calculate M^2
            m2 = arrayOld[0,i]*arrayOld[0,i]+arrayOld[1,i]*arrayOld[1,i]+arrayOld[2,i]*arrayOld[2,i]
            #Calculate second derivatives
            d2m[0] = (arrayOld[0,i-1]-2.0*arrayOld[0,i]+arrayOld[0,i+1])/(deltaZ*deltaZ)
            d2m[1] = (arrayOld[1,i-1]-2.0*arrayOld[1,i]+arrayOld[1,i+1])/(deltaZ*deltaZ)
            d2m[2] = (arrayOld[2,i-1]-2.0*arrayOld[2,i]+arrayOld[2,i+1])/(deltaZ*deltaZ)
            array[0,i] = arrayOld[0,i]-deltaT*arrayOld[1,i]*omega[i]+(deltaT/(1+uu*uu*m2))*(d2m[0]+uu*(arrayOld[1,i]*d2m[2]-arrayOld[2,i]*d2m[1]))
            array[1,i] = arrayOld[1,i]+deltaT*arrayOld[0,i]*omega[i]+(deltaT/(1+uu*uu*m2))*(d2m[1]+uu*(arrayOld[2,i]*d2m[0]-arrayOld[0,i]*d2m[2]))
            array[2,i] = arrayOld[2,i]+(deltaT/(1+uu*uu*m2))*(d2m[2]+uu*(arrayOld[0,i]*d2m[1]-arrayOld[1,i]*d2m[0]))
        arrayOld = numpy.copy(array)
        t = t+deltaT

    #End of program stuff
    t = t/timeScale
    print t
    for i in range(gridSize):
        arrayMag[i] = M0*math.sqrt(array[0,i]*array[0,i]+array[1,i]*array[1,i])
        zProj[i] = M0*array[2,i]
        xProj[i] = M0*array[0,i]
        yProj[i] = M0*array[1,i]
    writeData(outputFile,array,t,gridSize,spinDimension)
    writeMagData(outputFile2,arrayMag,t,gridSize)
    writeMagData(outputFile3,zProj,t,gridSize)
    writeMagData(outputFile5,xProj,t,gridSize)
    writeMagData(outputFile6,yProj,t,gridSize)
    print time.time()-t0
    #Print out the grid averaged magnetization and standard deviation
    gridAvg = numpy.sum(arrayMag[gridSize/2-gridSize/4:gridSize/2+gridSize/4])/(gridSize/2)
    gridStd = numpy.sqrt(numpy.var(arrayMag[gridSize/2-gridSize/4:gridSize/2+gridSize/4])/((gridSize/2)-1))
    writeMathOutput(outputFile4,gridAvg,gridStd,t)

    #Do Ramsey pulse
    for i in range(len(phaselist)):
        phase = phaselist[i]*180/math.pi
        ramseyArray[0,:] = arrayOld[0,:]*A[i] + arrayOld[1,:]*B[i] + arrayOld[2,:]*C[i]
        ramseyArray[1,:] = arrayOld[0,:]*Dnew[i]+arrayOld[1,:]*E[i]+arrayOld[2,:]*F[i]
        ramseyArray[2,:] = arrayOld[0,:]*Gnew[i]+arrayOld[1,:]*H[i]+arrayOld[2,:]*I[i]
        fraction = M0*numpy.sum(ramseyArray[2,gridSize/2-gridSize/4:gridSize/2+gridSize/4])/(gridSize/2)
        writeRamsey(outputFile7,t,phase,fraction)


#Close files
outputFile.close()
outputFile2.close()
outputFile3.close()
outputFile4.close()
outputFile5.close()
outputFile6.close()
outputFile7.close()
#End Time
print time.strftime("%H:%M:%S")    
