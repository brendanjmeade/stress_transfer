import math
import numpy as np
import scipy
from okada_wrapper import dc3d0wrapper, dc3dwrapper
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import cm

# load file and extract geometric coordiantes and slip distribution
fileName = 's1999HECTOR01SALI'
N = 100 # Number of grid points in x and y-directions for visualization
lambdaLame = 0.25 # First Lame parameter
muLame = 0.25 # shear modulus
coefficientOfFriction = 0.4 # Coefficient of friction
km2m = 1e3 # Convert kilometers to meters
cm2m = 1e-2 # Convert centimeters to meters
cfsUpperLimit = 5e-6; # for visualziation purposes
cfsLowerLimit = -5e-6; # for visualization purposes
obsDepth = -5e3; # depth of observation coordinates

# Load the .mat (HDF5-ish) version of the model geometry and slip distribution
F = scipy.io.loadmat(fileName + '.mat')
F = F[fileName]
F = F[0]
S = list() # Empty list

# Extract the fault geometry and slip into a single stucture of ungrouped patches
nPanel = int(F['invSEGM'][0][0][0]) # Count number of panels
for iPanel in range(1, nPanel + 1): # This index starts at 1 because of naming convention of the fields in the dictionary from the .mat file
    # Extract geometric parameters from this panel common to all patches
    strike = 0
    strike = F['seg' + str(iPanel) + 'AStke'][0][0][0]
    angle = -(strike-90)
    if angle < 0:
        angle = angle + 360

    # calculate the length and wide if individual patch elements in current panel
    L = F['seg' + str(iPanel) + 'DimWL'][0][0][1] / np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[1]
    W = F['seg' + str(iPanel) + 'DimWL'][0][0][1] / np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[0]
    rTemp = np.array([[math.cos(math.radians(angle)), -math.sin(math.radians(angle))], [math.sin(math.radians(angle)), math.cos(math.radians(angle))]])
    xTempOrig = np.array([[L/2.0], [0.0]])
    xTempRot = np.dot(rTemp, xTempOrig)
    xTopOffset = xTempRot[0];
    yTopOffset = xTempRot[1];
    zTopOffset = 0;
    xTopBottomOffset = F['seg' + str(iPanel) + 'geoX'][0][1][0] - F['seg' + str(iPanel) + 'geoX'][0][0][0]
    yTopBottomOffset = F['seg' + str(iPanel) + 'geoY'][0][1][0] - F['seg' + str(iPanel) + 'geoY'][0][0][0]
    zTopBottomOffset = F['seg' + str(iPanel) + 'geoZ'][0][1][0] - F['seg' + str(iPanel) + 'geoZ'][0][0][0]
    nDownDip = np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[0]
    nAlongStrike = np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[1]

    # Loops over the down-dip and along-strike patches of the current panel
    for iDownDip in range(0, nDownDip):
        for iAlongStrike in range(0, nAlongStrike):
            
            # Extract top center coordinates of current patch
            xTopCenter = F['seg' + str(iPanel) + 'geoX'][0][iDownDip][iAlongStrike]
            yTopCenter = F['seg' + str(iPanel) + 'geoY'][0][iDownDip][iAlongStrike]
            zTopCenter = F['seg' + str(iPanel) + 'geoZ'][0][iDownDip][iAlongStrike]
            
            # Calculate location of top corners
            panel = dict()
            panel['x1'] = xTopCenter + xTopOffset
            panel['y1'] = yTopCenter + yTopOffset
            panel['z1'] = zTopCenter - zTopOffset
            panel['x2'] = xTopCenter - xTopOffset
            panel['y2'] = yTopCenter - yTopOffset
            panel['z2'] = zTopCenter - zTopOffset

            # Calculate location of bottom corners
            panel['x3'] = xTopCenter + xTopBottomOffset + xTopOffset
            panel['y3'] = yTopCenter + yTopBottomOffset + yTopOffset
            panel['z3'] = zTopCenter + zTopBottomOffset - zTopOffset
            panel['x4'] = xTopCenter + xTopBottomOffset - xTopOffset
            panel['y4'] = yTopCenter + yTopBottomOffset - yTopOffset
            panel['z4'] = zTopCenter + zTopBottomOffset - zTopOffset
            
            # Extract fault slip
            panel['slip'] = F['seg' + str(iPanel) + 'SLIP'][0][iDownDip][iAlongStrike]

            # Extract patch dip, strike, width, and length
            panel['dip'] = F['seg' + str(iPanel) + 'DipAn']
            panel['strike'] = F['seg' + str(iPanel) + 'AStke']
            panel['rake'] = F['seg' + str(iPanel) + 'RAKE'][0][iDownDip][iAlongStrike]
            panel['angle'] = angle
            panel['width'] = W
            panel['length'] = L
            
            # Convert all distance measurements from kilometers to meters
            panel['x1'] = panel['x1'] * km2m
            panel['x2'] = panel['x2'] * km2m
            panel['x3'] = panel['x3'] * km2m
            panel['x4'] = panel['x4'] * km2m
            panel['y1'] = panel['y1'] * km2m
            panel['y2'] = panel['y2'] * km2m
            panel['y3'] = panel['y3'] * km2m
            panel['y4'] = panel['y4'] * km2m
            panel['z1'] = panel['z1'] * km2m
            panel['z2'] = panel['z2'] * km2m
            panel['z3'] = panel['z3'] * km2m
            panel['z4'] = panel['z4'] * km2m
            panel['width'] = panel['width'] * km2m
            panel['length'] = panel['length'] * km2m
            
            # Convert slip from centimeters to meters
            panel['slip'] = panel['slip'] * cm2m
            rTemp = np.array([[math.cos(math.radians(panel['rake'])), -math.sin(math.radians(panel['rake']))], [math.sin(math.radians(panel['rake'])), math.cos(math.radians(panel['rake']))]])
            xTempOrig = np.array([[panel['slip']], [0]])
            xTempRot = np.dot(rTemp, xTempOrig)
            panel['slipStrike'] = xTempRot[0]
            panel['slipDip'] = xTempRot[1]
            S.append(panel)

# Calculate elastic displacement field associated with one fault patch
xVec = np.linspace(-50*km2m, 50*km2m, N)
yVec = np.linspace(-50*km2m, 50*km2m, N)
xMat, yMat = np.meshgrid(xVec, yVec)
xVec = xMat.reshape(xMat.size, 1)
yVec = yMat.reshape(yMat.size, 1)
zVec = obsDepth*np.ones(xVec.size)
ux = np.zeros(xVec.size)
uy = np.zeros(xVec.size)
uz = np.zeros(xVec.size)
sxx = np.zeros(xVec.size)
sxy = np.zeros(xVec.size)
sxz = np.zeros(xVec.size)
syy = np.zeros(xVec.size)
syz = np.zeros(xVec.size)
szz = np.zeros(xVec.size)
cfs = np.zeros(xVec.size)
alpha = (lambdaLame+muLame) / (lambdaLame+2*muLame)

for iPatch in range(0, len(S)):
    print 'patch ' + str(iPatch+1) + ' of ' + str(len(S))
    # Loop over observation coordinates and calculate displacements and stresses for each source/observation pair
    for iObs in range(0, len(xVec)):
        # Translate and (un)rotate observation coordinates
        xTemp = xVec[iObs]-S[iPatch]['x1']
        yTemp = yVec[iObs]-S[iPatch]['y1']
        rTemp = np.array([[math.cos(math.radians(-S[iPatch]['angle'])), -math.sin(math.radians(-S[iPatch]['angle']))], [math.sin(math.radians(-S[iPatch]['angle'])), math.cos(math.radians(S[iPatch]['angle']))]])
        xTempOrig = np.array([xTemp, yTemp]) # no need for brackets around x/yTemp because they are already arrays
        xTempRot = np.dot(rTemp, xTempOrig)
        xTemp = xTempRot[0];
        yTemp = xTempRot[1];

        # Calculate elastic deformation using Okada's method
        # Seven arguments to DC3DWrapper are required:
        # alpha = (lambda + mu) / (lambda + 2 * mu)
        # xo = 3-vector representing the observation point (x, y, z in the original)
        # depth = the depth of the fault origin
        # dip = the dip-angle of the rectangular dislocation surface
        # strike_width = the along-strike range of the surface (al1,al2 in the original)
        # dip_width = the along-dip range of the surface (aw1, aw2 in the original)
        # dislocation = 3-vector representing the direction of motion on the surface (DISL1, DISL2, DISL3)
        success, u, uGrad = dc3dwrapper(alpha, [xTemp, yTemp, zVec[iObs]],
                                        S[iPatch]['z3'], S[iPatch]['dip'],
                                        [0.0, S[iPatch]['length']],
                                        [0.0, S[iPatch]['width']],
                                        [S[iPatch]['slipStrike'], S[iPatch]['slipDip'], 0.0])
        ux[iObs] = ux[iObs] + u[0]
        uy[iObs] = uy[iObs] + u[1]
        uz[iObs] = uz[iObs] + u[2]
        sxx[iObs] = sxx[iObs] + uGrad[0, 0]
        sxy[iObs] = sxy[iObs] + 0.5*(uGrad[0, 1] + uGrad[1, 0])
        sxz[iObs] = sxz[iObs] + 0.5*(uGrad[0, 2] + uGrad[2, 0])
        syy[iObs] = syy[iObs] + uGrad[1, 1]
        syz[iObs] = syz[iObs] + 0.5*(uGrad[1, 2] + uGrad[2, 1])
        szz[iObs] = szz[iObs] + uGrad[2, 2]
        
        # Resolve Coulomb failure stresses on reciever plane
        nVecInPlane = [0, 1, 0]
        nVecNormal = [1, 0, 0]
       
        deltaTau = (sxx[iObs] * nVecNormal[0] * nVecInPlane[0] + 
                    sxy[iObs] * nVecNormal[1] * nVecInPlane[0] + 
                    sxz[iObs] * nVecNormal[2] * nVecInPlane[0] + 
                    sxy[iObs] * nVecNormal[0] * nVecInPlane[1] + 
                    syy[iObs] * nVecNormal[1] * nVecInPlane[1] + 
                    syz[iObs] * nVecNormal[2] * nVecInPlane[1] + 
                    sxz[iObs] * nVecNormal[0] * nVecInPlane[2] + 
                    syz[iObs] * nVecNormal[1] * nVecInPlane[2] + 
                    szz[iObs] * nVecNormal[2] * nVecInPlane[2])
        deltaSigma = (sxx[iObs] * nVecNormal[0] * nVecNormal[0] + 
                      sxy[iObs] * nVecNormal[1] * nVecNormal[0] + 
                      sxz[iObs] * nVecNormal[2] * nVecNormal[0] + 
                      sxy[iObs] * nVecNormal[0] * nVecNormal[1] + 
                      syy[iObs] * nVecNormal[1] * nVecNormal[1] + 
                      syz[iObs] * nVecNormal[2] * nVecNormal[1] + 
                      sxz[iObs] * nVecNormal[0] * nVecNormal[2] + 
                      syz[iObs] * nVecNormal[1] * nVecNormal[2] + 
                      szz[iObs] * nVecNormal[2] * nVecNormal[2])
        cfs[iObs] = deltaTau - coefficientOfFriction * deltaSigma

# Reshape vectors as plain matrices for plotting
uxMat = np.reshape(ux, xMat.shape)
uyMat = np.reshape(uy, xMat.shape)
uzMat = np.reshape(uz, xMat.shape)
uHorizontalMagMat = np.sqrt(uxMat**2.0 + uyMat**2.0);
sxxMat = np.reshape(sxx, xMat.shape);
sxyMat = np.reshape(sxy, xMat.shape);
sxzMat = np.reshape(sxz, xMat.shape);
syyMat = np.reshape(syy, xMat.shape);
syzMat = np.reshape(syz, xMat.shape);
szzMat = np.reshape(szz, xMat.shape);

# Clip CFS values for plotting purposes
cfsHighIdx1 = (cfs>cfsUpperLimit).nonzero()
cfsHighIdx2 = (cfs>0).nonzero()
cfsHighIdx = np.intersect1d(np.array(cfsHighIdx1), np.array(cfsHighIdx2))
cfsLowIdx1 = (cfs<cfsLowerLimit).nonzero()
cfsLowIdx2 = (cfs<0).nonzero()
cfsLowIdx = np.intersect1d(np.array(cfsLowIdx1), np.array(cfsLowIdx2))
cfs[cfsHighIdx] = cfsUpperLimit
cfs[cfsLowIdx] = cfsLowerLimit
cfsMat = np.reshape(cfs, xMat.shape)

# Pick a fill color proportional to slip magnitude
# slipMin = min(S, key=lambda x: x['slip'])['slip'] # Don't need thus yet but could be useful
# slipMax = max(S, key=lambda x: x['slip'])['slip'] # Don't need thus yet but could be useful
# slipDiff = slipMax - slipMin
# nColors = 256;

fig = plt.figure()
ax2 = fig.gca()
origin = 'lower'
CS = plt.contourf(xMat, yMat, cfsMat, 10, cmap=cm.coolwarm, origin=origin, hold='on')
for iPatch in range(0, len(S)): # Plot the edges of each fault patch fault patches
    ax2.plot([S[iPatch]['x1'], S[iPatch]['x2']], [S[iPatch]['y1'], S[iPatch]['y2']], color='black')
    ax2.plot([S[iPatch]['x2'], S[iPatch]['x4']], [S[iPatch]['y2'], S[iPatch]['y4']], color='black')
    ax2.plot([S[iPatch]['x1'], S[iPatch]['x3']], [S[iPatch]['y1'], S[iPatch]['y3']], color='black')
    ax2.plot([S[iPatch]['x3'], S[iPatch]['x4']], [S[iPatch]['y3'], S[iPatch]['y4']], color='black')


plt.title(fileName)
plt.xlabel('x (m)')
plt.ylabel('y (m)')

# Make a colorbar for the ContourSet returned by the contourf call
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('CFS (Pa)')
plt.show()
