import math
import numpy as np
import code
from scipy import io as sio
from okada_wrapper import dc3d0wrapper, dc3dwrapper
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import cm

def readSrcmodFile(fileName):
    # Constants for changing dimensions
    km2m = 1e3 # Convert kilometers to meters
    cm2m = 1e-2 # Convert centimeters to meters

    # Load the .mat (HDF5-ish) version of the model geometry and slip distribution
    F = sio.loadmat(fileName + '.mat')
    F = F[fileName]
    F = F[0]

    # Declare fields to populate
    EventSrcmod = dict() # Empty dictionary
    EventSrcmod['angle'] = []
    EventSrcmod['dip'] = []
    EventSrcmod['length'] = []
    EventSrcmod['rake'] = []
    EventSrcmod['slip'] = []
    EventSrcmod['slipDip'] = []
    EventSrcmod['slipStrike'] = []
    EventSrcmod['strike'] = []
    EventSrcmod['width'] = []
    EventSrcmod['x1'] = []
    EventSrcmod['x2'] = []
    EventSrcmod['x3'] = []
    EventSrcmod['x4'] = []
    EventSrcmod['y1'] = []
    EventSrcmod['y2'] = []
    EventSrcmod['y3'] = []
    EventSrcmod['y4'] = []
    EventSrcmod['z1'] = []
    EventSrcmod['z2'] = []
    EventSrcmod['z3'] = []
    EventSrcmod['z4'] = []

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
            
                # Calculate location of top corners and convert from km to m
                EventSrcmod['x1'].append(km2m * (xTopCenter + xTopOffset))
                EventSrcmod['y1'].append(km2m * (yTopCenter + yTopOffset))
                EventSrcmod['z1'].append(km2m * (zTopCenter - zTopOffset)) # not sure this is right?!
                EventSrcmod['x2'].append(km2m * (xTopCenter - xTopOffset))
                EventSrcmod['y2'].append(km2m * (yTopCenter - yTopOffset))
                EventSrcmod['z2'].append(km2m * (zTopCenter - zTopOffset))

                # Calculate location of bottom corners and convert from km to m
                EventSrcmod['x3'].append(km2m * (xTopCenter + xTopBottomOffset + xTopOffset))
                EventSrcmod['y3'].append(km2m * (yTopCenter + yTopBottomOffset + yTopOffset))
                EventSrcmod['z3'].append(km2m * (zTopCenter + zTopBottomOffset - zTopOffset))
                EventSrcmod['x4'].append(km2m * (xTopCenter + xTopBottomOffset - xTopOffset))
                EventSrcmod['y4'].append(km2m * (yTopCenter + yTopBottomOffset - yTopOffset))
                EventSrcmod['z4'].append(km2m * (zTopCenter + zTopBottomOffset - zTopOffset))
            
                # Extract patch dip, strike, width, and length
                EventSrcmod['dip'].append(F['seg' + str(iPanel) + 'DipAn'])
                EventSrcmod['strike'].append(F['seg' + str(iPanel) + 'AStke'])
                EventSrcmod['rake'].append(F['seg' + str(iPanel) + 'RAKE'][0][iDownDip][iAlongStrike])
                EventSrcmod['angle'].append(angle)
                EventSrcmod['width'].append(km2m * W)
                EventSrcmod['length'].append(km2m * L)
            
                # Extract fault slip
                EventSrcmod['slip'].append(cm2m*(F['seg' + str(iPanel) + 'SLIP'][0][iDownDip][iAlongStrike]))
                rTemp = np.array([[math.cos(math.radians(EventSrcmod['rake'][-1])), 
                                   -math.sin(math.radians(EventSrcmod['rake'][-1]))],
                                  [math.sin(math.radians(EventSrcmod['rake'][-1])),
                                   math.cos(math.radians(EventSrcmod['rake'][-1]))]])
                xTempOrig = np.array([[EventSrcmod['slip'][-1]], [0]])
                xTempRot = np.dot(rTemp, xTempOrig)
                EventSrcmod['slipStrike'].append(xTempRot[0])
                EventSrcmod['slipDip'].append(xTempRot[1])
    return(EventSrcmod)

def calcCfs(StressTensor, nVecNormal, nVecInPlane, coefficientOfFriction):
    cfs = np.zeros(np.shape(StressTensor['sxx']))
    for iObs in range(0, len(StressTensor['sxx'])):
        deltaTau = (StressTensor['sxx'][iObs] * nVecNormal[0] * nVecInPlane[0] + 
                    StressTensor['sxy'][iObs] * nVecNormal[1] * nVecInPlane[0] + 
                    StressTensor['sxz'][iObs] * nVecNormal[2] * nVecInPlane[0] + 
                    StressTensor['sxy'][iObs] * nVecNormal[0] * nVecInPlane[1] + 
                    StressTensor['syy'][iObs] * nVecNormal[1] * nVecInPlane[1] + 
                    StressTensor['syz'][iObs] * nVecNormal[2] * nVecInPlane[1] + 
                    StressTensor['sxz'][iObs] * nVecNormal[0] * nVecInPlane[2] + 
                    StressTensor['syz'][iObs] * nVecNormal[1] * nVecInPlane[2] + 
                    StressTensor['szz'][iObs] * nVecNormal[2] * nVecInPlane[2])
        deltaSigma = (StressTensor['sxx'][iObs] * nVecNormal[0] * nVecNormal[0] + 
                      StressTensor['sxy'][iObs] * nVecNormal[1] * nVecNormal[0] + 
                      StressTensor['sxz'][iObs] * nVecNormal[2] * nVecNormal[0] + 
                      StressTensor['sxy'][iObs] * nVecNormal[0] * nVecNormal[1] + 
                      StressTensor['syy'][iObs] * nVecNormal[1] * nVecNormal[1] + 
                      StressTensor['syz'][iObs] * nVecNormal[2] * nVecNormal[1] + 
                      StressTensor['sxz'][iObs] * nVecNormal[0] * nVecNormal[2] + 
                      StressTensor['syz'][iObs] * nVecNormal[1] * nVecNormal[2] + 
                      StressTensor['szz'][iObs] * nVecNormal[2] * nVecNormal[2])
        cfs[iObs] = deltaTau - coefficientOfFriction * deltaSigma
    return(cfs)

def calcOkadaDisplacementStress(xVec, yVec, zVec, EventSrcmod, alpha):
    # Calculate elastic displacement field associated with one fault patch
    DisplacementVector = dict()
    DisplacementVector['ux'] = np.zeros(xVec.size)
    DisplacementVector['uy'] = np.zeros(xVec.size)
    DisplacementVector['uz'] = np.zeros(xVec.size)
    StressTensor = dict()
    StressTensor['sxx'] = np.zeros(xVec.size)
    StressTensor['sxy'] = np.zeros(xVec.size)
    StressTensor['sxz'] = np.zeros(xVec.size)
    StressTensor['syy'] = np.zeros(xVec.size)
    StressTensor['syz'] = np.zeros(xVec.size)
    StressTensor['szz'] = np.zeros(xVec.size)

    for iPatch in range(0, len(EventSrcmod['x1'])): # Loop over source patches
        print 'patch ' + str(iPatch+1) + ' of ' + str(len(EventSrcmod['x1']))
        # Loop over observation coordinates
        for iObs in range(0, len(xVec)):
            # Translate and (un)rotate observation coordinates
            xTemp = xVec[iObs]-EventSrcmod['x1'][iPatch]
            yTemp = yVec[iObs]-EventSrcmod['y1'][iPatch]
            rTemp = np.array([[math.cos(math.radians(-EventSrcmod['angle'][iPatch])), 
                              -math.sin(math.radians(-EventSrcmod['angle'][iPatch]))],
                              [math.sin(math.radians(-EventSrcmod['angle'][iPatch])),
                               math.cos(math.radians(EventSrcmod['angle'][iPatch]))]])
            xTempOrig = np.array([xTemp, yTemp])
            xTempRot = np.dot(rTemp, xTempOrig)
            xTemp = xTempRot[0]
            yTemp = xTempRot[1]

            # Calculate elastic deformation using Okada 1992 (BSSA)
            # Seven arguments to DC3DWrapper are required:
            # alpha = (lambda + mu) / (lambda + 2 * mu)
            # xo = 3-vector representing the observation point (x, y, z in the original)
            # depth = the depth of the fault origin
            # dip = the dip-angle of the rectangular dislocation surface
            # strike_width = the along-strike range of the surface (al1,al2 in the original)
            # dip_width = the along-dip range of the surface (aw1, aw2 in the original)
            # dislocation = 3-vector representing the direction of motion on the surface (DISL1, DISL2, DISL3)
            success, u, uGrad = dc3dwrapper(alpha, [xTemp, yTemp, zVec[iObs]],
                                            EventSrcmod['z3'][iPatch], EventSrcmod['dip'][iPatch],
                                            [0.0, EventSrcmod['length'][iPatch]],
                                            [0.0, EventSrcmod['width'][iPatch]],
                                            [EventSrcmod['slipStrike'][iPatch], EventSrcmod['slipDip'][iPatch], 0.0])
            DisplacementVector['ux'][iObs] = DisplacementVector['ux'][iObs] + u[0]
            DisplacementVector['uy'][iObs] = DisplacementVector['uy'][iObs] + u[1]
            DisplacementVector['uz'][iObs] = DisplacementVector['uz'][iObs] + u[2]
            StressTensor['sxx'][iObs] = StressTensor['sxx'][iObs] + uGrad[0, 0]
            StressTensor['sxy'][iObs] = StressTensor['sxy'][iObs] + 0.5*(uGrad[0, 1] + uGrad[1, 0])
            StressTensor['sxz'][iObs] = StressTensor['sxz'][iObs] + 0.5*(uGrad[0, 2] + uGrad[2, 0])
            StressTensor['syy'][iObs] = StressTensor['syy'][iObs] + uGrad[1, 1]
            StressTensor['syz'][iObs] = StressTensor['syz'][iObs] + 0.5*(uGrad[1, 2] + uGrad[2, 1])
            StressTensor['szz'][iObs] = StressTensor['szz'][iObs] + uGrad[2, 2]

    return(DisplacementVector, StressTensor)


def main():
    # Name of Srcmod file to read
    fileName = 's1999HECTOR01SALI'

    # Parameters for CFS calculation and visualization
    lambdaLame = 0.25 # First Lame parameter
    muLame = 0.25 # shear modulus
    alpha = (lambdaLame+muLame) / (lambdaLame+2*muLame)
    coefficientOfFriction = 0.4 # Coefficient of friction
    obsDepth = -5e3; # depth of observation coordinates
    nVecInPlane = [0, 1, 0]
    nVecNormal = [1, 0, 0]
    cfsUpperLimit = 5e-6; # for visualziation purposes
    cfsLowerLimit = -5e-6; # for visualization purposes

    # Observation coordinates
    N = 20 # Number of grid points in x and y-directions for visualization
    xVec = np.linspace(-50e3, 50e3, N)
    yVec = np.linspace(-50e3, 50e3, N)
    xMat, yMat = np.meshgrid(xVec, yVec)
    xVec = xMat.reshape(xMat.size, 1)
    yVec = yMat.reshape(yMat.size, 1)
    zVec = obsDepth*np.ones(xVec.size)

    # Read in Srcmod fault geometry and slip distribution for this representation of the event
    EventSrcmod = readSrcmodFile(fileName)

    # Calculate displacement vector and stress tensor at observation coordinates
    DisplacementVector, StressTensor = calcOkadaDisplacementStress(xVec, yVec, zVec, EventSrcmod, alpha)

    # Resolve Coulomb failure stresses on reciever plane
    cfs = calcCfs(StressTensor, nVecNormal, nVecInPlane, coefficientOfFriction)

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

    # Generate figure showing fault geometry and CFS feild
    fig = plt.figure(facecolor='white')
    ax = fig.gca()
    cs = plt.contourf(xMat, yMat, cfsMat, 10, cmap=cm.coolwarm, origin='lower', hold='on')
    for iPatch in range(0, len(EventSrcmod['x1'])): # Plot the edges of each fault patch fault patches
        ax.plot([EventSrcmod['x1'][iPatch], EventSrcmod['x2'][iPatch]], [EventSrcmod['y1'][iPatch], EventSrcmod['y2'][iPatch]], color='black')
        ax.plot([EventSrcmod['x2'][iPatch], EventSrcmod['x4'][iPatch]], [EventSrcmod['y2'][iPatch], EventSrcmod['y4'][iPatch]], color='black')
        ax.plot([EventSrcmod['x1'][iPatch], EventSrcmod['x3'][iPatch]], [EventSrcmod['y1'][iPatch], EventSrcmod['y3'][iPatch]], color='black')
        ax.plot([EventSrcmod['x3'][iPatch], EventSrcmod['x4'][iPatch]], [EventSrcmod['y3'][iPatch], EventSrcmod['y4'][iPatch]], color='black')

    plt.title(fileName)
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    #plt.axis('off') # turning off axes labels...will replace with scale bar
    cbar = plt.colorbar(cs) # Make a colorbar for the ContourSet returned by the contourf call
    cbar.ax.set_ylabel('CFS (Pa)')
    plt.show()

    # Provide keyboard control to interact with variables
    code.interact(local=locals())

if __name__ == "__main__":
   main()
