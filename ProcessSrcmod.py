import math
import numpy as np
import scipy
from okada_wrapper import dc3d0wrapper, dc3dwrapper

# load file and extract geometric coordiantes and slip distribution
eventName = 's1999HECTOR01SALI'
N = 50 # Number of grid points in x and y-directions for visualization
lambdaLame = 0.25 # First Lame parameter
muLame = 0.25 # shear modulus
coefficientOfFriction = 0.4 # Coefficient of friction
km2m = 1e3 # Convert kilometers to meters
cm2m = 1e-2 # Convert centimeters to meters
cfsUpperLimit = 5e-6; # for visualziation purposes
cfsLowerLimit = -5e-6; # for visualization purposes
obsDepth = -5e3; # depth of observation coordinates

# Load the .mat (HDF5-ish) version of the model geometry and slip distribution
F = scipy.io.loadmat(eventName + '.mat')
F = F[eventName]
F = F[0]
S = list() # Empty list

# Extract the fault geometry and slip into a single stucture of ungrouped patches
nPanel = int(F['invSEGM'][0][0][0]) # Count number of panels
for iPanel in range(1, nPanel + 1): # This index starts at 1 because of naming convention
    # Extract geometric parameters from this panel common to all patches
    strike = 0
    strike = F['seg' + str(iPanel) + 'AStke'][0][0][0]
    angle = -(strike-90)
    if angle < 0:
        angle = angle + 360

    # calculate the length and wide if individual patch elements in current panel
    L = F['seg' + str(iPanel) + 'DimWL'][0][0][1] / np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[1]
    W = F['seg' + str(iPanel) + 'DimWL'][0][0][1] / np.shape(F['seg' + str(iPanel) + 'geoX'][0][:][:])[0]
    rTemp = np.array([[math.cos(math.degrees(angle)), -math.sin(math.degrees(angle))], [math.sin(math.degrees(angle)), math.cos(math.degrees(angle))]])
    xTempOrig = np.array([[L/2], [0]])
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
            element = dict()
            element['x1'] = xTopCenter + xTopOffset
            element['y1'] = yTopCenter + yTopOffset
            element['z1'] = zTopCenter - zTopOffset
            element['x2'] = xTopCenter - xTopOffset
            element['y2'] = yTopCenter - yTopOffset
            element['z2'] = zTopCenter - zTopOffset

            # Calculate location of bottom corners
            element['x3'] = xTopCenter + xTopBottomOffset + xTopOffset
            element['y3'] = yTopCenter + yTopBottomOffset + yTopOffset
            element['z3'] = zTopCenter + zTopBottomOffset - zTopOffset
            element['x4'] = xTopCenter + xTopBottomOffset - xTopOffset
            element['y4'] = yTopCenter + yTopBottomOffset - yTopOffset
            element['z4'] = zTopCenter + zTopBottomOffset - zTopOffset
            
            # Extract fault slip
            element['slip'] = F['seg' + str(iPanel) + 'SLIP'][0][iDownDip][iAlongStrike]

            # Extract patch dip, strike, width, and length
            element['dip'] = F['seg' + str(iPanel) + 'DipAn']
            element['strike'] = F['seg' + str(iPanel) + 'AStke']
            element['rake'] = F['seg' + str(iPanel) + 'RAKE'][0][iDownDip][iAlongStrike]
            element['angle'] = angle
            element['width'] = W
            element['length'] = L
            
            # Convert all distance measurements from kilometers to meters
            element['x1'] = element['x1'] * km2m
            element['x2'] = element['x2'] * km2m
            element['x3'] = element['x3'] * km2m
            element['x4'] = element['x4'] * km2m
            element['y1'] = element['y1'] * km2m
            element['y2'] = element['y2'] * km2m
            element['y3'] = element['y3'] * km2m
            element['y4'] = element['y4'] * km2m
            element['z1'] = element['z1'] * km2m
            element['z2'] = element['z2'] * km2m
            element['z3'] = element['z3'] * km2m
            element['z4'] = element['z4'] * km2m
            element['width'] = element['width'] * km2m
            element['length'] = element['length'] * km2m
            
            # Convert slip from centimeters to meters
            element['slip'] = element['slip'] * cm2m
            rTemp = np.array([[math.cos(math.degrees(element['rake'])), -math.sin(math.degrees(element['rake']))], [math.sin(math.degrees(element['rake'])), math.cos(math.degrees(element['rake']))]])
            xTempOrig = np.array([[element['slip']], [0]])
            xTempRot = np.dot(rTemp, xTempOrig)
            element['slipStrike'] = xTempRot[0]
            element['slipDip'] = xTempRot[1]
            S.append(element)

# Calculate elastic displacement field associated with one fault patch
xVec = np.linspace(-50*km2m, 50*km2m, N)
yVec = np.linspace(-50*km2m, 50*km2m, N)
[xMat, yMat] = np.meshgrid(xVec, yVec)
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
    print iPatch
    # Loop over observation coordinates and calculate displacements and stresses for each source/observation pair
    for iObs in range(0, len(xVec)):
        # Translate and (un)rotate observation coordinates
        xTemp = xVec[iObs]-S[iPatch]['x1']
        yTemp = yVec[iObs]-S[iPatch]['y1']
        rTemp = np.array([[math.cos(math.degrees(-S[iPatch]['angle'])), -math.sin(math.degrees(-S[iPatch]['angle']))], [math.sin(math.degrees(-S[iPatch]['angle'])), math.cos(math.degrees(S[iPatch]['angle']))]])
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


# uxMat = reshape(ux, size(xMat));
# uyMat = reshape(uy, size(xMat));
# uzMat = reshape(uz, size(xMat));
# uHorizontalMagMat = sqrt(uxMat.^2 + uyMat.^2);
# sxxMat = reshape(uz, size(xMat));
# sxyMat = reshape(uz, size(xMat));
# sxzMat = reshape(uz, size(xMat));
# syyMat = reshape(uz, size(xMat));
# syzMat = reshape(uz, size(xMat));
# szzMat = reshape(szz, size(xMat));

# % Clip CFS values for plotting purposes
# cfsHighIdx1 = find(cfs>cfsUpperLimit);
# cfsHighIdx2 = find(cfs>0);
# cfsHighIdx = intersect(cfsHighIdx1, cfsHighIdx2);
# cfsLowIdx1 = find(cfs<cfsLowerLimit);
# cfsLowIdx2 = find(cfs<0);
# cfsLowIdx = intersect(cfsLowIdx1, cfsLowIdx2);
# cfs(cfsHighIdx) = cfsUpperLimit;
# cfs(cfsLowIdx) = cfsLowerLimit;
# cfsMat = reshape(cfs, size(xMat));


# Visualize the results
# figure;
# hold on;

# Pick a fill color proportional to slip magnitude
slipMin = min(S, key=lambda x: x['slip'])['slip']
slipMax = max(S, key=lambda x: x['slip'])['slip']
slipDiff = slipMax - slipMin
nColors = 256;
# cmap = flipud(gray(nColors));

# % Loop over all of the fault patches and plot
# for iPatch = 1:numel(S)
#     slipColorIdx = round((S(iPatch).slip - slipMin)/slipDiff * nColors);
#     if slipColorIdx < 1
#         slipColorIdx = 1;
#     end
#     if slipColorIdx > nColors
#         slipColorIdx = nColors
#     end
    
#     % Plot the individual fault patches
#     fh = fill3([S(iPatch).x1 S(iPatch).x2 S(iPatch).x4 S(iPatch).x3], [S(iPatch).y1 S(iPatch).y2 S(iPatch).y4 S(iPatch).y3], -[S(iPatch).z1 S(iPatch).z2 S(iPatch).z4 S(iPatch).z3], 'y');
#     set(fh, 'FaceColor', cmap(slipColorIdx, :));
# end

# % Add labels and set axes properties
# axis equal;
# box on;
# xlabel('x (m)');
# ylabel('y (m)');
# zlabel('z (m)');
# view(3);
# cameratoolbar;

# % plot a horizontal slice showing the magnitude of the horizontal displacement field
# sh = surf(xMat, yMat, obsDepth*ones(size(uxMat)), cfsMat);
# set(sh, 'EdgeColor', 'none');
# set(sh, 'FaceAlpha', 0.65)
# colormap(flipud(hot(20)));
# axis tight;

# % Add a small colorbar
# ch = colorbar('horizontal');
# set(ch, 'Position', [0.05 0.10 0.2 0.01]);
# colormap(bluewhitered);

