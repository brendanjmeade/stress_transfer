import numpy as np
import scipy.io as sio

# load file and extract geometric coordiantes and slip distribution
eventName = 's1999HECTOR01SALI'
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
F = sio.loadmat(eventName + '.mat')
F = F[eventName]

# extract the coordiantes for the rectangular fault patches
x = F[0]['seg1geoX']
y = F[0]['seg1geoY']
z = F[0]['seg1geoZ']

nPanel = F[0]['invSEGM'][0][0][0]

# % Extract the fault geometry and slip into a single stucture of ungrouped patches
# patchCount = 0;
# for iPanel = 1:nPanel
#     % Extract geometric parameters from this panel common to all patches
#     eval(sprintf('strike = F.seg%dAStke;', iPanel));
#     angle = -(strike-90);
#     angle(angle<0) = angle(angle<0)+360;
#     eval(sprintf('L = F.seg%dDimWL(2)/size(F.seg%dgeoX, 2);', iPanel, iPanel));
#     eval(sprintf('W = F.seg%dDimWL(2)/size(F.seg%dgeoX, 1);', iPanel, iPanel));
#     temp = [cosd(angle), -sind(angle); sind(angle), cosd(angle)]*[L/2; 0];
#     xTopOffset = temp(1);
#     yTopOffset = temp(2);
#     zTopOffset = 0;
#     eval(sprintf('xTopBottomOffset = F.seg%dgeoX(2, 1) - F.seg%dgeoX(1, 1);', iPanel, iPanel));
#     eval(sprintf('yTopBottomOffset = F.seg%dgeoY(2, 1) - F.seg%dgeoY(1, 1);', iPanel, iPanel));
#     eval(sprintf('zTopBottomOffset = F.seg%dgeoZ(2, 1) - F.seg%dgeoZ(1, 1);', iPanel, iPanel));

#     % For current panel get the number patches down-dip and along strike
#     eval(sprintf('[nDownDip, nAlongStrike] = size(F.seg%dgeoX);', iPanel));

#     % Loops over the down-dip and along-strike patches of the current panel
#     for iDownDip = 1:nDownDip
#         for iAlongStrike = 1:nAlongStrike
#             patchCount = patchCount + 1;
            
#             % Extract top center coordinates of current patch
#             eval(sprintf('xTopCenter = F.seg%dgeoX(iDownDip, iAlongStrike);', iPanel));
#             eval(sprintf('yTopCenter = F.seg%dgeoY(iDownDip, iAlongStrike);', iPanel));
#             eval(sprintf('zTopCenter = F.seg%dgeoZ(iDownDip, iAlongStrike);', iPanel));
            
#             % Calculate location of top corners
#             S(patchCount).x1 = xTopCenter + xTopOffset;
#             S(patchCount).y1 = yTopCenter + yTopOffset;
#             S(patchCount).z1 = zTopCenter - zTopOffset;
#             S(patchCount).x2 = xTopCenter - xTopOffset;
#             S(patchCount).y2 = yTopCenter - yTopOffset;
#             S(patchCount).z2 = zTopCenter - zTopOffset;
        
#             % Calculate location of bottom corners
#             S(patchCount).x3 = xTopCenter + xTopBottomOffset + xTopOffset;
#             S(patchCount).y3 = yTopCenter + yTopBottomOffset + yTopOffset;
#             S(patchCount).z3 = zTopCenter + zTopBottomOffset - zTopOffset;
#             S(patchCount).x4 = xTopCenter + xTopBottomOffset - xTopOffset;
#             S(patchCount).y4 = yTopCenter + yTopBottomOffset - yTopOffset;
#             S(patchCount).z4 = zTopCenter + zTopBottomOffset - zTopOffset;
            
#             % Extract fault slip
#             eval(sprintf('S(patchCount).slip = F.seg%dSLIP(iDownDip, iAlongStrike);', iPanel));

#             % Extract patch dip, strike, width, and length
#             eval(sprintf('S(patchCount).dip = F.seg%dDipAn;', iPanel));
#             eval(sprintf('S(patchCount).strike = F.seg%dAStke;', iPanel));
#             eval(sprintf('S(patchCount).rake = F.seg%dRAKE(iDownDip, iAlongStrike);', iPanel));
            
#             S(patchCount).angle = angle;
#             S(patchCount).width = W;
#             S(patchCount).length = L;
            
#             % Convert all distance measurements from kilometers to meters
#             S(patchCount).x1 = S(patchCount).x1 * km2m;
#             S(patchCount).x2 = S(patchCount).x2 * km2m;
#             S(patchCount).x3 = S(patchCount).x3 * km2m;
#             S(patchCount).x4 = S(patchCount).x4 * km2m;
#             S(patchCount).y1 = S(patchCount).y1 * km2m;
#             S(patchCount).y2 = S(patchCount).y2 * km2m;
#             S(patchCount).y3 = S(patchCount).y3 * km2m;
#             S(patchCount).y4 = S(patchCount).y4 * km2m;
#             S(patchCount).z1 = S(patchCount).z1 * km2m;
#             S(patchCount).z2 = S(patchCount).z2 * km2m;
#             S(patchCount).z3 = S(patchCount).z3 * km2m;
#             S(patchCount).z4 = S(patchCount).z4 * km2m;
#             S(patchCount).width = S(patchCount).width * km2m;
#             S(patchCount).length = S(patchCount).length * km2m;
            
#             % Convert slip from centimeters to meters
#             S(patchCount).slip = S(patchCount).slip * cm2m;
#             R = [cosd(S(patchCount).rake) , -sind(S(patchCount).rake) ; sind(S(patchCount).rake) , cosd(S(patchCount).rake)];
#             temp = R * [S(patchCount).slip ; 0];
#             S(patchCount).slipStrike = temp(1);
#             S(patchCount).slipDip = temp(2);
#         end
#     end
# end

# % Visualize the results
# figure;
# hold on;

# % Pick a fill color proportional to slip magnitude
# slipMin = min([S.slip]);
# slipMax = max([S.slip]);
# slipDiff = slipMax - slipMin;
# nColors = 256;
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

# % Calculate elastic displacement field associated with one fault patch
# xVec = linspace(-50*km2m, 50*km2m, N);
# yVec = linspace(-50*km2m, 50*km2m, N);
# [xMat, yMat] = meshgrid(xVec, yVec);
# xVec = xMat(:);
# yVec = yMat(:);
# zVec = obsDepth*ones(size(xVec));
# ux = zeros(size(xVec));
# uy = zeros(size(xVec));
# uz = zeros(size(xVec));
# sxx = zeros(size(xVec));
# sxy = zeros(size(xVec));
# sxz = zeros(size(xVec));
# syy = zeros(size(xVec));
# syz = zeros(size(xVec));
# szz = zeros(size(xVec));
# cfs = zeros(size(xVec));
# alpha = (lambda+mu) / (lambda+2*mu);

# for iPatch = 1:numel(S)
# % Loop over observation coordinates and calculate displacements and stresses for each source/observation pair
#     for iObs = 1:numel(xVec)
#         % Translate and (un)rotate observation coordinates
#         xtemp = xVec(iObs)-S(iPatch).x1;
#         ytemp = yVec(iObs)-S(iPatch).y1;
#         R = [cosd(-S(iPatch).angle) , -sind(-S(iPatch).angle) ; sind(-S(iPatch).angle) , cosd(-S(iPatch).angle)];
#         posTemp = R*[xtemp ; ytemp];
#         xtemp = posTemp(1);
#         ytemp = posTemp(2);
    
#         % Calculate elastic deformation using Okada's method
#         % Seven arguments to DC3DWrapper are required:
#         % alpha = (lambda + mu) / (lambda + 2 * mu)
#         % xo = 3-vector representing the observation point (x, y, z in the original)
#         % depth = the depth of the fault origin
#         % dip = the dip-angle of the rectangular dislocation surface
#         % strike_width = the along-strike range of the surface (al1,al2 in the original)
#         % dip_width = the along-dip range of the surface (aw1, aw2 in the original)
#         % dislocation = 3-vector representing the direction of motion on the surface (DISL1, DISL2, DISL3)
#         [success, u, uGrad] = DC3Dwrapper(alpha, ...
#                                           [xtemp, ytemp, zVec(iObs)], ...
#                                           S(iPatch).z3, S(iPatch).dip, ...
#                                           [0, S(iPatch).length], ...
#                                           [0, S(iPatch).width], ...
#                                           [S(iPatch).slipStrike, S(iPatch).slipDip, 0.0]);
#         ux(iObs) = ux(iObs) + u(1);
#         uy(iObs) = uy(iObs) + u(2);
#         uz(iObs) = uz(iObs) + u(3);
#         sxx(iObs) = sxx(iObs) + uGrad(1, 1);
#         sxy(iObs) = sxy(iObs) + 0.5*(uGrad(1, 2) + uGrad(2, 1));
#         sxz(iObs) = sxz(iObs) + 0.5*(uGrad(1, 3) + uGrad(3, 1));
#         syy(iObs) = syy(iObs) + uGrad(2, 2);
#         syz(iObs) = syz(iObs) + 0.5*(uGrad(2, 3) + uGrad(3, 2));
#         szz(iObs) = szz(iObs) + uGrad(3, 3);
        
#         % Resolve Coulomb failure stresses on reciever plane
#         nVecInPlane = [0 1 0];
#         nVecNormal = [1 0 0];
       
#         deltaTau = sxx(iObs) * nVecNormal(1) * nVecInPlane(1) + ...
#                    sxy(iObs) * nVecNormal(2) * nVecInPlane(1) + ...
#                    sxz(iObs) * nVecNormal(3) * nVecInPlane(1) + ...
#                    sxy(iObs) * nVecNormal(1) * nVecInPlane(2) + ...
#                    syy(iObs) * nVecNormal(2) * nVecInPlane(2) + ...
#                    syz(iObs) * nVecNormal(3) * nVecInPlane(2) + ...
#                    sxz(iObs) * nVecNormal(1) * nVecInPlane(3) + ...
#                    syz(iObs) * nVecNormal(2) * nVecInPlane(3) + ...
#                    szz(iObs) * nVecNormal(3) * nVecInPlane(3);
#         deltaSigma = sxx(iObs) * nVecNormal(1) * nVecNormal(1) + ...
#                      sxy(iObs) * nVecNormal(2) * nVecNormal(1) + ...
#                      sxz(iObs) * nVecNormal(3) * nVecNormal(1) + ...
#                      sxy(iObs) * nVecNormal(1) * nVecNormal(2) + ...
#                      syy(iObs) * nVecNormal(2) * nVecNormal(2) + ...
#                      syz(iObs) * nVecNormal(3) * nVecNormal(2) + ...
#                      sxz(iObs) * nVecNormal(1) * nVecNormal(3) + ...
#                      syz(iObs) * nVecNormal(2) * nVecNormal(3) + ...
#                      szz(iObs) * nVecNormal(3) * nVecNormal(3);
#         cfs(iObs) = deltaTau - coefficientOfFriction * deltaSigma;
#     end
# end

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

