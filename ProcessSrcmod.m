close all;
clear all;
faultName = 's1999HECTOR01SALI';
F = load(strcat(faultName, '.mat'));
F = F.s1999HECTOR01SALI;

% Count how many panels there are
names = fieldnames(F);
nPanels = 0;
for i = 1:numel(names)
    currentField = char(names(i));
    if numel(currentField) == 8
        if strcmp(currentField(5:8), 'geoZ')
            nPanels = nPanels + 1;
        end
    end
end

patchCount = 0;
for k = 1:nPanels
    % Extract geometric parameters from this panel common to all patches
    eval(sprintf('strike = F.seg%dAStke;', k));
    angle = -(strike-90);
    angle(angle<0) = angle(angle<0)+360;
    eval(sprintf('L = F.seg%dDimWL(2)/size(F.seg%dgeoX, 2);', k, k));
    eval(sprintf('W = F.seg%dDimWL(2)/size(F.seg%dgeoX, 1);', k, k));
    temp = [cosd(angle), -sind(angle); sind(angle), cosd(angle)]*[L/2; 0];
    xTopOffset = temp(1);
    yTopOffset = temp(2);
    zTopOffset = 0;
    eval(sprintf('xTopBottomOffset = F.seg%dgeoX(2, 1) - F.seg%dgeoX(1, 1);', k, k));
    eval(sprintf('yTopBottomOffset = F.seg%dgeoY(2, 1) - F.seg%dgeoY(1, 1);', k, k));
    eval(sprintf('zTopBottomOffset = F.seg%dgeoZ(2, 1) - F.seg%dgeoZ(1, 1);', k, k));

    % For current panel get the number patches down-dip and along strike
    eval(sprintf('[nDownDip, nAlongStrike] = size(F.seg%dgeoX);', k));

    % Loops over the down-dip and along-strike patches of the current panel
    for i = 1:nDownDip
        for j = 1:nAlongStrike
            patchCount = patchCount + 1;
            
            % Extract top center coordinates of current patch
            eval(sprintf('xTopCenter = F.seg%dgeoX(i, j);', k));
            eval(sprintf('yTopCenter = F.seg%dgeoY(i, j);', k));
            eval(sprintf('zTopCenter = F.seg%dgeoZ(i, j);', k));
            
            % Calculate location of top corners
            S(patchCount).x1 = xTopCenter + xTopOffset;
            S(patchCount).y1 = yTopCenter + yTopOffset;
            S(patchCount).z1 = zTopCenter - zTopOffset;
            S(patchCount).x2 = xTopCenter - xTopOffset;
            S(patchCount).y2 = yTopCenter - yTopOffset;
            S(patchCount).z2 = zTopCenter - zTopOffset;
        
            % Calculate location of bottom corners
            S(patchCount).x3 = xTopCenter + xTopBottomOffset + xTopOffset;
            S(patchCount).y3 = yTopCenter + yTopBottomOffset + yTopOffset;
            S(patchCount).z3 = zTopCenter + zTopBottomOffset - zTopOffset;
            S(patchCount).x4 = xTopCenter + xTopBottomOffset - xTopOffset;
            S(patchCount).y4 = yTopCenter + yTopBottomOffset - yTopOffset;
            S(patchCount).z4 = zTopCenter + zTopBottomOffset - zTopOffset;
            
            % Extract fault slip
            eval(sprintf('S(patchCount).slip = F.seg%dSLIP(i, j);', k));

            % Extract patch dip, strike, width, and length
            eval(sprintf('S(patchCount).dip = F.seg%dDipAn;', k));
            eval(sprintf('S(patchCount).strike = F.seg%dAStke;', k));
            S(patchCount).angle = angle;
            S(patchCount).width = W;
            S(patchCount).length = L;
            
        end
    end
end

% Visualize the results
figure;
hold on;

% Loop over all of the fault patches and plot
for i = 1:numel(S)
    % pick a fill color proportional to slip magnitude
    slipMin = min([S.slip]);
    slipMax = max([S.slip]);
    slipDiff = slipMax - slipMin;
    nColors = 256;
    cmap = flipud(gray(nColors));
    slipColorIdx = round((S(i).slip - slipMin)/slipDiff * nColors);
    if slipColorIdx < 1
        slipColorIdx = 1;
    end
    if slipColorIdx > nColors
        slipColorIdx = nColors
    end
    
    % Plot the individual fault patches
    fh = fill3([S(i).x1 S(i).x2 S(i).x4 S(i).x3], [S(i).y1 S(i).y2 S(i).y4 S(i).y3], -[S(i).z1 S(i).z2 S(i).z4 S(i).z3], 'y');
    set(fh, 'FaceColor', cmap(slipColorIdx, :))
end


% Add labels and set axes properties
axis equal;
box on;
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
view(3);
cameratoolbar;
camproj('perspective');

% Calculate elastic displacement field associated with one fault patch
N = 100;
xvec = linspace(-50, 50, N);
yvec = linspace(-50, 50, N);
[xmat, ymat] = meshgrid(xvec, yvec);
xvec = xmat(:);
yvec = ymat(:);
zvec = 0.0*ones(size(xvec));
ux = zeros(size(xvec));
uy = zeros(size(xvec));
uz = zeros(size(xvec));
sxx = zeros(size(xvec));
sxy = zeros(size(xvec));
sxz = zeros(size(xvec));
syy = zeros(size(xvec));
syz = zeros(size(xvec));
szz = zeros(size(xvec));
lambda = 0.25;
mu = 0.25;
alpha = (lambda+mu) / (lambda+2*mu);


for iPatch = 1:numel(S)
% Loop over observation coordinates and calculate displacements and stresses for each source/observation pair
    for iObs = 1:numel(xvec)
        % Translate and (un)rotate observation coordinates
        xtemp = xvec(iObs)-S(iPatch).x1;
        ytemp = yvec(iObs)-S(iPatch).y1;
        R = [cosd(-S(iPatch).angle) , -sind(-S(iPatch).angle) ; sind(-S(iPatch).angle) , cosd(-S(iPatch).angle)];
        posTemp = R*[xtemp ; ytemp];
        xtemp = posTemp(1);
        ytemp = posTemp(2);
    
        % Calculate elastic deformation using Okada's method
        % Seven arguments to DC3DWrapper are required:
        % alpha = (lambda + mu) / (lambda + 2 * mu)
        % xo = 3-vector representing the observation point (x, y, z in the original)
        % depth = the depth of the fault origin
        % dip = the dip-angle of the rectangular dislocation surface
        % strike_width = the along-strike range of the surface (al1,al2 in the original)
        % dip_width = the along-dip range of the surface (aw1, aw2 in the original)
        % dislocation = 3-vector representing the direction of motion on the surface (DISL1, DISL2, DISL3)
        [success, u, uGrad] = DC3Dwrapper(alpha, ...
                                          [xtemp, ytemp, zvec(iObs)], ...
                                          S(iPatch).z3, S(iPatch).dip, ...
                                          [0, S(iPatch).length], ...
                                          [0, S(iPatch).width], ...
                                          [1.0, 0.0, 0.0]);
        ux(iObs) = ux(iObs) + u(1);
        uy(iObs) = uy(iObs) + u(2);
        uz(iObs) = uz(iObs) + u(3);
        sxx(iObs) = sxx(iObs) + uGrad(1, 1);
        sxy(iObs) = sxy(iObs) + 0.5*(uGrad(1, 2) + uGrad(2, 1));
        sxz(iObs) = sxz(iObs) + 0.5*(uGrad(1, 3) + uGrad(3, 1));
        syy(iObs) = syy(iObs) + uGrad(2, 2);
        syz(iObs) = syz(iObs) + 0.5*(uGrad(2, 3) + uGrad(3, 2));
        szz(iObs) = szz(iObs) + uGrad(3, 3);
    end
end

uxmat = reshape(ux, size(xmat));
uymat = reshape(uy, size(xmat));
uzmat = reshape(uz, size(xmat));
sxxmat = reshape(uz, size(xmat));
sxymat = reshape(uz, size(xmat));
sxzmat = reshape(uz, size(xmat));
syymat = reshape(uz, size(xmat));
syzmat = reshape(uz, size(xmat));
szzmat = reshape(szz, size(xmat));

uHorizontalMag = sqrt(uxmat.^2 + uymat.^2);

% Plot a horizontal slice showing the magnitude of the horizontal displacement field
sh = surf(xmat, ymat, -zvec(1)*ones(size(uxmat)), log10(uHorizontalMag));
set(sh, 'EdgeColor', 'none');
set(sh, 'FaceAlpha', 0.85)
colormap(flipud(hot(20)));
camva(45);
campos([-143.5170  -87.9473   80.5196]);
axis tight

% Add a custom colorbar
ch = colorbar('horizontal');
set(ch, 'Position', [0.05 0.05 0.2 0.01]);
%set(get(ch, 'Title'), 'String', 'log_{10} (u_x^2 + u_y^2)^{1/2}')
%set(get(ch, 'Title'), 'VerticalAlignment', 'bottom');

