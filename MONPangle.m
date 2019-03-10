%% MONPangle
% Lloyd Russell 2018
% Calculate how to align optical axis of microscope opbjective to be 
% perpendicular to the tilted sample plane. Calculates the pitch and yaw
% values for use with the Bruker orbital nosepiece given the recorded tilt
% of window which is measured by recording the Z depth when three corners
% of the window/brain/sample comes into focus.

%% window coordinates (enter measurements here (negative for down))
topLeftZ    = 0;
topRightZ   = -80;
bottomLeftZ = +9;
FOVsize     = 433.8;  % distance between points measured (assuming square) % was 520.7

%% assemble coordinates
% x y coords of measured points
P0 = [-1 -1 0] * (FOVsize/2);  % top left
P1 = [1 -1 0] * (FOVsize/2);  % top right
P2 = [-1 1 0] * (FOVsize/2);  % bottom left

% Z coords
P0(3) = topLeftZ;  % top left
P1(3) = topRightZ;  % top right
P2(3) = bottomLeftZ;  % bottom left

%% calculate normal vector of plane
normal = cross(P0-P1, P0-P2);
normal = normal / norm( normal );  % unit vector

% plane coefficients
A = normal(1);
B = normal(2);
C = normal(3);
D = -((A*P0(1)) + (B*P0(2)) + (C*P0(3)));  % D = -dot(normal, P0);

% calculate the angle of tilt
u = normal;
v = [0 0 1];  % zaxis
angleFromStraight = atan2d(norm(cross(u,v)),dot(u,v));
v = [1 0 0];
angleX = 90-atan2d(norm(cross(u,v)),dot(u,v));
v = [0 1 0];
angleY = 90-atan2d(norm(cross(u,v)),dot(u,v));

%% model objective as 3 points, with 2 points of rotation (pitch and yaw)
% usable pitch and yaw range in practice
pitch_range = [-15 : 0.1 : 15];
yaw_range = [-180 : 0.1 : 0];

% objective model
O0 = [0 0 0];
O1 = [1 0 0];
O2 = [1 0 -1];

% rotate pitch
pitchTemplate = nan(numel(pitch_range), 3);
for i = 1:numel(pitch_range)
    theta = pitch_range(i);
    R = [1 0           0;
        0 cosd(theta) -sind(theta);
        0 sind(theta) cosd(theta)];
    pitchTemplate(i,:) = O2 * R;
end

% rotate all yaw
allO2 = nan(numel(yaw_range), numel(pitch_range), 3);
allO1 = nan(numel(yaw_range), 3);
for i = 1:numel(yaw_range)
    theta = yaw_range(i);
    R = [cosd(theta) -sind(theta) 0;
        sind(theta) cosd(theta)  0;
        0           0            1];
    allO2(i,:,:) = pitchTemplate * R;
    allO1(i,:) = O1 * R;
end


% calculate direction vectors or objective and test for parallel-ness
allDirectionVectors = nan(numel(yaw_range), numel(pitch_range), 3);
parallel = nan(numel(yaw_range), numel(pitch_range));
errorAngle = nan(numel(yaw_range), numel(pitch_range));
for i = 1:numel(yaw_range)
    for j = 1:numel(pitch_range)
        allDirectionVectors(i,j,:) = squeeze(allO2(i,j,:))' - allO1(i,:);
        v = squeeze(allDirectionVectors(i,j,:));
        parallel(i,j) = dot(v, normal);
        errorAngle(i,j) = 180 - atan2d(norm(cross(v,normal)), dot(v,normal));  % subtract from 180 because direction vectors are in opposite directions
    end
end

% direction of vectors irrelevant in this case, take absolute
parallel = abs(parallel);

% find vector closest to plane normal then take the corresponding pitch and yaw values
[accuracy,idx] = max(parallel(:));
[yaw_idx,pitch_idx] = ind2sub([numel(yaw_range), numel(pitch_range)],idx);

yaw = yaw_range(yaw_idx);
pitch = pitch_range(pitch_idx);
errorFinal = errorAngle(idx);


%% plot the window tilt and normal vector
% vector lines
normalPlot = (normal * (FOVsize/2));  % line starting at 0,0,0
straightPlot = (v * (FOVsize/2));

% plane surface
planeXLim = [-1 0 1] * (FOVsize/2);
planeYLim = [-1 0 1] * (FOVsize/2);
[X,Y] = meshgrid(planeXLim,planeYLim);
Z = (A*X + B*Y + D) / -C;    % a(x) + b(y) + c(z) + d = 0
middleZ = Z(2,2);

figure
surf(X,Y,Z, 'facecolor','w', 'edgecolor',[.7 .7 .7]);
hold on
patch([P0(1) P1(1) P2(1)], [P0(2) P1(2) P2(2)], [P0(3) P1(3) P2(3)], [1 1 1], 'edgecolor','r')
plot3([0 normalPlot(1)], [0 normalPlot(2)], [middleZ normalPlot(3)], 'r-')
plot3([normalPlot(1)], [normalPlot(2)], [normalPlot(3)], 'r^', 'markerfacecolor','r')
plot3([0 0], [0 0], [0 (FOVsize/2)], 'b-')
plot3([0], [0], [0], 'bv', 'markerfacecolor','b')
plot3(P0(1), P0(2), P0(3), 'ro', 'markerfacecolor','r')
plot3(P1(1), P1(2), P1(3), 'ro', 'markerfacecolor','r')
plot3(P2(1), P2(2), P2(3), 'ro', 'markerfacecolor','r')
text(0,0,FOVsize/2, 'objective ', 'HorizontalAlignment','Right', 'Color','b')
text(normalPlot(1), normalPlot(2), normalPlot(3), '  normal', 'HorizontalAlignment','left', 'Color','r')

axis equal
axis square
axis vis3d
grid off
axis off
xlabel('X (um)')
ylabel('Y (um)')
zlabel('Z (um)')

title([num2str(angleFromStraight, '%.1f') '° from vertical (X: ' num2str(angleX, '%.1f') '°, Y: ' num2str(angleY, '%.1f') '°)'])

% manually plot the axis lines and XY labels (for nicer animation)
xlims = xlim;
ylims = ylim;
zlims = zlim;
plot3(xlims([1 1]),ylim,zlims([1 1]), 'k')
plot3(xlim,ylims([1 1]),zlims([1 1]), 'k')
text(mean(xlim),ylims(1)-50,zlims(1)-50, 'X')
text(xlims(1)-50,mean(ylim),zlims(1)-50, 'Y')

%% animated rotating plot of plane and vectors
% OptionZ.FrameRate=30;
% OptionZ.Duration=10;
% OptionZ.Periodic=true;
% CaptureFigVid([0,30; 180,30; 360,30], 'RotatingPlaneNormal',OptionZ)

%% Plot the objective model
figure

subplot(1,4,1)
title('Objective')
plot3([O0(1) O1(1) ],[O0(2) O1(2) ],[O0(3) O1(3) ], '.-', 'color','k', 'LineWidth',2)
hold on
plot3([O1(1) O2(1) ],[O1(2) O2(2) ],[O1(3) O2(3) ], '.-', 'color','b', 'LineWidth',2)
plot3(O2(1), O2(2), O2(3), 'bv', 'markerfacecolor','b')
plot3(0,0,0, 'r.', 'MarkerSize',25)
plot3(1,0,0, 'r.', 'MarkerSize',25)
axis equal
axis vis3d
grid off
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
xticks([])
yticks([])
zticks([])
ax = gca;
ax.XColor = [.7 .7 .7];
ax.YColor = [.7 .7 .7];
ax.ZColor = [.7 .7 .7];



subplot(1,4,2)
title('Pitch')
plot3(1,0,0, 'r.', 'MarkerSize',25)
hold on
plot3(0,0,0, 'r.', 'MarkerSize',25)
for j = 1:30:numel(pitch_range)
    plot3([O0(1) allO1(1,1)], [O0(2) allO1(1,2)], [O0(3) allO1(1,3)], '-', 'color','k', 'LineWidth',2)
    plot3([allO1(1,1) allO2(1,j,1)], [allO1(1,2) allO2(1,j,2)], [allO1(1,3) allO2(1,j,3)], '-', 'color','b', 'LineWidth',1)
end
axis equal
axis vis3d
grid off
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
xticks([])
yticks([])
zticks([])
ax = gca;
ax.XColor = [.7 .7 .7];
ax.YColor = [.7 .7 .7];
ax.ZColor = [.7 .7 .7];



subplot(1,4,3)
title('Yaw')
plot3(0,0,0, 'r.', 'MarkerSize',25)
hold on
for i = 1:30:numel(yaw_range)
    plot3([allO1(i,1)], [allO1(i,2) ], [allO1(i,3) ], 'r.', 'MarkerSize',25)
    
    plot3([O0(1) allO1(i,1)], [O0(2) allO1(i,2)], [0 0], '-', 'color','k', 'LineWidth',1)
    plot3([allO1(i,1) allO1(i,1)], [allO1(i,2) allO1(i,2)], [0 -1], '-', 'color','b', 'LineWidth',1)
end
axis equal
axis vis3d
grid off
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
xticks([])
yticks([])
zticks([])
ax = gca;
ax.XColor = [.7 .7 .7];
ax.YColor = [.7 .7 .7];
ax.ZColor = [.7 .7 .7];



subplot(1,4,4)
title('Pitch and Yaw')
cmap = hsv(numel(yaw_range));
for i = 1:30:numel(yaw_range)
    for j = 1:30:numel(pitch_range)
        plot3([O0(1) allO1(i,1)], [O0(2) allO1(i,2)], [O0(3) allO1(i,3)], '.-', 'color','k')
        plot3([allO1(i,1) allO2(i,j,1)], [allO1(i,2) allO2(i,j,2)], [allO1(i,3) allO2(i,j,3)], '.-', 'color',cmap(i,:))
    end
    hold on
end

axis equal
axis vis3d
grid off
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2 2])
ylim([-2 2])
zlim([-2 2])
xticks([])
yticks([])
zticks([])
ax = gca;
ax.XColor = [.7 .7 .7];
ax.YColor = [.7 .7 .7];
ax.ZColor = [.7 .7 .7];


%% plot the final direction vector of objective and the plane normal

% objective plot (shifted on top of normal vector)
plotO0 = O0;
plotO1 = allO1(yaw_idx,:);
plotO2 = squeeze(allO2(yaw_idx,pitch_idx,:))';
offset = plotO2 - normal;
plotO0 = plotO0 - offset;
plotO1 = plotO1 - offset;
plotO2 = plotO2 - offset;

figure
plot3([plotO0(1) plotO1(1)], [plotO0(2) plotO1(2)], [plotO0(3) plotO1(3)], '.-', 'markersize',25, 'color','k', 'LineWidth',2)
hold on
plot3([plotO1(1) plotO2(1)], [plotO1(2) plotO2(2)], [plotO1(3) plotO2(3)], '-', 'color','b', 'LineWidth',2)
plot3([0 normal(1)], [0 normal(2)], [0 normal(3)], 'r-')
plot3(mean([plotO1(1) plotO2(1)]), mean([plotO1(2) plotO2(2)]), mean([plotO1(3) plotO2(3)]), 'bv', 'markerfacecolor','b')
plot3(mean([0 normal(1)]), mean([0 normal(2)]), mean([0 normal(3)]), 'r^', 'markerfacecolor','r')

[X,Y] = meshgrid([-1 0 1],[-1 0 1]);
Z = (A*X + B*Y + D) / -C;    % a(x) + b(y) + c(z) + d = 0
Z = Z - mean(Z(:));
surf(X,Y,Z, 'facecolor','w', 'edgecolor',[.7 .7 .7])

axis equal
axis vis3d
grid off
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-2 2])
ylim([-2 2])
% zlim([-2 2])
xticks([])
yticks([])
zticks([])
ax = gca;
ax.XColor = [.7 .7 .7];
ax.YColor = [.7 .7 .7];
ax.ZColor = [.7 .7 .7];


title(['Pitch: ' num2str(pitch) '°, Yaw: ' num2str(yaw), '° (Error: ' num2str(errorFinal, '%.1f') '°)'])
