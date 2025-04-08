%%% Define target structures %%%
% Define E -> entry zone as the vertices of the example circle %%

% Define the circle parameters
% 3 cm entry
r_circle = 3; % radius of the circle
theta_circle = linspace(0, 2*pi, 100); % parameter theta for the circle

% Parametric equations for the circle
x_circle = r_circle * cos(theta_circle);
y_circle = r_circle * sin(theta_circle);
z_circle = -6 * ones(size(theta_circle)); % z = -6

entry_zone = [x_circle.', y_circle.', z_circle.'];

% Find convex hull of E
convex_entry = convhull(entry_zone(:, 1:2));
convexE = [entry_zone(convex_entry, 1), entry_zone(convex_entry, 2), entry_zone(convex_entry,3)];

% Calculate T -> the target zone using the contour points and the vertices of the entry zone 

% Define the sphere parameters where the sphere is the target
% 3 cm target
r_sphere = 1.5; % radius of the sphere
[x_sphere, y_sphere, z_sphere] = sphere(50);
x_sphere = x_sphere * r_sphere;
y_sphere = y_sphere * r_sphere;
z_sphere = z_sphere * r_sphere;

% Add 5 mm margin to sphere to account for dose fall off
r_sphere_margin = (r_sphere + .5);
[x_sphere_margin, y_sphere_margin, z_sphere_margin] = sphere(50);
x_sphere_margin = x_sphere_margin .* r_sphere_margin;
y_sphere_margin = y_sphere_margin .* r_sphere_margin;
z_sphere_margin = z_sphere_margin .* r_sphere_margin;
figure
subplot(121)
    mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);
    contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, 20, 'LineWidth', 2);
    axis tight equal
subplot(122)
    mesh(x_sphere, y_sphere, z_sphere);
    contour3(x_sphere, y_sphere, z_sphere, 20, 'LineWidth', 2);
    axis tight equal


% Create the 3D plot
figure;
hold on;

% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the sphere
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot contours based on the Z-coordinates
contour_levels = 15; % number of contour levels
contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, contour_levels, 'LineWidth', 2);  % 20 contour levels


% Extract the coordinates of the contours
% Create contours on the sphere
[C, h] = contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, [-3.5, -3, -2.5, -2, -1.5, -1, -.5, 0, .5, 1, 1.5, 2, 2.5, 3, 3.5], 'k');
 
contour_data = struct(); % Structure to hold contour data

k = 1; % Index to keep track of contours

i = 1;
while i < size(C, 2)
    level = C(1, i); % Contour level
    num_points = C(2, i); % Number of points at this contour level
    i = i + 1; % Move to the next column for coordinates
    x_contour = C(1, i:i+num_points-1);
    y_contour = C(2, i:i+num_points-1);
    z_contour = level * ones(1, num_points); % The z coordinate is constant for a contour

    % Store the contour data
    contour_data(k).level = level;
    contour_data(k).x = x_contour;
    contour_data(k).y = y_contour;
    contour_data(k).z = z_contour;
    k = k + 1;
    
    i = i + num_points; % Move to the next contour level
end
contour_data_final = struct();
contour_data_final(1).level = -3.5;
contour_data_final(1).x = 0 * ones(1, 51);
contour_data_final(1).y = 0 * ones(1, 51);
contour_data_final(1).z = -3.5 * ones(1, 51);

for i = 1:length(contour_data)
    contour_data_final(i+1).level = contour_data(i).level;
    contour_data_final(i+1).x = contour_data(i).x;
    contour_data_final(i+1).y = contour_data(i).y;
    contour_data_final(i+1).z = contour_data(i).z;
end

z_target = 4.5;

% Compute T
% Contour points come from contour slice closest to entry zone 
% modify T to be more anatomically relevent
z_factor = ((z_target - entry_zone(1, 3)) / (contour_data(1).z(1) - entry_zone(1, 3)));
target_zone = [];
for i = 1:length(entry_zone(:,1))
    for j = 1:length(contour_data(1).x)
        t = [entry_zone(i, 1)/2, entry_zone(i, 2)/2, entry_zone(i, 3)] + z_factor * ([contour_data(1).x(j)/2, contour_data(1).y(j)/2, contour_data(1).z(j)] - [entry_zone(i, 1)/2, entry_zone(i, 2)/2, entry_zone(i, 3)]);
        target_zone = [target_zone; t];
    end
end


% Find convex hull of T
convex_target = convhull(target_zone(:, 1:2));
convexT = [target_zone(convex_target, 1), target_zone(convex_target, 2), target_zone(convex_target, 3)];

% Plot T
plot3(target_zone(:, 1), target_zone(:, 2), target_zone(:, 3),'b', 'LineWidth', 3);

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Set the axis properties
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('3D Plot of Entry zone, Target zone, and Target');
grid on;

%axis vis3d;
%view(3); % Sets the default 3D view

% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);

hold off;

%%% finish initializing target structure %%%
%%
x_range = -3:3;
y_range = -3:3;
z_range = -3:3;

% Define bounding box corners
X = [-3 3 3 -3 -3 3 3 -3];
Y = [-3 -3 3 3 -3 -3 3 3];
Z = [-3 -3 -3 -3 3 3 3 3];

% Define bounding box edges
edges = [1 2; 2 3; 3 4; 4 1; ...  % Bottom face edges
         5 6; 6 7; 7 8; 8 5; ...  % Top face edges
         1 5; 2 6; 3 7; 4 8];     % Vertical edges

% Create the 3D plot
figure;
hold on;

% Plot the sphere
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot bounding box edges
for i = 1:size(edges, 1)
    plot3(X(edges(i, :)), Y(edges(i, :)), Z(edges(i, :)), 'k', 'LineWidth', 2);
end

axis equal;
view(3);
grid on;

xlim([-5 5]);
ylim([-5 5]);
zlim([-5 5]);

hold off;


%%
%%% Generate needles through the target %%%

dtE = delaunayTriangulation(convexE(:,1), convexE(:,2));
dtT = delaunayTriangulation(convexT(:,1), convexT(:,2));

trianglesE = dtE.ConnectivityList;
trianglesT = dtT.ConnectivityList;

lines = struct('line', []);

for i=1:10000
    % Find random triangle 
    [random_TriangleE, random_idxE] = randTri(dtE.ConnectivityList, convexE);
    [random_TriangleT, random_idxT] = randTri(dtT.ConnectivityList, convexT);

    % Unifrmly sample to get points 
    pointE = uniformSample(random_TriangleE);
    pointT = uniformSample(random_TriangleT);

    line = [pointE; pointT];

    lines(i).line = line;
end

figure;
hold on;

% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot contours based on the Z-coordinates
contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, 20, 'LineWidth', 2);  % 20 contour levels

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the line connecting the points
for i=1:length(lines)
    plot3(lines(i).line(:,1), lines(i).line(:,2), lines(i).line(:,3), 'k-');
end

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Target Convex, Entry Convex, ROI, and Target and Entry Meshes, with 100 needles');

% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);
hold off;

%%% finish needle generation %%%
%%
dtE = delaunayTriangulation(convexE(:,1), convexE(:,2));
dtT = delaunayTriangulation(convexT(:,1), convexT(:,2));

trianglesE = dtE.ConnectivityList;
trianglesT = dtT.ConnectivityList;


[random_TriangleE, random_idxE] = randTri(dtE.ConnectivityList, convexE);
[random_TriangleT, random_idxT] = randTri(dtT.ConnectivityList, convexT);

highlightedTriangleE = trianglesE(random_idxE, :);
highlightedTriangleT = trianglesT(random_idxT, :);

pointE = uniformSample(random_TriangleE);
pointT = uniformSample(random_TriangleT);

line = [pointE; pointT];

figure;
hold on;

% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');
% Plot randomly selected traingle 
patch('Vertices', convexE, 'Faces', highlightedTriangleE, 'FaceColor', 'red');
% Plot uniformly sampled point 
plot3(pointE(1,1), pointE(1,2), pointE(1,3), 'o', 'Color', 'g', 'LineWidth', 5);

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');
% Plot randomly selected triangle
patch('Vertices', convexT, 'Faces', highlightedTriangleT, 'FaceColor', 'red');
% Plot uniformly sampled point
plot3(pointT(1,1), pointT(1,2), pointT(1,3), 'o', 'Color', 'g', 'LineWidth', 5);

% Plot the line connecting the points
plot3(line(:,1), line(:,2), line(:,3), 'k-', 'LineWidth', 2); 


axis equal;
% axis vis3d;
% view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Sample Needle Generation');
hold off;


%%
%%% Selectively color needles based on parity and remove if no intersection %%%
sphere_points = [x_sphere_margin(1,1), y_sphere_margin(1,1), z_sphere_margin(1,1);
    x_sphere_margin(2,2), y_sphere_margin(2,2), z_sphere_margin(2,2);
    x_sphere_margin(3,3), y_sphere_margin(3,3), z_sphere_margin(3,3);
    x_sphere_margin(4,4), y_sphere_margin(3,3), z_sphere_margin(4,4)];

% unit sphere with radius 1
sphere_margin_center = [0, 0, 0];
matching_indicies = [];
needle_tracks = struct('track',[]);

for i=1:length(lines)
    if checkIntersect(lines(i).line, sphere_margin_center, r_sphere_margin)
        if (lines(i).line(1,1) < 0 && lines(i).line(2,1) < 0) || (lines(i).line(1,1) > 0 && lines(i).line(2,1) > 0)
            matching_indicies = [matching_indicies, i];
        end
    end
end


for i=1:length(matching_indicies)
    needle_tracks(i).track = lines(matching_indicies(i)).line;
end

for i = 1:length(needle_tracks)
    needle_tracks(i).idx = i;
end


for i = 1:length(needle_tracks)
    v1 = needle_tracks(i).track(1, :);
    v2 = needle_tracks(i).track(2, :);
    u = v2 - v1;
    u_mag = norm(u);
    unit = u / u_mag;
    needle_tracks(i).unit = unit;
end

% select needles based on angle with z axis
z_unit = [0; 0; 1];
angles = [];
for i = 1:length(needle_tracks)
    ang = acosd(dot(needle_tracks(i).unit, z_unit));
    angles = [angles, ang];
end

ang_threshhold = 10; % degrees
needle_tracks = needle_tracks(angles <= ang_threshhold);


figure;
hold on;
% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot contours based on the Z-coordinates
contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, 20, 'LineWidth', 2);  % 20 contour levels

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the line connecting the points
% Plot the line connecting the points
for i=1:length(needle_tracks)
    if (needle_tracks(i).track(1,1) < 0 && needle_tracks(i).track(2,1) < 0)
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'r-');
    else
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'b-');
    end
end

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Target Convex, Entry Convex, ROI, and Target and Entry Meshes');
% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);
hold off;

%%% finish needle generation %%%
%%
%%% define needle dwells 
% Because not coming from image data, scale is arbitrary -> choose to define such that 1 = 1 cm
needle_offset = .7;

% Define first dwell position .7 cm from tip 
for i = 1:length(needle_tracks)

    % Get starting and ending points 
    start = needle_tracks(i).track(2,:);
    ending = needle_tracks(i).track(1,:);

        % Define unit vector
    unit = start - ending;
    unit_norm = unit / norm(unit);
    scale_unit_norm = unit_norm * needle_offset;

    % Compute first dwell position
    first_dwell = start - scale_unit_norm;

    % Store in needles
    needle_tracks(i).first_dwell = first_dwell;
end

% define dwells all the way to entry zone

dwells = [];
% Recall all distances in mm
for i = 1:length(needle_tracks)
    last_dwell = needle_tracks(i).track(1, :);
    dwells = generatePointsAlongLine(needle_tracks(i).first_dwell, last_dwell, .5);
    needle_tracks(i).dwells = dwells;
end

figure;
hold on;
% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot contours based on the Z-coordinates
contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, 20, 'LineWidth', 2);  % 20 contour levels

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the line connecting the points
% Plot the line connecting the points
for i=1:length(needle_tracks)
    if (needle_tracks(i).track(1,1) < 0 && needle_tracks(i).track(2,1) < 0)
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'r-');
        plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
        for j = 1:length(needle_tracks(i).dwells)
            plot3(needle_tracks(i).dwells(j, 1), needle_tracks(i).dwells(j, 2), needle_tracks(i).dwells(j, 3), 'k.', 'MarkerSize', 10);
        end
    elseif (needle_tracks(i).track(1,1) > 0 && needle_tracks(i).track(2,1) > 0)
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'b-');
        plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
        for j = 1:length(needle_tracks(i).dwells)
            plot3(needle_tracks(i).dwells(j, 1), needle_tracks(i).dwells(j, 2), needle_tracks(i).dwells(j, 3), 'k.', 'MarkerSize', 10);
        end
    else
        continue;
    end
end

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Target Convex, Entry Convex, ROI, and Target and Entry Meshes');
% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);
hold off;

%%% end definition of dwells %%%
%%
%%% create set of 'active dwells', dwells that are contained within the sphere_margin (target sphere + 5 mm margin) + an additional margin to allow
%%% for needle extension beyond the HR_CTV

for i = 1:length(needle_tracks)
    active_dwells = [];
    for j = 1:length(needle_tracks(i).dwells)
        if isPointInSphere(needle_tracks(i).dwells(j,:), sphere_margin_center, r_sphere_margin + .5)
            active_dwells = [active_dwells; needle_tracks(i).dwells(j,:)];
        end
    end
    needle_tracks(i).active_dwells = active_dwells;
end

figure;
hold on;
% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot contours based on the Z-coordinates
contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, 20, 'LineWidth', 2);  % 20 contour levels

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

for i=1:length(needle_tracks)
    if (needle_tracks(i).track(1,1) < 0 && needle_tracks(i).track(2,1) < 0)
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'r-');
        if isPointInSphere(needle_tracks(i).first_dwell, sphere_margin_center, r_sphere_margin)
            plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
        end
        for j = 1:size(needle_tracks(i).active_dwells)
            plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'k.', 'MarkerSize', 10);
        end
        elseif (needle_tracks(i).track(1,1) > 0 && needle_tracks(i).track(2,1) > 0)
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'b-');
        if isPointInSphere(needle_tracks(i).first_dwell, sphere_margin_center, r_sphere_margin)
            plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
        end
        for j = 1:size(needle_tracks(i).active_dwells)
            plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'k.', 'MarkerSize', 10);
        end
    else
        continue;
    end
end

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Target Convex, Entry Convex, ROI, and Target and Entry Meshes');
% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);
hold off;

%%
% Remove needles within 5mm of eachother
% Assume needles is a cell array where each cell contains an Mx3 matrix
% representing the 3D coordinates of the dwell positions for one needle.

% % Remove any needles with 1 or no active dwells
% for i = 1:length(needle_tracks)
%     num_dwells = size(needle_tracks(i).active_dwells, 1);
%     if (isempty(needle_tracks(i).active_dwells) || num_dwells == 1)
%         needle_tracks(i) = [];
%     end
% end

numNeedles = length(needle_tracks);

% Initialize a logical array to mark selected needles
selected = true(numNeedles, 1);

% Iterate through each needle
for i = 1:numNeedles
    if selected(i) % Only check if the current needle is still selected
        % Get the dwell positions of the current needle
        currentNeedleDwellPositions = needle_tracks(i).active_dwells;
        
        % Iterate over the other needles
        for j = i+1:numNeedles
            if selected(j) % Only check if the other needle is still selected
                % Get the dwell positions of the other needle
                otherNeedleDwellPositions = needle_tracks(j).active_dwells;
                
                % Compute pairwise distances between all dwell positions
                distances = pdist2(currentNeedleDwellPositions, otherNeedleDwellPositions);
                
                % Check if any distance is less than 0.5 cm
                if any(distances(:) < 0.5)
                    % Mark one of the needles as not selected
                    % For this example, we remove needle `j` (can be adjusted based on priority)
                    selected(j) = false;
                end
            end
        end
    end
end

% Filter the needle set
needle_tracks = needle_tracks(selected);

figure;
hold on;
% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot contours based on the Z-coordinates
contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, 20, 'LineWidth', 2);  % 20 contour levels

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

for i=1:length(needle_tracks)
    if (needle_tracks(i).track(1,1) < 0 && needle_tracks(i).track(2,1) < 0)
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'r-');
        if isPointInSphere(needle_tracks(i).first_dwell, sphere_margin_center, r_sphere_margin)
            plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
        end
        for j = 1:size(needle_tracks(i).active_dwells)
            plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'k.', 'MarkerSize', 10);
        end
        elseif (needle_tracks(i).track(1,1) > 0 && needle_tracks(i).track(2,1) > 0)
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'b-');
        if isPointInSphere(needle_tracks(i).first_dwell, sphere_margin_center, r_sphere_margin)
            plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
        end
        for j = 1:size(needle_tracks(i).active_dwells)
            plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'k.', 'MarkerSize', 10);
        end
    else
        continue;
    end
end

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Target Convex, Entry Convex, ROI, and Target and Entry Meshes');
% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);
hold off;

%%
figure;
hold on;


% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin,'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', .5);


% for i=1:length(needle_tracks)
%     if (needle_tracks(i).track(1,1) < 0 && needle_tracks(i).track(2,1) < 0)
%         plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'r-');
%         if isPointInSphere(needle_tracks(i).first_dwell, sphere_margin_center, r_sphere_margin)
%             plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
%         end
%         for j = 1:size(needle_tracks(i).active_dwells)
%             plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'k.', 'MarkerSize', 10);
%         end
%         elseif (needle_tracks(i).track(1,1) > 0 && needle_tracks(i).track(2,1) > 0)
%         plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'b-');
%         if isPointInSphere(needle_tracks(i).first_dwell, sphere_margin_center, r_sphere_margin)
%             plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
%         end
%         for j = 1:size(needle_tracks(i).active_dwells)
%             plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'k.', 'MarkerSize', 10);
%         end
%     else
%         continue;
%     end
% end

idx = randi([1, length(needle_tracks)]);
plot3(needle_tracks(idx).track(:,1), needle_tracks(idx).track(:,2), needle_tracks(idx).track(:,3), 'b-');
plot3(needle_tracks(idx).first_dwell(1,1), needle_tracks(idx).first_dwell(1,2), needle_tracks(idx).first_dwell(1,3), 'g.', 'MarkerSize', 20);
for j = 1:size(needle_tracks(idx).active_dwells, 1)
    plot3(needle_tracks(idx).active_dwells(j, 1), needle_tracks(idx).active_dwells(j, 2), needle_tracks(idx).active_dwells(j, 3), 'c.', 'MarkerSize', 20);
end
for k = 1:size(needle_tracks(idx).dwells)
    if k == 1
        continue;
    end 
    plot3(needle_tracks(idx).dwells(k, 1), needle_tracks(idx).dwells(k, 2), needle_tracks(idx).dwells(k, 3), 'r.', 'MarkerSize', 15);
end

axis equal;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Target Convex, Entry Convex, ROI, and Target and Entry Meshes');
% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);
hold off;
%% 
% Add the central needle to act as a tandem 
tandem = needle_tracks(1);
tandem.idx = 1;
tandem.track = [];
tandem.unit = [];
tandem.first_dwell = [];
tandem.dwells = [];

track = [0,0, convexE(1, 3); 0, 0, convexT(1, 3)];

% track and index
tandem.track = track;
tandem.idx = 1;

% unit vector
v1 = tandem.track(1, :);
v2 = tandem(end).track(2, :);
u = v2 - v1;
u_mag = norm(u);
unit = u / u_mag;
tandem.unit = unit;

% first dwell
start = tandem(end).track(1,:);
ending = tandem(end).track(2,:);

scale_unit_norm = tandem(end).unit * 0.7;

first_dwell = ending + scale_unit_norm;

tandem.first_dwell = first_dwell;

% all dwells 
dwells = [];
last_dwell = tandem.track(1, :);
dwells = generatePointsAlongLine(tandem.first_dwell, last_dwell, .5);
tandem.dwells = dwells;

% active dwells
active_dwells = [];
for j = 1:length(tandem.dwells)
    if isPointInSphere(tandem.dwells(j,:), sphere_margin_center, r_sphere_margin + .5)
        active_dwells = [active_dwells; tandem.dwells(j,:)];
    end
end
tandem.active_dwells = active_dwells;

% Check if any needle has idx == 1
has_idx = any([needle_tracks(:).idx] == 1);

% Check if any needle has unit == [0, 0, 1]
has_unit = any(arrayfun(@(s) isequal(s.unit, [0, 0, 1]), needle_tracks));

% If both conditions are met, stop execution
if has_idx && has_unit
    disp('Tandem already inserted into needle tracks, quitting to not duplicate...');
    return;
end

needle_tracks(end+1) = tandem;

% rearrange so 'tandem' is at the front
needle_tracks_temp = needle_tracks;

needle_tracks(1) = needle_tracks_temp(end);
needle_tracks(2:end) = needle_tracks_temp(1:end-1);

figure;
hold on;
% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot contours based on the Z-coordinates
contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, 20, 'LineWidth', 2);  % 20 contour levels

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

for i=1:length(needle_tracks)
    if i == 1
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'r-', 'LineWidth', 2);
        for j = 1:size(needle_tracks(i).active_dwells, 1)
            plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'm.', 'MarkerSize', 10);
        end
    else
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'k-', 'LineWidth', 2);
        for j = 1:size(needle_tracks(i).active_dwells, 1)
            plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'm.', 'MarkerSize', 10);
        end
    end
end

axis equal;

view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Target Convex, Entry Convex, ROI, and Target and Entry Meshes');
% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);
hold off;

%%
% Calculate unit dose rate contribution per needle

% Initialize global variables
global F_table gLr_table G_L_r_0 L S_k lambda;

consensus_files = dir('Consensus Data/*.mat');
for i = 1:length(consensus_files)
    load("Consensus Data/" + consensus_files(i).name + ''); % Load each file
end

% Initialize tables and constants for dose calculation 
F_table = Frtheta;
gLr_table = gLr;
L = .35;
G_L_r_0 = G_L(1, L, 90);

% Nominal 10 Ci source strength + source specific dose rate constant 
S_k = 40700; 
lambda = Lambda;

% Create Dose Calculation Grid 
% 1mm pixel spacing in x and y, 2.5 slice thickness in z
% xRange = -3:.1:3;
% yRange = -3:.1:3;
% zRange = -3:.25:3;

% Mid range 
xRange = -3.5:.25:3.5;
yRange = -3.5:.25:3.5;
zRange = -3.5:.25:3.5;

% Old range -> more coarse, less accurate 
 % xRange = -3.5:.5:3.5;
 % yRange = -3.5:.5:3.5;
 % zRange = -3.5:.5:3.5;

[xGrid, yGrid, zGrid] = meshgrid(xRange, yRange, zRange);

% Create target mask
target_center = sphere_margin_center;
target_r = r_sphere_margin;

d_squared2center = (xGrid - target_center(1)).^2 + (yGrid - target_center(2)).^2 + ...
                   (zGrid - target_center(3)).^2;
target_mask = d_squared2center <= target_r^2;

% Get grid size and total number of voxels
gridSize = size(xGrid);
num_voxels = prod(gridSize);

% Initialize the total dose distribution
dose_distribution = zeros(gridSize);

% Loop over each needle
for i = 1:length(needle_tracks)
    % Get the number of active dwell positions for this needle
    num_dwells = size(needle_tracks(i).active_dwells, 1);

    % Initialize the unit dose matrix for this needle (V_total x J)
    needle_unit_dose = zeros(num_voxels, num_dwells);

    % % Get the dwell times for this needle
    % dwell_end_idx = dwell_start_idx + num_dwells - 1;
    % dwells = dwell_times(dwell_start_idx:dwell_end_idx);
    % dwell_start_idx = dwell_end_idx + 1;

    % Loop over each dwell position for this needle
    for j = 1:num_dwells
        % Get source position
        sourcePosition = needle_tracks(i).active_dwells(j, :);

        % Loop over each voxel in the entire volume
        for xi = 1:gridSize(1)
            for yi = 1:gridSize(2)
                for zi = 1:gridSize(3)
                    % Linear index of the current voxel
                    voxel_idx = sub2ind(gridSize, xi, yi, zi);

                    % Get the grid point position
                    gridPoint = [xGrid(xi, yi, zi), yGrid(xi, yi, zi), zGrid(xi, yi, zi)];

                    % Compute distance and angle
                    r_vec = gridPoint - sourcePosition;
                    r = norm(r_vec);

                    % Skip zero-distance to avoid division by zero
                    if r > 0
                        % Calculate theta using dot product
                        dot = sum(conj(needle_tracks(i).unit) .* r_vec);
                        theta = acosd(abs(dot) / r);

                        % Compute the dose contribution for this voxel and dwell
                        G_Lr_s = abs(G_L(r, L, theta) / G_L_r_0);
                        g_Lr_s = g_Lr(gLr_table, r);
                        F_rt_s = F_rt(F_table, r, theta);

                        doseRate = S_k * lambda * G_Lr_s * g_Lr_s * F_rt_s;

                        % Store in the unit dose matrix
                        needle_unit_dose(voxel_idx, j) = doseRate;
                    end
                end
            end
        end
    end

    % Store the unit dose matrix in the needle structure
    needle_tracks(i).unit_dose = needle_unit_dose;

    % Display progress
    disp(['Needle ', num2str(i), '/', num2str(length(needle_tracks)), ' (', num2str(i / length(needle_tracks) * 100), '%) complete']);
end


for i = 1:length(needle_tracks)
    % Logical indices of the target voxels (flattened into a vector)
    target_voxel_indices = find(target_mask);

    % Extract only the rows (voxels) corresponding to the target
    needle_tracks(i).target_dose = needle_tracks(i).unit_dose(target_voxel_indices, :);
end

%% Test Dose Calcs to make sure the same 
% num_dc_dwells = 0;
% for i = 1:length(needle_tracks)
%     for j = 1 : length(needle_tracks(i).active_dwells)
%         num_dc_dwells = num_dc_dwells + 1;
%     end
% end
% dc_dwell_times = 100*ones(num_dc_dwells, 1);
% 
% dose_calc_dist = dose_calc(needle_tracks, dc_dwell_times);
% 
% % Original dose calc loops not in function
% dose_loop = zeros(size(xGrid));
% 
% % Calculate dose rate contribution for each needle to the whole grid 
% for i = 1:length(needle_tracks) % random chosen just to generate grid
%     for j = 1:length(needle_tracks(i).active_dwells)
%         % Get source position
%         sourcePosition = needle_tracks(i).active_dwells(j, :);
% 
%         for xi = 1:gridSize(1)
%             for yi = 1:gridSize(2)
%                 for zi = 1:gridSize(3)
%                     % Grid point position
%                     gridPoint = [xGrid(xi, yi, zi), yGrid(xi, yi, zi), zGrid(xi, yi, zi)];
% 
%                     % Compute distance and angle
%                     r_vec = gridPoint - sourcePosition;
%                     r = norm(r_vec);
% 
%                     % Calculate theta using dot product
%                     if r > 0
%                         dot = sum(conj(needle_tracks(i).unit) .* r_vec);
%                         theta = acosd(abs(dot) / r);
%                     else
%                         theta = 0; % Avoid NaN for zero distance
%                     end
% 
%                     % Dose calculation (if distance > 0)
%                     if r > 0
%                         G_Lr_s = abs(G_L(r, L, theta) / G_L_r_0);
%                         g_Lr_s = g_Lr(gLr_table, r);
%                         F_rt_s = F_rt(F_table, r, theta);
% 
%                         doseRate = S_k * lambda * G_Lr_s * g_Lr_s * F_rt_s;
% 
%                         dose_loop(xi, yi, zi) = dose_loop(xi, yi, zi) + (doseRate * 100/3600) ; % rand number between 0 and 2 seconds
%                     end
%                 end
%             end
%         end
%     end
% 
%     disp(['Needle ', num2str(i), '/', num2str(length(needle_tracks)), '(', num2str(i/length(needle_tracks) * 100), '%) complete']);
% end
% 
% disp('Dose calculation finished')

%%

Rx = 680; % cGy
dose_min = Rx; % 600 cGy
dose_max = 1.5 * Rx; % 1200 cGy
target_dose_range = [dose_min, dose_max];

num_dwells = 0;
for i = 1 : length(needle_tracks)
    for j = 1 : length(needle_tracks(i).active_dwells)
        num_dwells = num_dwells + 1;
    end
end 

% Dose Calculation Grid and Ranges defined above during needle unit dose calculation  

% Create target mask
target_center = sphere_margin_center;
target_r = r_sphere_margin;

d_squared2center = (xGrid - target_center(1)).^2 + (yGrid - target_center(2)).^2 + ...
                   (zGrid - target_center(3)).^2;
target_mask = d_squared2center <= target_r^2;

for i = 1:length(needle_tracks)
    num_dwells = length(needle_tracks(i).active_dwells);
    if num_dwells == 1
        needle_tracks(i) = [];
    end
end


%[selected_milp, times_milp] = milp_needle_optimization(needle_tracks, target_dose_range, 15, target_mask);
[selected_needles_opt, dwell_times_opt] = greedy_needle_selection(needle_tracks, target_dose_range, 15, target_mask, xGrid, yGrid, zGrid);
disp("Optimization Done, starting final dose calc...")
dose_needles_opt = dose_calc(selected_needles_opt, dwell_times_opt, xGrid);
disp("Final dose calc done, continue to plotting");
%%
save('Sphere_vol_results/selected_needles_opt_5.mat', 'selected_needles_opt');
save('Sphere_vol_results/dwell_times_opt_5.mat', 'dwell_times_opt');
save('Sphere_vol_results/dose_needles_opt_5.mat', 'dose_needles_opt');

%%
% Extract the coordinates of the contours
% Create contours on the sphere
[C, h] = contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, zRange, 'k');

contour_data = struct(); % Structure to hold contour data

k = 1; % Index to keep track of contours

i = 1;
while i < size(C, 2)
    level = C(1, i); % Contour level
    num_points = C(2, i); % Number of points at this contour level
    i = i + 1; % Move to the next column for coordinates
    x_contour = C(1, i:i+num_points-1);
    y_contour = C(2, i:i+num_points-1);
    z_contour = level * ones(1, num_points); % The z coordinate is constant for a contour

    % Store the contour data
    contour_data(k).level = level;
    contour_data(k).x = x_contour;
    contour_data(k).y = y_contour;
    contour_data(k).z = z_contour;
    k = k + 1;
    
    i = i + num_points; % Move to the next contour level
end
contour_data_final = struct();
contour_data_final(1).level = -2;
contour_data_final(1).x = 0 * ones(1, 51);
contour_data_final(1).y = 0 * ones(1, 51);
contour_data_final(1).z = -2 * ones(1, 51);

for i = 1:length(contour_data)
    contour_data_final(i+1).level = contour_data(i).level;
    contour_data_final(i+1).x = contour_data(i).x;
    contour_data_final(i+1).y = contour_data(i).y;
    contour_data_final(i+1).z = contour_data(i).z;
end

%%
dwell_times_opt_redo = optimize_dwell_times(selected_needles_opt, [588, 750], target_mask, xGrid);
dose_needles_opt_redo = dose_calc(selected_needles_opt, dwell_times_opt_redo, xGrid);

%%
target_dose = dose_needles_opt(target_mask == 1);
target_dose = sort(target_dose, 'ascend');

% Compute volume corresponding to each dose level
numVoxels = numel(target_dose);
totalVolume = numVoxels * ((.25*.25*.25)); % Total target volume

volumeRemaining = (numVoxels:-1:1) * (.25*.25*.25); % Absolute volume
percentvolumeRemaining = volumeRemaining / totalVolume;

D150 = 1.5 * Rx;
D200 = 2 * Rx;

numVoxels_V150 = sum(target_dose(:) >= D150);
numVoxels_V200 = sum(target_dose(:) >= D200);

% Compute volumes
V150_abs = numVoxels_V150 * (.25*.25*.25); % Absolute volume (cc)
V150_rel = (V150_abs / totalVolume) * 100; % Relative volume (%)

V200_abs = numVoxels_V200 * (.25*.25*.25); % Absolute volume (cc)
V200_rel = (V200_abs / totalVolume) * 100; % Relative volume (%)

% Get D90 (dose at 90% volume)
index = round(0.1 * numVoxels);  % 10% of the smallest doses
D90 = target_dose(index);

 % Plot DVH
figure;
plot(target_dose, percentvolumeRemaining, 'b-', 'LineWidth', 2);
xlabel('Dose (Gy)');
ylabel('Relative Volume (%)');
title('Dose-Volume Histogram (DVH)');
grid on;
xlim([min(target_dose), 2000]);
ylim([0, 1]);

disp(['D90: ', num2str(D90)]);
disp(['V150: ', num2str(V150_abs), ' (absolute), ', num2str(V150_rel), ' (relative)']);
disp(['V200: ', num2str(V200_abs), ' (absolute), ', num2str(V200_rel), ' (relative)']);

%%
% Mid range 
xRange = -3.5:.25:3.5;
yRange = -3.5:.25:3.5;
zRange = -3.5:.25:3.5;

[xGrid, yGrid, zGrid] = meshgrid(xRange, yRange, zRange);

loaded_data = load("Sphere_vol_results/dose_needles_opt_3.mat");
dose_needles_opt = loaded_data.dose_needles_opt;

loaded_data = load("Sphere_vol_results/selected_needles_opt_3.mat");
selected_needles_opt = loaded_data.selected_needles_opt;

%%
% Extract isosurface (3D isodose contour)
figure;
hold on;

isosurf = isosurface(xGrid, yGrid, zGrid, dose_needles_opt, 600);
p = patch(isosurf, 'FaceColor', 'green', 'EdgeColor', 'none');

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin, 'FaceAlpha', .5);

axis equal;
hold off;

%% Plotting

% Initialize variables to store slice data
%%% Isodose levels %%%
isodose_levels = [1260.6, 600, 588, 500, 400, 300, 200, 100, 50, 20];
colors = [
    0.3871, 0.9056, 0.7088;
    0.5888, 0.1904, 0.1904;
    0.1023, 0.8296, 0.5910;
    0.6873, 0.0685, 0.9229;
    0.7992, 0.2411, 0.2136;
    0.2151, 0.3238, 0.5223;
    0.4388, 0.3121, 0.6007;
    0.1755, 0.3129, 0.3797;
    0.4605, 0.7567, 0.2297;
    0.5128, 0.5832, 0.0918];
slice_figures = struct();
c = 0;

% Store dose grid slices with isodose lines
for s = 1:length(zRange)
    % Create a new figure for each slice
    fig = figure('Visible', 'off'); 
    hold on;

    if -r_sphere_margin <= zRange(s) && zRange(s) <= r_sphere_margin
        c = c + 1;

        % Plot target contour
        plot3(contour_data_final(c).x, contour_data_final(c).y, contour_data_final(c).z, 'w-', 'LineWidth', 2);
    end

    % Plot dose slice
    slice(xGrid, yGrid, zGrid, dose_needles_opt, [], [], zRange(s));
    shading interp;
    colormap(jet);
    colorbar;

    % Track which isodose levels have already been added to the legend
    added_to_legend = false(1, length(isodose_levels));
    visible_isodose_labels = {};
    legend_handles = [];

    for i = 1:length(isodose_levels)
        % Extract the current dose slice
        dose_slice = squeeze(dose_needles_opt(:, :, s));

        % Compute the contour data for the current isodose level
        temp_fig = figure('Visible', 'off');
        [C, ~] = contour(xRange, yRange, dose_slice, [isodose_levels(i) isodose_levels(i)]);
        close(temp_fig);

        % Check if the contour data is valid
        if isempty(C)
            continue; % Skip this isodose level if no contour is found
        end

        % Parse the contour matrix
        k = 1; % Starting index in C
        while k < size(C, 2)
            level = C(1, k); % Isodose level (should match isodose_levels(i))
            num_points = C(2, k); % Number of points in the segment

            % Extract the segment points
            x_contour = C(1, k+1:k+num_points);
            y_contour = C(2, k+1:k+num_points);

            % Skip the contour if it lies entirely outside the current plane
            if all(x_contour < min(xRange) | x_contour > max(xRange)) || ...
               all(y_contour < min(yRange) | y_contour > max(yRange))
                k = k + num_points + 1; % Move to the next contour segment
                continue;
            end

            % Plot the contour at the correct Z level
            z_contour = zRange(s) * ones(size(x_contour)); % Set Z to the current slice level
            handle = plot3(x_contour, y_contour, z_contour, 'Color', colors(i, :), 'LineWidth', 2);

            % Add to the legend only if it hasn't been added yet
            if ~added_to_legend(i)
                visible_isodose_labels{end+1} = [num2str(isodose_levels(i)),' cGy'];
                legend_handles = [legend_handles, handle];
                added_to_legend(i) = true; % Mark this level as added
            end

            % Move to the next segment
            k = k + num_points + 1;
        end
    end

    % Plot needles and dwell positions
    for i = 1:length(selected_needles_opt)
        track = selected_needles_opt(i).track;
        plot3(track(:, 1), track(:, 2), track(:, 3), 'r-', 'LineWidth', 2.5);
        plot3(selected_needles_opt(i).active_dwells(:, 1), selected_needles_opt(i).active_dwells(:, 2), ...
              selected_needles_opt(i).active_dwells(:, 3), 'ro', 'MarkerSize', 6);
    end
    % for j = 1:length(needle_tracks)
    %     track1 = needle_tracks(j).track;
    %     plot3(track1(:, 1), track1(:, 2), track1(:,3), 'k-', 'LineWidth', 1.5);
    % end

    % Customize axis and labels
    hold off;
    xlim([-3, 3]); ylim([-3, 3]); zlim([-3, 3]);
    xlabel('X'); ylabel('Y'); zlabel('Z');

    % Set the legend dynamically
    if ~isempty(legend_handles)
        legend(legend_handles, visible_isodose_labels, 'FontSize', 16, 'Location', 'best');
    end

    view(3);
    rotate3d on;  % Enable 3D rotation

    % Store the figure handle in the struct (instead of the image)
    slice_figures(s).fig_handle = fig;
    slice_figures(s).title = ['Dose Slice at Z = ', num2str(zRange(s))];

    % Pause to allow for figure rendering
    pause(0.1);
end

% Display the first slice initially
current_index = 1;
current_fig = slice_figures(current_index).fig_handle;  % Retrieve the first figure handle
figure(current_fig); % Bring the figure to the front

% Add a slider to navigate the slices
slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(zRange), ...
    'Value', 1, 'SliderStep', [1/(length(zRange)-1), 10/(length(zRange)-1)], ...
    'Position', [20, 20, 300, 20], ...
    'Callback', @(src, event) updateSliceWithIsoDose(src, slice_figures));

% Callback function to update the displayed slice
function updateSliceWithIsoDose(slider, slice_figures)
    % Get the current index from the slider
    index = round(slider.Value);

    % Get the figure handle for the selected slice
    selected_fig = slice_figures(index).fig_handle;

    % Bring the selected slice figure to the front
    figure(selected_fig); 

    % Update the title
    title(slice_figures(index).title);

    % Refresh the display (this allows for 3D rotation to be retained)
    drawnow;
end
%%
% Testing uniform needle distribution

test_tracks = struct();
test_tracks(1).track = [0, 0, -6; 0, 0, 4.5];
test_tracks(2).track = [1, 1, -6; 1, 1, 4.5];
test_tracks(3).track = [0, 1.5, -6; 0, 1.5, 4.5];
test_tracks(4).track = [-1, 1, -6; -1, 1, 4.5];
test_tracks(5).track = [-1, -1, -6; -1, -1, 4.5];
test_tracks(6).track = [0, -1.5, -6; 0, -1.5, 4.5];
test_tracks(7).track = [1, -1, -6; 1, -1, 4.5];
test_tracks(8).track = [-1.5, 0, -6; -1.5, 0, 4.5];
test_tracks(9).track = [1.5, 0, -6; 1.5, 0, 4.5];


for i = 1:length(test_tracks)
    test_tracks(i).idx = i;
end


for i = 1:length(test_tracks)
    v1 = test_tracks(i).track(1, :);
    v2 = test_tracks(i).track(2, :);
    u = v2 - v1;
    u_mag = norm(u);
    unit = u / u_mag;
    test_tracks(i).unit = unit;
end

needle_offset = .7;

% Define first dwell position .7 cm from tip 
for i = 1:length(test_tracks)

    % Get starting and ending points 
    start = test_tracks(i).track(1,:);
    ending = test_tracks(i).track(2,:);

        % Define unit vector
    unit = start - ending;
    unit_norm = unit / norm(unit);
    scale_unit_norm = unit_norm * needle_offset;

    % Compute first dwell position
    first_dwell = ending + scale_unit_norm;

    % Store in needles
    test_tracks(i).first_dwell = first_dwell;
end

% define dwells all the way to entry zone

dwells = [];
% Recall all distances in mm
for i = 1:length(test_tracks)
    last_dwell = test_tracks(i).track(1, :);
    dwells = generatePointsAlongLine(test_tracks(i).first_dwell, last_dwell, .5);
    test_tracks(i).dwells = dwells;
end

for i = 1:length(test_tracks)
    active_dwells = [];
    for j = 1:length(test_tracks(i).dwells)
        if isPointInSphere(test_tracks(i).dwells(j,:), sphere_margin_center, r_sphere_margin + .5)
            active_dwells = [active_dwells; test_tracks(i).dwells(j,:)];
        end
    end
    test_tracks(i).active_dwells = active_dwells;
end

figure;
hold on;
% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(entry_zone(convex_entry,1), entry_zone(convex_entry,2), entry_zone(convex_entry, 3), 'r-', 'LineWidth', 2);

% Plot the target (target includes margin)
mesh(x_sphere_margin, y_sphere_margin, z_sphere_margin);

% Plot contours based on the Z-coordinates
contour3(x_sphere_margin, y_sphere_margin, z_sphere_margin, 20, 'LineWidth', 2);  % 20 contour levels

% Plot the convex hull of T
plot3(target_zone(convex_target,1), target_zone(convex_target,2), target_zone(convex_target, 3), 'r', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

for i=1:length(test_tracks)
    plot3(test_tracks(i).track(:,1), test_tracks(i).track(:,2), test_tracks(i).track(:,3), 'k-');
    for j = 1:size(test_tracks(i).active_dwells, 1)
        plot3(test_tracks(i).active_dwells(j, 1), test_tracks(i).active_dwells(j, 2), test_tracks(i).active_dwells(j, 3), 'r.', 'MarkerSize', 10);
    end

end

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Target Convex, Entry Convex, ROI, and Target and Entry Meshes');
% Zoom out by adjusting axis limits
padding = 5; % Define padding around objects
xlim([-5*r_sphere-padding, 5*r_sphere+padding]);
ylim([-5*r_sphere-padding, 5*r_sphere+padding]);
zlim([-5*r_sphere-padding, 5*r_sphere+padding]);
hold off;
%%
% Test tracks needle unit dose 
% Initialize tables and constants for dose calculation 
F_table = Frtheta;
gLr_table = gLr;
L = .35;
G_L_r_0 = G_L(1, L, 90);

% Nominal 10 Ci source strength + source specific dose rate constant 
S_k = 40700; 
lambda = Lambda;

% % Create Dose Calculation Grid 
% % 1mm pixel spacing in x and y, 2.5 slice thickness in z
% xRange = -7.5:.5:7.5;
% yRange = -7.5:.5:7.5;
% zRange = -7.5:.5:7.5;

[xGrid, yGrid, zGrid] = meshgrid(xRange, yRange, zRange);

% Create target mask
target_center = sphere_margin_center;
target_r = r_sphere_margin;

d_squared2center = (xGrid - target_center(1)).^2 + (yGrid - target_center(2)).^2 + ...
                   (zGrid - target_center(3)).^2;
target_mask = d_squared2center <= target_r^2;

% Get grid size and total number of voxels
gridSize = size(xGrid);
num_voxels = prod(gridSize);

% Initialize the total dose distribution
dose_distribution = zeros(gridSize);

% Loop over each needle
for i = 1:length(test_tracks)
    % Get the number of active dwell positions for this needle
    num_dwells = size(test_tracks(i).active_dwells, 1);

    % Initialize the unit dose matrix for this needle (V_total x J)
    test_unit_dose = zeros(num_voxels, num_dwells);

    % % Get the dwell times for this needle
    % dwell_end_idx = dwell_start_idx + num_dwells - 1;
    % dwells = dwell_times(dwell_start_idx:dwell_end_idx);
    % dwell_start_idx = dwell_end_idx + 1;

    % Loop over each dwell position for this needle
    for j = 1:num_dwells
        % Get source position
        sourcePosition = test_tracks(i).active_dwells(j, :);

        % Loop over each voxel in the entire volume
        for xi = 1:gridSize(1)
            for yi = 1:gridSize(2)
                for zi = 1:gridSize(3)
                    % Linear index of the current voxel
                    voxel_idx = sub2ind(gridSize, xi, yi, zi);

                    % Get the grid point position
                    gridPoint = [xGrid(xi, yi, zi), yGrid(xi, yi, zi), zGrid(xi, yi, zi)];

                    % Compute distance and angle
                    r_vec = gridPoint - sourcePosition;
                    r = norm(r_vec);

                    % Skip zero-distance to avoid division by zero
                    if r == 0
                        doseRate = (3.34 + 2.17)/2; % average of adjacent positions from QA along and away table

                    else
                        % Calculate theta using dot product
                        dot = sum(conj(test_tracks(i).unit) .* r_vec);
                        theta = acosd(abs(dot) / r);

                        % Compute the dose contribution for this voxel and dwell
                        G_Lr_s = abs(G_L(r, L, theta) / G_L_r_0);
                        g_Lr_s = g_Lr(gLr_table, r);
                        F_rt_s = F_rt(F_table, r, theta);

                        doseRate = S_k * lambda * G_Lr_s * g_Lr_s * F_rt_s;

                        % Store in the unit dose matrix
                        test_unit_dose(voxel_idx, j) = doseRate;
                    end
                end
            end
        end
    end

    % Store the unit dose matrix in the needle structure
    test_tracks(i).unit_dose = test_unit_dose;

    % Display progress
    disp(['Needle ', num2str(i), '/', num2str(length(test_tracks)), ' (', num2str(i / length(test_tracks) * 100), '%) complete']);
end
%%
% Test tracks dwell time optimization 
dwell_times_test = optimize_dwell_times(test_tracks, [600, 1200], target_mask);
dose_test = dose_calc(test_tracks, dwell_times_test);


% %% Plotting
% 
% % Initialize variables to store slice data
% %%% Isodose levels %%%
% isodose_levels = [1260.6, 600, 588, 500, 400, 300, 200, 100, 50, 20];
% colors = [
%     0.3871, 0.9056, 0.7088;
%     0.5888, 0.1904, 0.1904;
%     0.1023, 0.8296, 0.5910;
%     0.6873, 0.0685, 0.9229;
%     0.7992, 0.2411, 0.2136;
%     0.2151, 0.3238, 0.5223;
%     0.4388, 0.3121, 0.6007;
%     0.1755, 0.3129, 0.3797;
%     0.4605, 0.7567, 0.2297;
%     0.5128, 0.5832, 0.0918];
% slice_figures = struct();
% c = 0;
% 
% % Store dose grid slices with isodose lines
% for s = 1:length(zRange)
%     % Create a new figure for each slice
%     fig = figure('Visible', 'off'); 
%     hold on;
% 
%     if -r_sphere_margin <= zRange(s) && zRange(s) <= r_sphere_margin
%         c = c + 1;
% 
%         % Plot target contour
%         plot3(contour_data_final(c).x, contour_data_final(c).y, contour_data_final(c).z, 'k-', 'LineWidth', 2);
%     end
% 
%     % Plot dose slice
%     slice(xGrid, yGrid, zGrid, dose_test, [], [], zRange(s));
%     shading interp;
%     colormap(jet);
%     colorbar;
% 
%     % Track which isodose levels have already been added to the legend
%     added_to_legend = false(1, length(isodose_levels));
%     visible_isodose_labels = {};
%     legend_handles = [];
% 
%     for i = 1:length(isodose_levels)
%         % Extract the current dose slice
%         dose_slice = squeeze(dose_test(:, :, s));
% 
%         % Compute the contour data for the current isodose level
%         temp_fig = figure('Visible', 'off');
%         [C, ~] = contour(xRange, yRange, dose_slice, [isodose_levels(i) isodose_levels(i)]);
%         close(temp_fig);
% 
%         % Check if the contour data is valid
%         if isempty(C)
%             continue; % Skip this isodose level if no contour is found
%         end
% 
%         % Parse the contour matrix
%         k = 1; % Starting index in C
%         while k < size(C, 2)
%             level = C(1, k); % Isodose level (should match isodose_levels(i))
%             num_points = C(2, k); % Number of points in the segment
% 
%             % Extract the segment points
%             x_contour = C(1, k+1:k+num_points);
%             y_contour = C(2, k+1:k+num_points);
% 
%             % Skip the contour if it lies entirely outside the current plane
%             if all(x_contour < min(xRange) | x_contour > max(xRange)) || ...
%                all(y_contour < min(yRange) | y_contour > max(yRange))
%                 k = k + num_points + 1; % Move to the next contour segment
%                 continue;
%             end
% 
%             % Plot the contour at the correct Z level
%             z_contour = zRange(s) * ones(size(x_contour)); % Set Z to the current slice level
%             handle = plot3(x_contour, y_contour, z_contour, 'Color', colors(i, :), 'LineWidth', 2);
% 
%             % Add to the legend only if it hasn't been added yet
%             if ~added_to_legend(i)
%                 visible_isodose_labels{end+1} = [num2str(isodose_levels(i))];
%                 legend_handles = [legend_handles, handle];
%                 added_to_legend(i) = true; % Mark this level as added
%             end
% 
%             % Move to the next segment
%             k = k + num_points + 1;
%         end
%     end
% 
%     % Plot needles and dwell positions
%     for i = 1:length(test_tracks)
%         track = test_tracks(i).track;
%         plot3(track(:, 1), track(:, 2), track(:, 3), 'r-', 'LineWidth', 1.5);
%         plot3(test_tracks(i).active_dwells(:, 1), test_tracks(i).active_dwells(:, 2), ...
%               test_tracks(i).active_dwells(:, 3), 'ro', 'MarkerSize', 6);
%     end
% 
%     % Customize axis and labels
%     hold off;
%     xlim([-8, 8]); ylim([-8, 8]); zlim([-8, 8]);
%     xlabel('X'); ylabel('Y'); zlabel('Z');
% 
%     % Set the legend dynamically
%     if ~isempty(legend_handles)
%         legend(legend_handles, visible_isodose_labels, 'FontSize', 16, 'Location', 'best');
%     end
% 
%     view(3);
%     rotate3d on;  % Enable 3D rotation
% 
%     % Store the figure handle in the struct (instead of the image)
%     slice_figures(s).fig_handle = fig;
%     slice_figures(s).title = ['Dose Slice at Z = ', num2str(zRange(s))];
% 
%     % Pause to allow for figure rendering
%     pause(0.1);
% end
% 
% % Display the first slice initially
% current_index = 1;
% current_fig = slice_figures(current_index).fig_handle;  % Retrieve the first figure handle
% figure(current_fig); % Bring the figure to the front
% 
% % Add a slider to navigate the slices
% slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', length(zRange), ...
%     'Value', 1, 'SliderStep', [1/(length(zRange)-1), 10/(length(zRange)-1)], ...
%     'Position', [20, 20, 300, 20], ...
%     'Callback', @(src, event) updateSliceWithIsoDose(src, slice_figures));
% 
% % Callback function to update the displayed slice
% function updateSliceWithIsoDose(slider, slice_figures)
%     % Get the current index from the slider
%     index = round(slider.Value);
% 
%     % Get the figure handle for the selected slice
%     selected_fig = slice_figures(index).fig_handle;
% 
%     % Bring the selected slice figure to the front
%     figure(selected_fig); 
% 
%     % Update the title
%     title(slice_figures(index).title);
% 
%     % Refresh the display (this allows for 3D rotation to be retained)
%     drawnow;
% end