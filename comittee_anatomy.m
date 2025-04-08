%%% Load anatomy in %%%

rtStruct_info = dicominfo('9_6_24_MR_Struct_ang/RS.qCRQMi6g4AYjASZTXUYnXMY3z.9424_axt2 3d.dcm');
mr_info = dicominfo("9_6_24_MR_Struct_ang/MR.qCRQMi6g4AYjASZTXUYnXMY3z.Image 1.dcm");
image_height = mr_info.Height;

loadedData = load('cerv_hrctv.mat');
cerv_hrctv = loadedData.cerv_hrctv;

loadedData = load('uterus.mat');
uterus = loadedData.uterus;

loadedData = load('rectum.mat');
rectum = loadedData.rectum;

loadedData = load('vagina.mat');
vagina = loadedData.vagina;

loadedData = load('hrctv_marg.mat');
hrctv_marg = loadedData.hrctv_marg;

% Plot Anatomy
figure;
hold on;
trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', 1);
%trisurf(hrctv_marg.faces, hrctv_marg.verts(1,:), hrctv_marg.verts(2,:), hrctv_marg.verts(3,:), 'FaceColor', 'c', 'facealpha', .25);
trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', 1);
axis equal;
view(0, 270)
hold off;

%%
% Create figure
fig = figure;
ax = axes(fig);

hold on;

trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', 1);
trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', 1);

xlim([-20 40])
ylim([15 130])
zlim([-90 90])
axis equal;

% Set up video writer
v = VideoWriter('rotation_video.mp4', 'MPEG-4');
v.FrameRate = 30; % Adjust frame rate as needed
open(v);

% Rotate the view and capture frames
for angle = 1:2:360
    view(ax, angle, 30); % Change camera angle
    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);
end

% Close the video writer
close(v);
disp('Video saved as rotation_video.mp4');

%%
% Create entry zone

% define entry zone at halfway along vaginal canal
half_way_inner = sprintf('Item_%d', round(length(fieldnames(rtStruct_info.ROIContourSequence.Item_6.ContourSequence))/2));
z_coord = getfield(rtStruct_info.ROIContourSequence.Item_6.ContourSequence, half_way_inner);
z_coord = reshape(z_coord.ContourData, 3, []).';
z_coord = z_coord(1, 3);

entry_coords = getfield(rtStruct_info.ROIContourSequence.Item_1.ContourSequence, 'Item_1');

entry_coords = reshape(entry_coords.ContourData, 3, []);
entry_coords = entry_coords.';

z_coord = repmat(z_coord, length(entry_coords(:, 1)), 1);

entry_zone = [entry_coords(:,1), entry_coords(:,2), z_coord];

% Find convex hull of E, the entry zone 
convex_entry = convhull(entry_zone(:, 1:2));
convexE = [entry_zone(convex_entry, 1), entry_zone(convex_entry, 2), entry_zone(convex_entry,3)];

% Plot the entry zone and convex hull
figure;
hold on;
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);
plot3(convexE(:, 1), convexE(:, 2), convexE(:, 3), 'r-', 'LineWidth', 2);

trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', 1);
trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', 1);

% Add labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D View of Entry Zone and Convex Hull');

% Set axis properties
axis equal;
grid on;

% Set a 3D view angle (optional)
view(3); % Default 3D view

hold off;

%%
% Create target zone 

% set target zone 1 cm above last slice of target 
last = sprintf('Item_%d', length(fieldnames(rtStruct_info.ROIContourSequence.Item_1.ContourSequence)));

z_target = getfield(rtStruct_info.ROIContourSequence.Item_1.ContourSequence, last);
z_target = (reshape(z_target.ContourData, 3, [])).';
z_target = z_target(1,3) + 5;

ROI = rtStruct_info.ROIContourSequence.Item_1.ContourSequence.Item_1.ContourData;
ROI = reshape(ROI, 3, []).';

z_factor = ((z_target - entry_zone(1, 3)) / (ROI(1,3) - entry_zone(1, 3)));
target_zone = [];
for i = 1:length(entry_zone(:,1))
    for j = 1:length(ROI(:,1))
        %t = entry_zone(i, :) + z_factor * ([(ROI(1).layer(j,1)), (ROI(1).layer(j,2)), ROI(1).layer(j,3)] - entry_zone(i, :));
        t = [entry_zone(i, 1)/7, entry_zone(i, 2)/7, entry_zone(i, 3)] + z_factor * ([(ROI(j,1))/7, (ROI(j,2))/7, ROI(j,3)] - [entry_zone(i, 1)/7, entry_zone(i, 2)/7, entry_zone(i, 3)]);
        target_zone = [target_zone; t];
    end
end

% Find convex hull of T
convex_target = convhull(target_zone(:, 1:2));
convexT = [target_zone(convex_target, 1) + 5, target_zone(convex_target, 2) + 53, target_zone(convex_target, 3)];

figure;
hold on;
%plot3(target_zone(:, 1), target_zone(:, 2), target_zone(:, 3), 'b-', 'LineWidth', 3);
plot3(convexT(:, 1), convexT(:, 2), convexT(:, 3), 'r-', 'LineWidth', 2);
plot3(convexE(:, 1), convexE(:, 2), convexE(:, 3), 'b-', 'LineWidth', 2);

trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', 1);
% trisurf(hrctv_marg.faces, hrctv_marg.verts(1,:), hrctv_marg.verts(2,:), hrctv_marg.verts(3,:), 'FaceColor', 'c', 'facealpha', .25);
% trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
% trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', 1);

% Add labels and title
xlabel('X');
ylabel('Y');
zlabel('Z');
%title('3D View of Entry Zone and Target Zone Relative to the HR-CTV');

% Set axis properties
axis equal;
grid on;

% Set a 3D view angle (optional)
view(3); % Default 3D view

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
plot3(convexE(:, 1), convexE(:, 2), convexE(:, 3), 'r-', 'LineWidth', 2);

% Plot the convex hull of T
plot3(convexT(:, 1), convexT(:, 2), convexT(:, 3), 'r-', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the line connecting the points
for i=1:length(lines)
    plot3(lines(i).line(:,1), lines(i).line(:,2), lines(i).line(:,3), 'k-');
end

trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'EdgeColor', 'none', 'facealpha', 1);
trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'EdgeColor', 'none', 'facealpha', .5);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'EdgeColor', 'none', 'facealpha', .75);
trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'EdgeColor', 'none', 'facealpha', .75);

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Entry Convex, Target Convex, Anatomy, Target and Entry Triangulation, and 10000 needles');

%%
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
%plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(convexE(:, 1), convexE(:, 2), convexE(:, 3), 'b-', 'LineWidth', 2);

% Plot the target (target includes margin)
trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'c', 'EdgeColor', 'none', 'facealpha', 0.5);

% Plot the convex hull of T
plot3(convexT(:, 1), convexT(:, 2), convexT(:, 3), 'r-', 'LineWidth', 2);

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
% define needle tracks and important components

needle_tracks = struct('track', []);

for i = 1:length(lines)
    needle_tracks(i).track = lines(i).line;
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

ang_threshhold = 12.5; % degrees
needle_tracks = needle_tracks(angles <= ang_threshhold);

numFaces = size(cerv_hrctv.faces, 1);

cerv_hrctv.verts = cerv_hrctv.verts';
cerv_hrctv.verts(:, 4) = [];

% Prepare triangle vertices (vert0, vert1, vert2)
vert0 = cerv_hrctv.verts(cerv_hrctv.faces(:, 1), :);  % Mx3 matrix (vertices of each triangle)
vert1 = cerv_hrctv.verts(cerv_hrctv.faces(:, 2), :);  % Mx3 matrix (vertices of each triangle)
vert2 = cerv_hrctv.verts(cerv_hrctv.faces(:, 3), :);  % Mx3 matrix (vertices of each triangle)

keepTracks = true(length(needle_tracks), 1);
for i = 1:length(needle_tracks)
    orig = needle_tracks(i).track(1,:);
    direc = needle_tracks(i).unit;

    orig = repmat(orig, numFaces, 1);
    direc = repmat(direc, numFaces, 1);

    [intersect, t, u, v, xcoor] = TriangleRayIntersection(orig, direc, vert0, vert1, vert2);

   % If no intersection, mark the track for removal
    if ~any(intersect)
        keepTracks(i) = false;
    end
end 

needle_tracks = needle_tracks(keepTracks);

figure;
hold on;

% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(convexE(:, 1), convexE(:, 2), convexE(:, 3), 'r-', 'LineWidth', 2);

% Plot the convex hull of T
plot3(convexT(:, 1), convexT(:, 2), convexT(:, 3), 'r-', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

for i=1:length(needle_tracks)
    plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'k-');
end

% turn verts back
cerv_hrctv.verts = cerv_hrctv.verts.';

trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', 1);
trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', 1);

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Entry Convex, Target Convex, Anatomy, Target and Entry Triangulation, and 10000 needles');

%%
%%% define needle dwells 
% Recall all distances in mm
needle_offset = 7;

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
    dwells = generatePointsAlongLine(needle_tracks(i).first_dwell, last_dwell, 5);
    needle_tracks(i).dwells = dwells;
end

for i = 1:length(needle_tracks)
    in = inpolyhedron(hrctv_marg.faces, hrctv_marg.verts(1:3, :).', needle_tracks(i).dwells);
    needle_tracks(i).active_dwells = needle_tracks(i).dwells(in, :);
end

% remove dwells below start of cervix
thresh = rtStruct_info.ROIContourSequence.Item_4.ContourSequence.Item_1.ContourData(3); % Define the threshold value

for i = 1:length(needle_tracks)
    % Get the active dwell positions for the current needle
    active_dwells = needle_tracks(i).active_dwells;
    
    % Keep only the rows where the third column (z-coordinate) is >= x
    needle_tracks(i).active_dwells = active_dwells(active_dwells(:, 3) >= thresh, :);
end

% remove dwells that extend past 5mm of end of cervix
thresh = rtStruct_info.ROIContourSequence.Item_4.ContourSequence.Item_24.ContourData(3) + 5;
for i = 1:length(needle_tracks)
    % Get the active dwell positions for the current needle
    active_dwells = needle_tracks(i).active_dwells;
    
    % Keep only the rows where the third column (z-coordinate) is <= x
    needle_tracks(i).active_dwells = active_dwells(active_dwells(:, 3) <= thresh, :);
end

with_adwells = false(length(needle_tracks), 1);
for i = 1:length(needle_tracks)
    if ~isempty(needle_tracks(i).active_dwells)
        with_adwells(i) = 1;
    end
end

needle_tracks = needle_tracks(with_adwells);

numNeedles = length(needle_tracks);


% Remove any needles within 5 mm of eachother
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
                if any(distances(:) < 5)
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
plot3(convexE(:, 1), convexE(:, 2), convexE(:, 3), 'r-', 'LineWidth', 2);

% Plot the convex hull of T
plot3(convexT(:, 1), convexT(:, 2), convexT(:, 3), 'r-', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the line connecting the points
for i=1:length(needle_tracks)
    plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'k-', 'LineWidth', 2);
    plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
    for j = 1:size(needle_tracks(i).active_dwells, 1)
        plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'm.', 'MarkerSize', 10);
    end
end

trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', 1);
%trisurf(hrctv_marg.faces, hrctv_marg.verts(1,:), hrctv_marg.verts(2,:), hrctv_marg.verts(3,:), 'FaceColor', 'c', 'facealpha', .25);
trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', .5, 'EdgeColor', 'none');

axis equal;
axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Entry Convex, Target Convex, Anatomy, Target and Entry Triangulation, and 10000 needles');
hold off;

%% 
% Add central tandem to the needle set 
% get relevant needle structure from another needle 
tandem = needle_tracks(1);
tandem.idx = 1;
tandem.track = [];
tandem.unit = [];
tandem.first_dwell = [];
tandem.dwells = [];
tandem.active_dwells = [];

center = [mean(convexE(:, 1)), mean(convexE(:, 2)), convexE(1, 3);...
    mean(convexE(:, 1)), mean(convexE(:, 2)), max(cerv_hrctv.verts(3, :)) + 5];

tandem.track = center;

% Get starting and ending points 
start = tandem.track(2,:);
ending = tandem.track(1,:);

% Define unit vector
unit = start - ending;
unit_norm = unit / norm(unit);
% store
tandem.unit = unit_norm;

scale_unit_norm = tandem.unit * 7;

% Compute first dwell position
first_dwell = start - scale_unit_norm;
% store
tandem.first_dwell = first_dwell;

% compute all dwells 
dwells = [];
% Recall all distances in mm
last_dwell = tandem.track(1, :);
dwells = generatePointsAlongLine(tandem.first_dwell, last_dwell, 5);
tandem.dwells = dwells;

% remove dwells that dont "intersect" hrctv
tandem.dwells(tandem.dwells(:, 3) < min(cerv_hrctv.verts(3, :)), :) = [];

tandem.active_dwells = tandem.dwells(2:8, :);

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

needle_tracks_temp = needle_tracks;

needle_tracks(1) = needle_tracks_temp(end);
needle_tracks(2:end) = needle_tracks_temp(1:end-1);

figure;
hold on;
% Plot the circle
plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
plot3(convexE(:, 1), convexE(:, 2), convexE(:, 3), 'r-', 'LineWidth', 2);

% Plot the convex hull of T
plot3(convexT(:, 1), convexT(:, 2), convexT(:, 3), 'r-', 'LineWidth', 2);

% Plot the entry triangulated mesh
trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the line connecting the points
for i=1:length(needle_tracks)
    if i == 1
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'r-', 'LineWidth', 2);
        plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
        for j = 1:size(needle_tracks(i).active_dwells, 1)
            plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'm.', 'MarkerSize', 10);
        end
    else
        plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'k-', 'LineWidth', 2);
        plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
        for j = 1:size(needle_tracks(i).active_dwells, 1)
            plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'm.', 'MarkerSize', 10);
        end
    end
end

trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', 1);
trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', 1);

axis equal;
%axis vis3d;

xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Entry Zone, Entry Convex, Target Convex, Anatomy, Target and Entry Triangulation, and 10000 needles');
hold off;

%%
figure;
hold on;
% Plot the circle
% plot3(entry_zone(:, 1), entry_zone(:, 2), entry_zone(:, 3), 'b-', 'LineWidth', 3);

% Plot the entry convex hull
% plot3(convexE(:, 1), convexE(:, 2), convexE(:, 3), 'k-', 'LineWidth', 2);

% Plot the convex hull of T
% plot3(convexT(:, 1), convexT(:, 2), convexT(:, 3), 'k-', 'LineWidth', 2);

% Plot the entry triangulated mesh
% trisurf(dtE.ConnectivityList, convexE(:,1), convexE(:,2), convexE(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% Plot the target triangulated mesh
% trisurf(dtT.ConnectivityList, convexT(:,1), convexT(:,2), convexT(:,3), 'FaceAlpha', 0.5, 'EdgeColor', 'black');

% % Plot the line connecting the points
% for i=1:length(needle_tracks)
%     if i == 1
%         plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'r-', 'LineWidth', 2);
%         plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
%         for j = 1:size(needle_tracks(i).active_dwells, 1)
%             plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'm.', 'MarkerSize', 10);
%         end
%     else
%         plot3(needle_tracks(i).track(:,1), needle_tracks(i).track(:,2), needle_tracks(i).track(:,3), 'k-', 'LineWidth', 2);
%         plot3(needle_tracks(i).first_dwell(1,1), needle_tracks(i).first_dwell(1,2), needle_tracks(i).first_dwell(1,3), 'go', 'MarkerSize', 10);
%         for j = 1:size(needle_tracks(i).active_dwells, 1)
%             plot3(needle_tracks(i).active_dwells(j, 1), needle_tracks(i).active_dwells(j, 2), needle_tracks(i).active_dwells(j, 3), 'm.', 'MarkerSize', 10);
%         end
%     end
% end
idx = 21;
plot3(needle_tracks(idx).track(:,1), needle_tracks(idx).track(:,2), needle_tracks(idx).track(:,3), 'k-', 'LineWidth', 2);
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


trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:),'EdgeColor', 'none', 'FaceColor', 'b', 'facealpha', .1);
% trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
% trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
% trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', 1);

axis equal;
% axis vis3d;
view(3);
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
title('Plot Showing Reduced Candidate Set Relative to the Anatomy');
hold off;






%%
% Dose Calculation Grid Initialization

% Creat Dose grid based on size of masks
% find min and max extent of structures
loaded_data = load('in2di_rotMat.mat');
in2di_rotMat = loaded_data.in2di_rotationMat;

% Extract relevant info from the transformation matrix
dx = in2di_rotMat(1,1);
dy = in2di_rotMat(2,2);
dz = in2di_rotMat(3,3);
patient_position = in2di_rotMat(1:3, 4);
margin = 10;

% Define bounding boxes for each of the structures
[rowIdx, colIdx, sliceIdx] = ind2sub(size(cerv_hrctv.mask), find(cerv_hrctv.mask));

xBound = [min(colIdx) - margin, max(colIdx) + margin];
yBound = [min(rowIdx) - margin, max(rowIdx) + (margin+20)];
zBound = [min(sliceIdx) - margin, max(sliceIdx) + margin];


[xGrid, yGrid, zGrid] = meshgrid(...
    patient_position(1) + (xBound(1) - 1:xBound(2) - 1)*dx, ...
    patient_position(2) + (yBound(1) - 1:yBound(2) - 1)*dy, ...
    patient_position(3) + (zBound(1) - 1:zBound(2) - 1)*dz);
%%
x_range = patient_position(1) + (xBound(1) - 1:xBound(2) - 1)*dx;
y_range = patient_position(2) + (yBound(1) - 1:yBound(2) - 1)*dy;
z_range = patient_position(3) + (zBound(1) - 1:zBound(2) - 1)*dz;

figure; 
    hold on; 
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Anatomical Meshes with Bounding Box');

    trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', 1);
    trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
    trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
    trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', 1);

    % Define bounding box corners
    X = [x_range(1) x_range(end) x_range(end) x_range(1) x_range(1) x_range(end) x_range(end) x_range(1)];
    Y = [y_range(1) y_range(1) y_range(end) y_range(end) y_range(1) y_range(1) y_range(end) y_range(end)];
    Z = [z_range(1) z_range(1) z_range(1) z_range(1) z_range(end) z_range(end) z_range(end) z_range(end)];

    % Define bounding box edges
    edges = [1 2; 2 3; 3 4; 4 1; ...  % Bottom face edges
             5 6; 6 7; 7 8; 8 5; ...  % Top face edges
             1 5; 2 6; 3 7; 4 8];     % Vertical edges

    % Plot bounding box edges
    for i = 1:size(edges, 1)
        plot3(X(edges(i, :)), Y(edges(i, :)), Z(edges(i, :)), 'k', 'LineWidth', 2);
    end

    view(3); grid on;
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
        sourcePosition = .1*needle_tracks(i).active_dwells(j, :);

        % Loop over each voxel in the entire volume
        for xi = 1:gridSize(1)
            for yi = 1:gridSize(2)
                for zi = 1:gridSize(3)
                    % Linear index of the current voxel
                    voxel_idx = sub2ind(gridSize, xi, yi, zi);

                    % Get the grid point position
                    gridPoint = .1*[xGrid(xi, yi, zi), yGrid(xi, yi, zi), zGrid(xi, yi, zi)];

                    % Compute distance and angle
                    r_vec = double(gridPoint - sourcePosition);
                    r = double(norm(r_vec));

                    % Skip zero-distance to avoid division by zero
                    if r > 0
                        % Calculate theta using dot product
                        dot = sum(conj(needle_tracks(i).unit) .* r_vec);
                        theta = double(acosd(abs(dot) / r));

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

%%
% Needle Selection and Optimization 

Rx = 600; % cGy
dose_min = Rx; % 600 cGy
dose_max = 750; % 750 cGy 
target_dose_range = [dose_min, dose_max];

target_mask = cerv_hrctv.mask(yBound(1):yBound(2), xBound(1):xBound(2), zBound(1):zBound(2));
rectum_mask = rectum.mask(yBound(1):yBound(2), xBound(1):xBound(2), zBound(1):zBound(2));

[selected_needles_opt, dwell_times_opt] = greedy_needle_selection(needle_tracks, target_dose_range, 15, target_mask, rectum_mask, xGrid, yGrid, zGrid);
disp("Optimization Done, starting final dose calc...")
dose_needles_opt = dose_calc(selected_needles_opt, dwell_times_opt, xGrid);
disp("Final dose calc done, continue to plotting");

%%
save('selected_needles_BEST_anat_22025.mat', 'selected_needles_opt');
save('dwell_times_BEST_anat_22025.mat', 'dwell_times_opt');
save('dose_needles_BEST_anat_22025.mat', 'dose_needles_opt');

%%
dwell_times_opt_redo = optimize_dwell_times(selected_needles_opt3, [692, 700], target_mask, rectum_mask, xGrid);
dose_needles_opt_redo = dose_calc(selected_needles_opt3, dwell_times_opt_redo, xGrid);

%%
save('Ang_vol_results/selected_needles_opt_3.mat', 'selected_needles_opt3');
save('Ang_vol_results/dwell_times_opt_3.mat', 'dwell_times_opt_redo');
save('Ang_vol_results/dose_needles_opt_3.mat', 'dose_needles_opt_redo');

%%
target_mask = cerv_hrctv.mask(yBound(1):yBound(2), xBound(1):xBound(2), zBound(1):zBound(2));
rectum_mask = rectum.mask(yBound(1):yBound(2), xBound(1):xBound(2), zBound(1):zBound(2));
Rx = 600;

target_dose = dose_needles_opt_redo(target_mask == 1);
target_dose = sort(target_dose, 'ascend');

% Compute volume corresponding to each dose level
numVoxels = numel(target_dose);
voxelVolume = (dx/10*dy/10*dz/10);
totalVolume = numVoxels * ((dx/10*dy/10*dz/10)); % Total target volume

volumeRemaining = (numVoxels:-1:1) * (dx/10*dy/10*dz/10); % Absolute volume
percentvolumeRemaining = volumeRemaining / totalVolume;

D150 = 1.5 * Rx;
D200 = 2 * Rx;

numVoxels_V150 = sum(target_dose(:) >= D150);
numVoxels_V200 = sum(target_dose(:) >= D200);

% Compute volumes
V150_abs = numVoxels_V150 * ((dx/10*dy/10*dz/10)); % Absolute volume (cc)
V150_rel = (V150_abs / totalVolume) * 100; % Relative volume (%)

V200_abs = numVoxels_V200 * ((dx/10*dy/10*dz/10)); % Absolute volume (cc)
V200_rel = (V200_abs / totalVolume) * 100; % Relative volume (%)

% Get D90 (dose at 90% volume)
index = round(0.1 * numVoxels);  % 10% of the smallest doses
D90 = target_dose(index);

% rectum d2cc
N_2cc = round(2 / voxelVolume);
rectum_dose = dose_needles_opt_redo(rectum_mask == 1);
rectum_dose = sort(rectum_dose, 'descend');

rectum_d2cc = rectum_dose(N_2cc);

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
disp(['D2cc, rectum:', num2str(rectum_d2cc)]);

%%
loaded_data = load("Ang_vol_results/dose_needles_opt_1.mat");
dose_needles_opt1 = loaded_data.dose_needles_opt;

loaded_data = load("Ang_vol_results/selected_needles_opt_1.mat");
selected_needles_opt1 = loaded_data.selected_needles_opt;

%%
figure;
hold on;

color1 = [0.3871, 0.9056, 0.7088];
color2 = [0.4605, 0.7567, 0.2297];

isosurf = isosurface(xGrid, yGrid, zGrid, dose_needles_opt1, 600);
isosurf1 = isosurface(xGrid, yGrid, zGrid, dose_needles_opt1, 390);
p = patch(isosurf, 'FaceColor', color1, 'EdgeColor', 'none', 'FaceAlpha', .75);
p1 = patch(isosurf1, 'FaceColor', color2, 'EdgeColor', 'none', 'FaceAlpha', .25);

trisurf(cerv_hrctv.faces, cerv_hrctv.verts(1,:), cerv_hrctv.verts(2,:), cerv_hrctv.verts(3,:), 'FaceColor', 'b', 'facealpha', .25);
%trisurf(hrctv_marg.faces, hrctv_marg.verts(1,:), hrctv_marg.verts(2,:), hrctv_marg.verts(3,:), 'FaceColor', 'c', 'facealpha', .25);
trisurf(uterus.faces, uterus.verts(1,:), uterus.verts(2,:), uterus.verts(3,:), 'FaceColor', 'r', 'facealpha', 1);
trisurf(rectum.faces, rectum.verts(1,:), rectum.verts(2,:), rectum.verts(3,:), 'FaceColor', 'g', 'facealpha', 1);
trisurf(vagina.faces, vagina.verts(1,:), vagina.verts(2,:), vagina.verts(3,:), 'FaceColor', 'y', 'facealpha', .5);

for i = 1:length(selected_needles_opt3)
        track = selected_needles_opt1(i).track;
        plot3(track(:, 1), track(:, 2), track(:, 3), 'r-', 'LineWidth', 2.5);
        plot3(selected_needles_opt1(i).active_dwells(:, 1), selected_needles_opt1(i).active_dwells(:, 2), ...
              selected_needles_opt1(i).active_dwells(:, 3), 'ro', 'MarkerSize', 6);
end

axis equal;
grid on;
hold off;
%%
%%% Isodose levels %%%
xRange = linspace(min(xGrid(:)), max(xGrid(:)), size(dose_needles_opt1, 2));
yRange = linspace(min(yGrid(:)), max(yGrid(:)), size(dose_needles_opt1, 1));
zRange = linspace(min(zGrid(:)), max(zGrid(:)), size(dose_needles_opt1, 3));
isodose_levels = [1260.6, 600, 588, 500, 400, 390, 200, 100, 50, 20];
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

    % Plot dose slice
    slice(xGrid, yGrid, zGrid, dose_needles_opt1, [], [], zRange(s));
    shading interp;
    colormap(jet);
    colorbar;


    % Overlay anatomical structure contours
    hrctv_mask = cerv_hrctv.mask(yBound(1):yBound(2), xBound(1):xBound(2), zBound(1):zBound(2));
    mask_slice_hrctv = squeeze(hrctv_mask(:, :, s));
    [B_h, ~] = bwboundaries(mask_slice_hrctv, 'noholes');

    for k = 1:length(B_h)
        boundary_h = B_h{k};
        plot3(xRange(boundary_h(:,2)), yRange(boundary_h(:,1)), zRange(s)*ones(size(boundary_h,1),1), 'y-', 'LineWidth', 2);
        
    end

    uterus_mask = uterus.mask(yBound(1):yBound(2), xBound(1):xBound(2), zBound(1):zBound(2));
    mask_slice_uterus = squeeze(uterus_mask(:, :, s));
    [B_u, ~] = bwboundaries(mask_slice_uterus, 'noholes');

    for k = 1:length(B_u)
        boundary_u = B_u{k};
        plot3(xRange(boundary_u(:,2)), yRange(boundary_u(:,1)), zRange(s)*ones(size(boundary_u,1),1), 'w-', 'LineWidth', 2);
    end

    rectum_mask = rectum.mask(yBound(1):yBound(2), xBound(1):xBound(2), zBound(1):zBound(2));
    mask_slice_rectum = squeeze(rectum_mask(:, :, s));
    [B_r, ~] = bwboundaries(mask_slice_rectum, 'noholes');

    for k = 1:length(B_r)
        boundary_r = B_r{k};
        plot3(xRange(boundary_r(:,2)), yRange(boundary_r(:,1)), zRange(s)*ones(size(boundary_r,1),1), 'g-', 'LineWidth', 2);
    end

    % Track which isodose levels have already been added to the legend
    added_to_legend = false(1, length(isodose_levels));
    visible_isodose_labels = {};
    legend_handles = [];

    dose_slice = squeeze(dose_needles_opt1(:, :, s));
    %[C, ~] = contour(xRange, yRange, dose_slice, isodose_levels);
    for i = 1:length(isodose_levels)
        % Extract the current dose slice
        dose_slice = squeeze(dose_needles_opt1(:, :, s));

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
            handle = plot3(x_contour, y_contour, z_contour, 'Color', colors(i, :), 'LineWidth', 1);

            % Add to the legend only if it hasn't been added yet
            if ~added_to_legend(i)
                visible_isodose_labels{end+1} = [num2str(isodose_levels(i)), ' cGy'];
                legend_handles = [legend_handles, handle];
                added_to_legend(i) = true; % Mark this level as added
            end

            % Move to the next segment
            k = k + num_points + 1;
        end
    end

    % Plot needles and dwell positions
    for i = 1:length(selected_needles_opt1)
        track = selected_needles_opt1(i).track;
        plot3(track(:, 1), track(:, 2), track(:, 3), 'r-', 'LineWidth', 2.5);
        plot3(selected_needles_opt1(i).active_dwells(:, 1), selected_needles_opt1(i).active_dwells(:, 2), ...
              selected_needles_opt1(i).active_dwells(:, 3), 'ro', 'MarkerSize', 6);
    end
    % for j = 1:length(needle_tracks)
    %     track1 = needle_tracks(j).track;
    %     plot3(track1(:, 1), track1(:, 2), track1(:,3), 'k-', 'LineWidth', 1.5);
    % end

    % Customize axis and labels

    xlim([min(xGrid(:)), max(xGrid(:))]); 
    ylim([min(yGrid(:)), max(yGrid(:))]); 
    zlim([min(zGrid(:)), max(zGrid(:))]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal
    hold off;

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


