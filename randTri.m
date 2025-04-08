%% Function to randomly select a triangle given a triangulate polygon
%  conList is the connectivities produced by delaunayTriangulation()
%  convexHull is the convex hull of the triangulated polygon
function [randTri, idxTri] = randTri(conList, convexHull)

triangles = conList;
numTriangles = size(triangles, 1);
areas = zeros(numTriangles, 1);

for i = 1:numTriangles
    % Get the vertices of the triangle
    v1 = convexHull(triangles(i,1), :);
    v2 = convexHull(triangles(i,2), :);
    v3 = convexHull(triangles(i,3), :);
    
    % Compute the area using the cross product
    area = norm(cross(v2-v1, v3-v1)) / 2;
    areas(i) = area;
end

% Normalize the areas to get a probability distribution
totalArea = sum(areas);
probabilities = areas / totalArea;

% Compute the cumulative distribution function (CDF)
cdf = cumsum(probabilities);

% Randomly select a triangle based on the CDF
randomValue = rand();
selectedTriangleIndex = find(cdf >= randomValue, 1, 'first');

% Get the randomly selectred triangle
randTri = convexHull(triangles(selectedTriangleIndex, :), :);
idxTri = selectedTriangleIndex;
end