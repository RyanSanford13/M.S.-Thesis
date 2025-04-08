%% Function to uniformly sample a given triangle
% tri_verts is a matrix of triangle vertices
function uniformSample = uniformSample(tri_verts)

% Generate random numbers for point selection 
% generate random numbers
p_x = rand();
p_y = rand();
p_z = tri_verts(1,3);

% restrictions on random variables 
if (p_x + p_y > 1)
    p_x = 1 - p_x;
    p_y = 1 - p_y;
end

% Find random point via uniform sampling 
point = [];
for i = 1:3
    p = (1 - sqrt(p_x))*tri_verts(1,i) + (sqrt(p_x)*(1-p_y))*tri_verts(2,i)+(sqrt(p_x)*p_y)*tri_verts(3,i);
    point = [point, p];
end

uniformSample = point;
end