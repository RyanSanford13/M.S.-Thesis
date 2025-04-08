function points = generatePointsAlongLine(P1, P2, increment)
    % generatePointsAlongLine: Generates a set of points along a line
    % defined by two points, spaced by a given increment.
    %
    % INPUTS:
    %   P1        - A 1x3 vector representing the starting point [x1, y1, z1].
    %   P2        - A 1x3 vector representing the ending point [x2, y2, z2].
    %   increment - The distance between consecutive points.
    %
    % OUTPUT:
    %   points - An Nx3 matrix where each row is a point [x, y, z] along the line.

    % Calculate the direction vector and its magnitude
    direction_vector = P2 - P1;
    line_length = norm(direction_vector);

    % Normalize the direction vector
    unit_vector = direction_vector / line_length;

    % Number of points (from start to end with increments)
    num_points = floor(line_length / increment) + 1;

    % Initialize array to store the points
    points = zeros(num_points, 3);
    
    % Generate points along the line
    for i = 0:num_points-1
        if (P1 - i * increment * unit_vector) < 0
            points(i+1, :) = P1 - i * increment * unit_vector;
        else
            points(i+1, :) = P1 + i * increment * unit_vector;
        end
    end
    
    % Add the ending point explicitly if not already included
    if norm(points(end, :) - P2) > 1e-5  % Check if it's already included
        points = [points; P2];
    end
end