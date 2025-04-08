%% Function checks if needle tracks intersects target
% line is the endpoints of the neeedle segment stored as a 2x3 matrix
% 
function checkIntersect = checkIntersect(line, C, R)
    % Extract coordinates
    % point 1
    x1 = line(1,1); 
    y1 = line(1,2); 
    z1 = line(1,3);
    
    % point 2
    x2 = line(2,1); 
    y2 = line(2,2); 
    z2 = line(2,3);
    
    % sphere center
    xc = C(1); 
    yc = C(2); 
    zc = C(3);
    
    % Calculate coefficients of the quadratic equation
    % line equation: P(t) = P1 + t(P2-P1) -> gives a point along the line defined by t in [0,1]
    % sphere equation: (x-xc)^2 + (y-yc)^2 + (z-zc)^2 = r^2 -> gives point on the sphere given r and the center 
    % substitute x,y,z for P(t) and solve for t
    % left with at^2 + bt + c = 0
    
    dx = x2 - x1; 
    dy = y2 - y1; 
    dz = z2 - z1;
    a = dx^2 + dy^2 + dz^2;
    b = 2 * (dx * (x1 - xc) + dy * (y1 - yc) + dz * (z1 - zc));
    c = (x1 - xc)^2 + (y1 - yc)^2 + (z1 - zc)^2 - R^2;
    
    % Solve the quadratic equation a*t^2 + b*t + c = 0
    discriminant = b^2 - 4*a*c;
    
    if discriminant < 0
        % No real solutions, the line does not intersect the sphere
        checkIntersect = false;
    else
        % Compute the two solutions for t
        t1 = (-b + sqrt(discriminant)) / (2 * a);
        t2 = (-b - sqrt(discriminant)) / (2 * a);
        
        % Check if either solution is within the range [0, 1]
        if (t1 >= 0 && t1 <= 1) || (t2 >= 0 && t2 <= 1)
            checkIntersect = true;
        else
            checkIntersect = false;
        end
    end
end