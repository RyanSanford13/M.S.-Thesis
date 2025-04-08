function isInside = isPointInSphere(point, center, radius)
    % Calculate the squared distance between the point and the center of the sphere
    distanceSquared = (point(1) - center(1))^2 + (point(2) - center(2))^2 + (point(3) - center(3))^2;
    
    % Compare the squared distance with the squared radius (to avoid sqrt computation)
    isInside = distanceSquared <= radius^2;
end