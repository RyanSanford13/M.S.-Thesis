%% Function to calculate the resultant dose distribution 
function dose_distribution = dose_calc(needles, dwell_times, xGrid)


dose_distribution = zeros(size(xGrid));
dwell_start = 1;

for i = 1 : length(needles)
    % Get the number of active dwells for the current needle
    num_active_dwells = size(needles(i).active_dwells, 1);
    dwell_end = dwell_start + num_active_dwells - 1;
    
    % % Debug information to verify indices
    % fprintf('Needle: %d, dwell_start: %d, dwell_end: %d, num_active_dwells: %d\n', ...
    %         i, dwell_start, dwell_end, num_active_dwells);
    
    % Ensure indices do not exceed bounds of dwell_times
    if dwell_start > length(dwell_times) || dwell_end > length(dwell_times)
        error('Index exceeds bounds: dwell_start = %d, dwell_end = %d, length(dwell_times) = %d', ...
              dwell_start, dwell_end, length(dwell_times));
    end
    
    % Extract dwell times
    needle_times = dwell_times(dwell_start:dwell_end);
    dwell_start = dwell_end + 1;
    
    % Compute dose contribution
    needle_dose = needles(i).unit_dose * (needle_times / 3600);
    dose_distribution = dose_distribution + reshape(needle_dose, size(xGrid));
end
