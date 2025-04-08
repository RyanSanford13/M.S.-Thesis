%% Function to evaluate a given set of needles' satisfaction of the dose criteria
% update with other objectives as needed 
function satisfaction = dose_eval(needles, dwell_times, target_mask, rectum_mask, target_dose_range, xGrid, yGrid, zGrid, weights)

    % dose criteria
    dose = dose_calc(needles, dwell_times, xGrid);
    
    target_dose = dose(target_mask == 1);
    target_dose = sort(target_dose, 'ascend');
    %dose_range = max(target_dose) - min(target_dose);

    num_targ_voxels = numel(target_dose);

    rectum_dose = dose(rectum_mask == 1);


    % Comment out when doing anatomy 
    % weighting criteria (for sphere)
    % radius = 2; % hard coded could pass this in
    % distance_map = sqrt(xGrid.^2 + yGrid.^2 + zGrid.^2);
    % weight_map = 1 - abs(distance_map/radius - .5);
    % weight_map(~target_mask) = 0;
    % weight_map = weight_map(weight_map > 0);
    % 
    % weighted_dose = target_dose .* weight_map;

    %/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    % dose satisfaction -> evaluate normalized dose confromity such that more dose conformity => better satsifaction
    % point dose satsifaction 
    % norm_dose_satisfaction = sum(target_dose(:) >= target_dose_range(1) & target_dose(:) <= target_dose_range(2))/num_targ_voxels;
    % Dose volume satisfaction 
    D90_idx = round(0.1*num_targ_voxels);
    D90 = target_dose(D90_idx);
    
    % Sigmoid function for D90 satisfaction
    k_D90 = 0.1;  % Steepness factor (adjustable)
    tau_D90 = 50; % Tolerance threshold (adjustable, in cGy)
    D90_target = target_dose_range(1);
    
    % Satisfaction grows smoothly as D90 approaches the target range
    norm_D90_satisfaction = 1 / (1 + exp(k_D90 * (abs(D90 - D90_target) - tau_D90)));

    %///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % volume satisfaction 
    k_V200 = 20;
    tau_V200 = .5;
    
    V200_threshold = 2 * target_dose_range(1); % 200% of prescription dose
    V200_voxels = sum(target_dose >= V200_threshold);
    V200 = (V200_voxels / num_targ_voxels) * 100; % V200 as percentage of target volume

    norm_V200_satisfaction = 1 / (1 + exp(k_V200 * (V200 - tau_V200)));

    %///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    % rectum satisfaction 
    voxelVolume = .9375/10 * .9375/10 * 3/10;
    N2cc = round(2/voxelVolume);

    rectum_dose = sort(rectum_dose, 'descend');
    rectum_d2cc = rectum_dose(N2cc);

    alpha = 1.5; % controls steepness of fall off

    % 395 = midpoint at which satisfaction = .5
    norm_rectum_satisfaction = 1 / (1 + exp(alpha * (rectum_d2cc - 395)));
    %///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    % homogeneity satisfaction -> evaluate normalized variance (inverted) such that less variance => better satisfaction
    % homogeneity = var(target_dose, 0, 'all');
    % homogeneity_max = max(dose_range, 1e-6);
    % homogeneity_satisfaction = 1 - (homogeneity / homogeneity_max);
    % norm_homogeneity_satisfaction = max(0, min(1, homogeneity_satisfaction));

    % weighting satisfaction
    %norm_weight_satsifaction = sum(weighted_dose(:) >= target_dose_range(1) & weighted_dose(:) <= target_dose_range(2))/num_targ_voxels;

    % overall satisfaction weighted for each metric 
    w_d = weights(1);
    %w_h = weights(2);
    w_v = weights(2);
    %w_w = weights(3);
    w_r = weights(3);
    satisfaction = w_d * norm_D90_satisfaction+ w_v * norm_V200_satisfaction + w_r * norm_rectum_satisfaction; %+ w_w * norm_weight_satsifaction;
 
end