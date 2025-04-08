function dwell_times_opt = optimize_dwell_times(needles, target_dose_range, target_mask, rectum_mask, xGrid)
    % Optimize dwell times to satisfy the dose constraints using fmincon
    % INPUTS:
    %   - needles: array of needle structures with dose influence data
    %   - target_dose_range: [min_dose, max_dose] for target
    %   - target_mask: binary mask of the target volume
    %
    % OUTPUT:
    %   - dwell_times_opt: optimized dwell times for the needles


    % Number of dwells across all needles
    num_dwells = 0;
    for i = 1:length(needles)
        num_dwells = num_dwells + size(needles(i).active_dwells, 1); % Number of dwells for this needle
    end
    

    % Define the objective function
    function penalty = objective(dwell_times, xGrid)
        
        dose = dose_calc(needles, dwell_times, xGrid);
        if any(isnan(dose(:))) || any(isinf(dose(:)))
            disp('Warning: NaN or Inf encountered in dose calculation');
        end

        % Dose volume penalty function 
        target_dose = dose(target_mask == 1);

        numVoxels = numel(target_dose);
        voxelVolume = .9375/10 * .9375/10 * 3/10;
        totalVolume = voxelVolume * numVoxels;

        v200_dose = 1200;
        numVoxels_v200 = sum(target_dose(:) > v200_dose);

        v200_abs = numVoxels_v200 * voxelVolume;
        v200_rel = (v200_abs / totalVolume) * 100;

        D90_idx = round(0.1 * numVoxels); % D90 is the dose to the lowest 10% of voxels
        D90 = target_dose(D90_idx);

        % rectum penalty 
        N_2cc = round(2 / voxelVolume);
        rectum_dose = dose(rectum_mask == 1);
        rectum_dose = sort(rectum_dose, 'descend');

        rectum_d2cc = rectum_dose(N_2cc);

        % Penalize deviations from the target dose range
        underdose_penalty = max(0, target_dose_range(1) - D90); % Penalize if D90 < min target dose
        %overdose_penalty = max(0, D90 - target_dose_range(1));  % Penalize if D90 > max target dose
        volume_penalty = max(0, v200_rel - 10.0); % Penalize if V200 > 5 %
        rectum_penalty = max(0, rectum_d2cc - 390); % Penalize if rectum D2cc > 3.9 Gy

        % non_targ_mask = ~target_mask;
        % non_target_dose = dose(non_targ_mask == 1);
        % 
        % overdose_penalty_outside = sum(max(0, non_target_dose - target_dose_range(1)));
        
        penalty = underdose_penalty^2 + volume_penalty^2 + rectum_penalty^2; % Squared penalty for smoother optimization
        
        % Point dose penalty function
        % % Focus on target voxels only
        % target_dose = dose(target_mask == 1);
        % 
        % non_targ_mask = ~target_mask;
        % non_target_dose = dose(non_targ_mask == 1);
        % 
        % % % Penalize deviations from the target dose range
        % underdose_penalty = sum(max(0, target_dose_range(1) - target_dose)); % Penalize underdosing
        % overdose_penalty = sum(max(0, target_dose - target_dose_range(2)));  % Penalize overdosing
        % 
        % % Penalize overdoses outside the target
        % overdose_penalty_outside = sum(max(0, non_target_dose - target_dose_range(1)));
        % 
        % % Total penalty
        % penalty = underdose_penalty + overdose_penalty + overdose_penalty_outside;
    end

    % Define the constraints
    function [c, ceq] = constraints(dwell_times, xGrid)
        % Calculate dose distribution
        dose = dose_calc(needles, dwell_times, xGrid);
        
        % Focus on target voxels only
        target_dose = dose(target_mask ==1);
        target_dose = sort(target_dose, 'ascend');
        %non_targ_dose = dose(~target_mask);

        numVoxels = numel(target_dose);
        voxelVolume = .9375/10 * .9375/10 * 3/10;
        totalVolume = voxelVolume * numVoxels;

        v200_dose = 1200;
        numVoxels_v200 = sum(target_dose(:) > v200_dose);

        v200_abs = numVoxels_v200 * voxelVolume;
        v200_rel = (v200_abs / totalVolume) * 100;

        D90_idx = round(0.1 * numVoxels);
        D90 = target_dose(D90_idx);

        % rectum constraint
        % rectum penalty 
        N_2cc = round(2 / voxelVolume);
        rectum_dose = dose(rectum_mask == 1);
        rectum_dose = sort(rectum_dose, 'descend');

        rectum_d2cc = rectum_dose(N_2cc);


        % Enforce dose constraints: min_dose >= D90 
        c_lower = target_dose_range(1) - D90;
        %c_upper = D90 - target_dose_range(1);
        c_vol = v200_rel - 5.0;
        c_rectum = rectum_d2cc - 390;
        c =  [c_lower, c_vol, c_rectum];  % Ensure D90 ≥ min_dose and D90 ≤ max_dose
        
        % Hard constraints: Dose must stay within the range [min, max]
        % Point dose constraints
        % c_lower_targ = target_dose_range(1) - target_dose; % Ensure dose >= min_dose
        % c_upper_targ = target_dose - target_dose_range(2); % Ensure dose <= max_dose
        % %c_upper_nontarg = non_targ_dose - target_dose_range(1); % Ensure dose(non target) <= min_dose
        % c =  [c_lower_targ; c_upper_targ]; %c_upper_nontarg]; % Combine constraints

        %  % Time constraint: t_i <= 1.5 * t_j for adjacent dwell times
        % max_diff_penalty = 1.5;  % This corresponds to the "1.5" multiplier
        % 
        % % Loop through adjacent dwell times and add inequality constraints
        % for i = 1:length(dwell_times)-1
        %     % Add constraint: t_i - 1.5 * t_{i+1} <= 0 (i.e., t_i <= 1.5 * t_{i+1})
        %     c = [c; dwell_times(i) - max_diff_penalty * dwell_times(i+1)];
        % end

        ceq = [];
    end

    % Initial dwell times
    dwell_times_initial = ones(num_dwells, 1); % Uniform initial guess

    % Bounds for dwell times
    lb = 0.2 * ones(num_dwells, 1); % Minimum dwell time of 0.1 s
    ub = 20 * ones(num_dwells, 1);   % Maximum dwell time of 20 s

    % Optimization options
    options = optimoptions('fmincon', 'Display', 'iter', 'MaxIter', 300, 'MaxFunEvals', inf, 'TolFun', 1e-6, 'TolX', 1e-6);

    % Run the optimization
    dwell_times_opt = fmincon(@(dwell_times) objective(dwell_times, xGrid), dwell_times_initial, [], [], [], [], lb, ub, @(dwell_times) constraints(dwell_times, xGrid), options);

    % Display results
    disp('Optimized dwell times:');
    disp(dwell_times_opt);
end