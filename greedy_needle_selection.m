function [selected_needles, optimized_dwell_times] = greedy_needle_selection(all_needles, target_dose_range, max_selected_needles, target_mask, rectum_mask, xGrid, yGrid, zGrid)
    selected_needles = [];
    best_needle = [];
    best_dwell_times = [];
    best_satisfaction = -inf;
    best_marginal_gain = -inf;
    metrics = [];
    satisfaction_history = [];
    weights =  [.7, .2, .1]; % hr_ctv || % sphere [.4, .3, .3];  
    
    remaining_needles = all_needles;

    current_satisfaction = 0; % Initialize satisfaction
    dose_improvement_threshold = 0.01;  % Small threshold to consider as "significant improvement"

    % STEP 1: Find the Best Initial Needle
    % Commented out to test predefining first needle 
    % for i = 1:length(remaining_needles)
    %     temp_needles = remaining_needles(i);  % Single needle evaluation
    % 
    %     % Optimize dwell times for this needle alone
    %     temp_dwell_times = optimize_dwell_times(temp_needles, target_dose_range, target_mask, xGrid);
    % 
    %     % Compute satisfaction for this needle alone
    %     temp_satisfaction  = dose_eval(temp_needles, temp_dwell_times, target_mask, target_dose_range, xGrid, yGrid, zGrid, weights);
    % 
    %     % Keep track of the best single needle
    %     if temp_satisfaction > best_satisfaction
    %         best_satisfaction = temp_satisfaction;
    %         best_needle = remaining_needles(i);
    %         best_dwell_times = temp_dwell_times;
    %     end
    % end
    
    % Select the best initial needle
    best_needle = all_needles(1); % Select first needle 
    best_dwell_times = optimize_dwell_times(best_needle, target_dose_range, target_mask, rectum_mask, xGrid);

    % Compute satisfaction for this needle
    best_satisfaction = dose_eval(best_needle, best_dwell_times, target_mask, rectum_mask, target_dose_range, xGrid, yGrid, zGrid, weights);
    
    selected_needles = [selected_needles, best_needle];
    current_satisfaction = best_satisfaction;
    

    % Remove the selected needle from the remaining list
    best_needle_index = find([remaining_needles.idx] == best_needle.idx);
    remaining_needles(best_needle_index) = [];

    while length(selected_needles) < max_selected_needles && ~isempty(remaining_needles)
        best_needle = [];
        best_marginal_gain = -inf;
        temp_marginal_gain = 0; % Initialize temp marginal gain
        %best_contirbution = [];
        %interation_contributions = zeros(length(remaining_needles), 3);
        

        for i = 1:length(remaining_needles)
            temp_needles = [selected_needles, remaining_needles(i)];
            

            % Optimize dwell times for the current subset
            temp_dwell_times = optimize_dwell_times(temp_needles, target_dose_range, target_mask, rectum_mask, xGrid);

            % Calculate satisfaction for the current subset of needles
            temp_satisfaction = dose_eval(temp_needles, temp_dwell_times, target_mask, rectum_mask, target_dose_range, xGrid, yGrid, zGrid, weights);
            %iteration_contributions(i, :) = metrics;

            % Calculate marginal gain for subsequent needles
            temp_marginal_gain = temp_satisfaction - current_satisfaction;

            % Log the marginal gain for debugging
            disp(['Needle ', num2str(i), ': Marginal gain = ', num2str(temp_marginal_gain)]);

            % Update best marginal gain if the current gain is better
            if temp_marginal_gain > best_marginal_gain
                best_marginal_gain = temp_marginal_gain;
                best_dwell_times = temp_dwell_times;
                best_needle = remaining_needles(i);
                %best_contribution = iteration_contributions(i, :);
            end
            
        end

        % store best contribution metrics for iterative weight updates 
        %satisfaction_history = [satisfaction_history, best_contribution];

        % % Update satisfaction weights by boosting the weakest metric
        % if size(satisfaction_history, 1) >= 5  % Start tuning after a few iterations
        %     avg_contributions = mean(satisfaction_history, 1); 
        %     [~, min_index] = min(avg_contributions);  % Find weakest contributor
        % 
        %     % Adjust weight by shifting from strongest metric
        %     [~, max_index] = max(avg_contributions);
        %     weights(min_index) = weights(min_index) + 0.05; % Increase weight for weak metric
        %     weights(max_index) = max(0.05, weights(max_index) - 0.05); % Decrease weight for strong metric + ensure weights don't go negative
        %     weights = weights / sum(weights); % Normalize
        % end

        % If the best marginal gain is above the threshold, select the needle
        if ~isempty(best_needle) && (best_marginal_gain > 0 || isempty(selected_needles))
            selected_needles = [selected_needles, best_needle];
            current_satisfaction = current_satisfaction + best_marginal_gain;

            % Remove the selected needle from the remaining list
            best_needle_index = find([remaining_needles.idx] == best_needle.idx);
            remaining_needles(best_needle_index) = [];

            % Calculate the dose distribution for the selected needles
            dose_distribution = dose_calc(selected_needles, optimize_dwell_times(selected_needles, target_dose_range, target_mask, rectum_mask, xGrid), xGrid);

            dose2targ = dose_distribution(target_mask == 1);
            
            % Check if all target voxels have met the minimum dose requirement
            if sum(dose2targ(:) >= target_dose_range(1)) == sum(target_mask(:))
                disp('All target voxels meet the minimum dose requirement.');
                break;
            end
        else
            disp('No significant improvement found. Stopping the selection process.');
            break; % Stop if no improvement is possible
        end
    end
    optimized_dwell_times = best_dwell_times;
end