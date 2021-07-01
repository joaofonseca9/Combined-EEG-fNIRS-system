function [power_theta, power_alpha, power_beta, power_gamma, freq_theta,...
    freq_alpha, freq_beta, freq_gamma] =...
    calculateAveragePowerBandsAllTrials(EEG, event_samp, startTask, endTask)
% Loop through the power from the individual trials and average them.

for trial=1:length(startTask)
    
        title = char(strcat('Trial_', string(trial)));
        startTask_times = event_samp(startTask(trial));
        endTask_times = event_samp(endTask(trial));
        
        if startTask_times < endTask_times
        EEG_trial = pop_select(EEG, 'point', [startTask_times endTask_times]);
        trial_data = EEG_trial.data;
        else
        EEG_trial = pop_select(EEG, 'point', [endTask_times startTask_times]);
        trial_data = EEG_trial.data;
        end

        [power_theta_oneTrial, power_alpha_oneTrial, power_beta_oneTrial,...
            power_gamma_oneTrial, freq_theta, freq_alpha, freq_beta,...
            freq_gamma] = calculatePowerBandsPerTrial(EEG_trial, trial_data);
            
        power_theta_allTrials(:, trial) = power_theta_oneTrial;
        power_alpha_allTrials(:, trial) = power_alpha_oneTrial;
        power_beta_allTrials(:, trial) = power_beta_oneTrial;
        power_gamma_allTrials(:, trial) = power_gamma_oneTrial;
    
    end
    
    power_theta = mean(power_theta_allTrials, 2);
    power_alpha = mean(power_alpha_allTrials, 2);
    power_beta = mean(power_beta_allTrials, 2);
    power_gamma = mean(power_gamma_allTrials, 2);

end