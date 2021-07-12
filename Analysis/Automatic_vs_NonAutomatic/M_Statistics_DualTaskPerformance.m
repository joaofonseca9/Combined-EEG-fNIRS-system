%% Statistics for performance under dual-task conditions.

clear; clc; close all;

results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Dual Task Performance';

subrec = ["28"; "02"; "76"];

% Load the struct from all subs.
load(strcat(results_path, '\allsubs.mat'))

% Run through every subject.
i=1;
for subject=1:size(subrec)
    sub = subrec(subject);
    s = allsubs.(genvarname(strcat('sub', char(sub))));

    %% Average mistakes on letter counting.
    auto_avgMistakes(i, 1) = s.autodual_avgMistakes_uncued;
    auto_avgMistakes(i, 2) = s.autodual_avgMistakes_cued;
    nonauto_avgMistakes(i, 1) = s.nonautodual_avgMistakes_uncued;
    nonauto_avgMistakes(i, 2) = s.nonautodual_avgMistakes_cued;
    
    %% Average mistakes on finger tapping.
    auto_incorrectSeqs(i, 1) = s.autodual_incorrectSeqs_uncued;
    auto_incorrectSeqs(i, 2) = s.autodual_incorrectSeqs_cued;
    nonauto_incorrectSeqs(i, 1) = s.nonautodual_incorrectSeqs_uncued;
    nonauto_incorrectSeqs(i, 2) = s.nonautodual_incorrectSeqs_cued;
    
    %% Delay on finger tapping.
    auto_delay(i, 1) = s.autodual_delay_uncued;
    auto_delay(i, 2) = s.autodual_delay_cued;
    nonauto_delay(i, 1) = s.nonautodual_delay_uncued;
    nonauto_delay(i, 2) = s.nonautodual_delay_cued;
    
    i = i+1;
end

%% Paired t-test.

% Average mistakes on letter counting.
[h, p] = ttest(auto_avgMistakes(:, 1), auto_avgMistakes(:, 2));
stats.lettercounting_auto = p;
[h, p] = ttest(nonauto_avgMistakes(:, 1), nonauto_avgMistakes(:, 2));
stats.lettercounting_nonauto = p;

% Average mistakes on finger tapping.
[h, p] = ttest(auto_incorrectSeqs(:, 1), auto_incorrectSeqs(:, 2));
stats.fingertapping_auto = p;
[h, p] = ttest(nonauto_incorrectSeqs(:, 1), nonauto_incorrectSeqs(:, 2));
stats.fingertapping_nonauto = p;

% Delay on finger tapping.
[h, p] = ttest(auto_delay(:, 1), auto_delay(:, 2));
stats.delay_auto = p;
[h, p] = ttest(nonauto_delay(:, 1), nonauto_delay(:, 2));
stats.delay_nonauto = p;