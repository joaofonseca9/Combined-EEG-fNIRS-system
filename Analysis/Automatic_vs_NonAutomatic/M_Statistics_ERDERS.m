%% Statistic for ERD/ERS.

clear; clc; close all;

results_path = 'C:\Users\maria\OneDrive\Ambiente de Trabalho\Automaticity Results\Topoplots';

load(fullfile(results_path, 'erders_allsubs.mat'), 'allsubs');

autouncued_ERD_ERS_theta = allsubs.avg.autouncued_ERD_ERS_theta;
autocued_ERD_ERS_theta = allsubs.avg.autocued_ERD_ERS_theta;

% Anderson-Darling test.
% Returns a test decision for the null hypothesis that the data in vector
% x is from a population with a normal distribution. 
% The alternative hypothesis is that x is not from a population with a
% normal distribution.
% The result h is 1 if the test rejects the null hypothesis at the 5%
% significance level, or 0 otherwise.
[h, p] = adtest(autouncued_ERD_ERS_theta);
[h, p] = adtest(autocued_ERD_ERS_theta);