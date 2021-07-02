%% Analysis of the Behavioural Data 
clear; clc; close all;

%% Initialize data
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis');
addpath('C:\Users\catar\OneDrive - Universidade do Porto\Twente\Combined-EEG-fNIRS-system\Analysis\Dual Task');
laptop = 'laptopCatarina';
[mainpath_in, mainpath_out, eeglab_path] = addFolders(laptop);
results_path = 'C:\Users\catar\OneDrive - Universidade do Porto\Twente\Data Analysis\behavioural';

subrec = ["28" "04" "A"; "64" "01" "A"; "02" "02" "B"; "76" "01" "B"];

A = '243413412132';
B = '413241423213';

numRemovedSubs = 0;

array_dual_finalAverageMistakes_cued = zeros(1, size(subrec, 1));
array_dual_finalAverageMistakes_uncued = zeros(1, size(subrec, 1));
array_dual_finalIncorrectSequences_cued = zeros(1, size(subrec, 1));
array_dual_finalIncorrectSequences_uncued = zeros(1, size(subrec, 1));
array_single_finalIncorrectSequences_cued = zeros(1, size(subrec, 1));
array_single_finalIncorrectSequences_uncued = zeros(1, size(subrec, 1));
array_dual_finalDelay_cued = zeros(1, size(subrec, 1));
array_dual_finalDelay_uncued = zeros(1, size(subrec, 1));
array_single_finalDelay_cued = zeros(1, size(subrec, 1));
array_single_finalDelay_uncued = zeros(1, size(subrec, 1));

%% Go through every subject
for subject = 1:size(subrec, 1)
    sub = subrec(subject, 1);
    rec = subrec(subject, 2);
    
    if subrec(subject, 3) == "A"
        sequence = B;
    else
        sequence = A;
    end
    
    subvar = genvarname(sub);

    % Load the subject's results
    load([mainpath_in, '\source\sub-', char(sub), '\stim\results_sub-',...
        char(sub), '_rec-', char(rec), '.mat']);
    
    % Make the changes by hand
    [events_dual, events_single] =...
        byHandChanges(sub, events_nonautodual, events_nonautosingle);
    
    %% Dual Task
    % Check if counting answers from sequence were correct
    dual_correct = checkCorrectCountingPerTrial(events_dual);

    % Check average number of mistakes per condition (cued and uncued)
    [dual_averageMistakes_cued, dual_averageMistakes_uncued] = ...
        averageMistakesPerCondition(events_dual, dual_correct);
    
    % Check number of incorrectly typed sequences per number of trials
    [dual_incorrectSequences_cued, dual_incorrectSequences_uncued] = ...
        checkIncorrectSequencePerTrial(events_dual, sequence);
       
    % Get the delay in the fingertapping performance
    [dual_delay_cued, dual_delay_uncued] = ...
        calculateDelayPerformance(events_dual, sequence);
    
    % Add values to array of all subjects   
    array_dual_finalAverageMistakes_cued(subject-numRemovedSubs) =...
        dual_averageMistakes_cued;
    array_dual_finalAverageMistakes_uncued(subject-numRemovedSubs) =...
        dual_averageMistakes_uncued;
    array_dual_finalIncorrectSequences_cued(subject-numRemovedSubs) =...
        dual_incorrectSequences_cued;
    array_dual_finalIncorrectSequences_uncued(subject-numRemovedSubs) =...
        dual_incorrectSequences_uncued;
    array_dual_finalDelay_cued(subject-numRemovedSubs) =...
        dual_delay_cued;
    array_dual_finalDelay_uncued(subject-numRemovedSubs) =...
        dual_delay_uncued;
    
    %% Single Task
    % Check number of incorrectly typed sequences per number of trials
    [single_incorrectSequences_cued, single_incorrectSequences_uncued] = ...
        checkIncorrectSequencePerTrial(events_single, sequence);
    
    % Get the delay in the fingertapping performance
    [single_delay_cued, single_delay_uncued] = ...
        calculateDelayPerformance(events_single, sequence);
    
    % Add values to array all subjects
    array_single_finalIncorrectSequences_cued(subject-numRemovedSubs) =...
        single_incorrectSequences_cued;
    array_single_finalIncorrectSequences_uncued(subject-numRemovedSubs) =...
        single_incorrectSequences_uncued;
    array_single_finalDelay_cued(subject-numRemovedSubs) =...
        single_delay_cued;
    array_single_finalDelay_uncued(subject-numRemovedSubs) =...
        single_delay_uncued;
    
    %% Put values of the current subject into final struct
    s.dual_avgMistakes_cued = dual_averageMistakes_cued;
    s.dual_avgMistakes_uncued = dual_averageMistakes_uncued;
    s.dual_incorrectSeqs_cued = dual_incorrectSequences_cued;
    s.dual_incorrectSeqs_uncued = dual_incorrectSequences_uncued;
    s.single_incorrectSeqs_cued = single_incorrectSequences_cued;
    s.single_incorrectSeqs_uncued = single_incorrectSequences_uncued;
    s.dual_delay_cued = dual_delay_cued;
    s.dual_delay_uncued = dual_delay_uncued;
    s.single_delay_cued = single_delay_cued;
    s.single_delay_uncued = single_delay_uncued;
    
    % Add struct of current subject to all subjects struct
    allsubs.(genvarname(strcat('sub', char(sub)))) = s;
      
end

%% Compare mistakes in dual cued/uncued and between motor/cognitive tasks
mean_dual_finalAverageMistakes_cued = mean(array_dual_finalAverageMistakes_cued);
std_dual_finalAverageMistakes_cued = std(array_dual_finalAverageMistakes_cued);
mean_dual_finalAverageMistakes_uncued = mean(array_dual_finalAverageMistakes_uncued);
std_dual_finalAverageMistakes_uncued = std(array_dual_finalAverageMistakes_uncued);
mean_dual_finalAverageMistakes = mean([array_dual_finalAverageMistakes_cued array_dual_finalAverageMistakes_uncued]);
std_dual_finalAverageMistakes = std([array_dual_finalAverageMistakes_cued array_dual_finalAverageMistakes_uncued]);

mean_dual_finalIncorrectSequences_cued = mean(array_dual_finalIncorrectSequences_cued);
std_dual_finalIncorrectSequences_cued = std(array_dual_finalIncorrectSequences_cued);
mean_dual_finalIncorrectSequences_uncued = mean(array_dual_finalIncorrectSequences_uncued);
std_dual_finalIncorrectSequences_uncued = std(array_dual_finalIncorrectSequences_uncued);
mean_dual_finalIncorrectSequences = mean([array_dual_finalIncorrectSequences_cued array_dual_finalIncorrectSequences_uncued]);
std_dual_finalIncorrectSequences = std([array_dual_finalIncorrectSequences_cued array_dual_finalIncorrectSequences_uncued]);

X1 = categorical({'Cognitive Task'; 'Motor Task'});
Y1 = [mean_dual_finalAverageMistakes; mean_dual_finalIncorrectSequences];
error1 = [std_dual_finalAverageMistakes; std_dual_finalIncorrectSequences];

figure;
b = bar(X1, Y1, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.3];
b.CData(2,:) = [0.5 1 0.5];
ylim([0 1]);
hold on;
errorbar(X1, Y1, error1, error1);    
hold off;

saveas(gcf, fullfile(results_path,...
    'Dual_LetterCountingvsFingerTappingMistakes'),'png');

X2 = categorical({'Dual Uncued'; 'Dual Cued'});
X2 = reordercats(X2, {'Dual Uncued'; 'Dual Cued'});
Y2 = [mean_dual_finalAverageMistakes_uncued; mean_dual_finalAverageMistakes_cued];
error2 = [std_dual_finalAverageMistakes_uncued; std_dual_finalAverageMistakes_cued];

X3 = categorical({'Dual Uncued'; 'Dual Cued'});
X3 = reordercats(X3, {'Dual Uncued'; 'Dual Cued'});
Y3 = [mean_dual_finalIncorrectSequences_uncued; mean_dual_finalIncorrectSequences_cued];
error3 = [std_dual_finalIncorrectSequences_uncued; std_dual_finalIncorrectSequences_cued];

figure;
subplot(1, 2, 1);
b = bar(X2, Y2, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.67 0.52];
b.CData(2,:) = [1 0.37 0.10];
ylim([0 1]);
hold on;
errorbar(X2, Y2, error2, error2);    
hold off;
subplot(1, 2, 2);
b = bar(X3, Y3, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [0.71 1 0.71];
b.CData(2,:) = [0.35 1 0.35];
ylim([0 1]);
hold on;
errorbar(X3, Y3, error3, error3);    
hold off;

saveas(gcf, fullfile(results_path,...
    'Dual_CuedvsUncued_LetterCountingvsFingerTappingMistakes'),'png');

%% Compare dual task performance with single task performance
mean_single_finalIncorrectSequences_cued = mean(array_single_finalIncorrectSequences_cued);
std_single_finalIncorrectSequences_cued = std(array_single_finalIncorrectSequences_cued);
mean_single_finalIncorrectSequences_uncued = mean(array_single_finalIncorrectSequences_uncued);
std_single_finalIncorrectSequences_uncued = std(array_single_finalIncorrectSequences_uncued);
mean_single_finalIncorrectSequences = mean([array_single_finalIncorrectSequences_cued array_single_finalIncorrectSequences_uncued]);
std_single_finalIncorrectSequences = std([array_single_finalIncorrectSequences_cued array_single_finalIncorrectSequences_uncued]);

X1 = categorical({'Dual'; 'Single'});
Y1 = [mean_dual_finalIncorrectSequences; mean_single_finalIncorrectSequences];
error1 = [std_dual_finalIncorrectSequences; std_single_finalIncorrectSequences];

figure;
b = bar(X1, Y1, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.3];
b.CData(2,:) = [0.5 1 0.5];
ylim([0 1]);
hold on;
errorbar(X1, Y1, error1, error1);    
hold off;

saveas(gcf, fullfile(results_path,...
    'DualvsSingle_FingerTappingMistakes'),'png');

X2 = categorical({'Dual Uncued'; 'Dual Cued'});
X2 = reordercats(X2, {'Dual Uncued'; 'Dual Cued'});
Y2 = [mean_dual_finalIncorrectSequences_uncued; mean_dual_finalIncorrectSequences_cued];
error2 = [std_dual_finalIncorrectSequences_uncued; std_dual_finalIncorrectSequences_cued];

X3 = categorical({'Single Uncued'; 'Single Cued'});
X3 = reordercats(X3, {'Single Uncued'; 'Single Cued'});
Y3 = [mean_single_finalIncorrectSequences_uncued; mean_single_finalIncorrectSequences_cued];
error3 = [std_single_finalIncorrectSequences_uncued; std_single_finalIncorrectSequences_cued];

figure;
subplot(1, 2, 1);
b = bar(X2, Y2, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.67 0.52];
b.CData(2,:) = [1 0.37 0.10];
ylim([0 1]);
hold on;
errorbar(X2, Y2, error2, error2);    
hold off;
subplot(1, 2, 2);
b = bar(X3, Y3, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [0.71 1 0.71];
b.CData(2,:) = [0.35 1 0.35];
ylim([0 1]);
hold on;
errorbar(X3, Y3, error3, error3);    
hold off;

saveas(gcf, fullfile(results_path,...
    'DualvsSingle_CuedvsUncued_FingerTappingMistakes'),'png');
 
%% Compare dual task performance with single task performance (delay)
mean_dual_finalDelay_cued = mean(array_dual_finalDelay_cued);
std_dual_finalDelay_cued = std(array_dual_finalDelay_cued);
mean_dual_finalDelay_uncued = mean(array_dual_finalDelay_uncued);
std_dual_finalDelay_uncued = std(array_dual_finalDelay_uncued);
mean_dual_finalDelay = mean([array_dual_finalDelay_cued array_dual_finalDelay_uncued]);
std_dual_finalDelay = std([array_dual_finalDelay_cued array_dual_finalDelay_uncued]);

mean_single_finalDelay_cued = mean(array_single_finalDelay_cued);
std_single_finalDelay_cued = std(array_single_finalDelay_cued);
mean_single_finalDelay_uncued = mean(array_single_finalDelay_uncued);
std_single_finalDelay_uncued = std(array_single_finalDelay_uncued);
mean_single_finalDelay = mean([array_single_finalDelay_cued array_single_finalDelay_uncued]);
std_single_finalDelay = std([array_single_finalDelay_cued array_single_finalDelay_uncued]);

X1 = categorical({'Dual'; 'Single'});
Y1 = [mean_dual_finalDelay; mean_single_finalDelay];
error1 = [std_dual_finalDelay; std_single_finalDelay];

figure;
b = bar(X1, Y1, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.5 0.3];
b.CData(2,:) = [0.5 1 0.5];
ylim([0 0.2]);
hold on;
errorbar(X1, Y1, error1, error1);    
hold off;

saveas(gcf, fullfile(results_path,...
    'DualvsSingle_FingerTappingDelay'),'png');

X2 = categorical({'Dual Uncued'; 'Dual Cued'});
X2 = reordercats(X2, {'Dual Uncued'; 'Dual Cued'});
Y2 = [mean_dual_finalDelay_uncued; mean_dual_finalDelay_cued];
error2 = [std_dual_finalDelay_uncued; std_dual_finalDelay_cued];

X3 = categorical({'Single Uncued'; 'Single Cued'});
X3 = reordercats(X3, {'Single Uncued'; 'Single Cued'});
Y3 = [mean_single_finalDelay_uncued; mean_single_finalDelay_cued];
error3 = [std_single_finalDelay_uncued; std_single_finalDelay_cued];

figure;
subplot(1, 2, 1);
b = bar(X2, Y2, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0.67 0.52];
b.CData(2,:) = [1 0.37 0.10];
ylim([0 0.2]);
hold on;
errorbar(X2, Y2, error2, error2);    
hold off;
subplot(1, 2, 2);
b = bar(X3, Y3, 0.4);
b.FaceColor = 'flat';
b.CData(1,:) = [0.71 1 0.71];
b.CData(2,:) = [0.35 1 0.35];
ylim([0 0.2]);
hold on;
errorbar(X3, Y3, error3, error3);    
hold off;

saveas(gcf,fullfile(results_path,...
    'DualvsSingle_CuedvsUncued_FingerTappingDelay'),'png');

%% Add average values of all subjects to final struct
average.dual_avgMistakes_cued = mean_dual_finalAverageMistakes_cued;
average.dual_avgMistakes_uncued = mean_dual_finalAverageMistakes_uncued;
average.dual_incorrectSeqs_cued = mean_dual_finalIncorrectSequences_cued;
average.dual_incorrectSeqs_uncued = mean_dual_finalIncorrectSequences_uncued;
average.single_incorrectSeqs_cued = mean_single_finalIncorrectSequences_cued;
average.single_incorrectSeqs_uncued = mean_single_finalIncorrectSequences_uncued;
average.dual_avgDelay_cued = mean_dual_finalDelay_cued;
average.dual_avgDelay_uncued = mean_dual_finalDelay_uncued;
average.single_avgDelay_cued = mean_single_finalDelay_cued;
average.single_avgDelay_uncued = mean_single_finalDelay_uncued;
allsubs.avg = average;

% Save the struct from all subs
save(strcat(results_path, '\allsubs.mat'), 'allsubs')

