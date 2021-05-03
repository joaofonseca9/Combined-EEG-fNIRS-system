data=load_xdf('C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Combined-EEG-fNIRS-system\sub-P001\ses-S001\eeg\sub-P001_ses-S001_task-Default_run-003_eeg.xdf');

time=nma_rescale(data{1}.time_stamps,0,max(data{1}.time_stamps)-min(data{1}.time_stamps));

scatter(time, data{1}.time_series==1700,  '.')
hold on
scatter(time, data{1}.time_series==1701,  '.')
legend('Start Cue','End Cue')

scatter(time,data{1}.time_series==1777,  '.')
legend('Keypress')
 
scatter(time,data{1}.time_series==1702,  '.')
legend('Start-AutomaticSequence_Cued')
 
scatter(time,data{1}.time_series==1703,  '.')
legend('Start-AutomaticSequence_Uncued')
 
scatter(time,data{1}.time_series==1704,  '.')
legend('Start-NonAutomaticSequence_Cued')
 
scatter(time,data{1}.time_series==1705,  '.')
legend('Start - NonAutomaticSequence_Uncued')
 
scatter(time,data{1}.time_series==1706,  '.')
legend('Start - AutomaticSequence_Dual_Cued')
 
scatter(time,data{1}.time_series==1707,  '.')
legend('Start - AutomaticSequence_Dual_Unued')
 
scatter(time,data{1}.time_series==1708,  '.')
legend('Start - NonAutomaticSequence_Dual_Cued')
 
scatter(time,data{1}.time_series==1709,  '.')
legend('Start - NonAutomaticSequence_Dual_Uncued')
 
scatter(time,data{1}.time_series==1710,  '.')
legend('End - AutomaticSequence_Cued')
 
scatter(time,data{1}.time_series==1711,  '.')
legend('End - AutomaticSequence_Uncued')
 
scatter(time,data{1}.time_series==1712,  '.')
legend('End - NonAutomaticSequence_Cued')
 
scatter(time,data{1}.time_series==1713,  '.')
legend('End - NonAutomaticSequence_Uncued')
 
scatter(time,data{1}.time_series==1715,  '.')
legend('End - AutomaticSequence_Dual_Cued')
 
scatter(time,data{1}.time_series==1716,  '.')
legend('End - AutomaticSequence_Dual_Uncued')
 
scatter(time,data{1}.time_series==1717,  '.')
legend('End - NonAutomaticSequence_Dual_Cued')
 
scatter(time,data{1}.time_series==1718,  '.')
legend('End - NonAutomaticSequence_Dual_Uncued')
 
scatter(data{1}.time_series==1255, time, '.')
legend('Check flip')
 
scatter(time,data{1}.time_series==1555,  '.')
legend('Start - Checkerboard')
 
scatter(time,data{1}.time_series==1500,  '.')
legend('End - Checkerboard')
 
scatter( time,data{1}.time_series==1798, '.')
legend('Start Movie')
 
scatter( time,data{1}.time_series==1799, '.')
legend('End Movie')
 
scatter(time, data{1}.time_series==1797, '.')
legend('Letter')
hold off



function C = nma_rescale(A,new_min,new_max)
%Nasser M. Abbasi 011212
%NO ERROR CHECKING DONE ON INPUT. Rescale a matrix or a vector A
current_max = max(A(:));
current_min = min(A(:));
C =((A-current_min)*(new_max-new_min))/(current_max-current_min) + new_min;
end