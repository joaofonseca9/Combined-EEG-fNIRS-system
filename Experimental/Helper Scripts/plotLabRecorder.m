clear all; %#ok<CLALL>
[baseName, folder] = uigetfile();
fullFileName = fullfile(folder, baseName);
addpath('C:\Users\joaop\Downloads\xdf-Matlab-master');
% fullFileName='C:\Users\joaop\OneDrive - Universidade do Porto\Erasmus\Internship\Experiment\Combined-EEG-fNIRS-system\sub-P006\ses-S001\eeg\sub-P006_ses-S001_task-Default_run-001_eeg.xdf';
data=load_xdf(fullFileName);
time=data{1}.time_stamps;
load('markers.mat')
n=0;

%colors
colors=rand(length(markers.id),3); % Cell array of colros.

figure
for m=markers.id
    n=n+1;
    y=data{1}.time_series==m;
    time=nma_rescale(data{1}.time_stamps,0,max(data{1}.time_stamps)-min(data{1}.time_stamps));
    time(y==0)=[];
    y(y==0)=[];
    y=double(y);
    if ~isempty(time)
        if length(time)>1
            for ii=1:length(time)
                xline(time(ii), 'Color',colors(n,:))
            end
        else
            xline(time, 'Color',colors(n,:))
        end
    end
    h=text(time,y,markers.labels(n));
    set(h,'Rotation',90);
    hold on
end
axis([0 max(data{1}.time_stamps)-min(data{1}.time_stamps) 0.99985 1.0001])


%% Movie delay
% start_movie=data{1}.time_series==1798;
% start_auto_dual_no_cue=data{1}.time_series==1707;
% start_auto_dual_cue=data{1}.time_series==1706;
% 
% start_movie=data{1}.time_stamps(start_movie);
% start_auto_dual_no_cue=data{1}.time_stamps(start_auto_dual_no_cue);
% start_auto_dual_cue=data{1}.time_stamps(start_auto_dual_cue);
% 
% first_movie_delay=start_movie(1)-start_auto_dual_no_cue
% second_movie_delay=start_movie(2)-start_auto_dual_cue

%%
% legend('StartBlock_Metronome = 1698' , ...
% 'EndBlock_Metronome   = 1699' , ...
% 'StartBlock_Cue       = 1700' , ...         
% 'EndBlock_Cue         = 1701' , ...
% 'StartBlock_AutomaticSequence_Cued            = 1702' , ...
% 'StartBlock_AutomaticSequence_Uncued          = 1703' , ...
% 'StartBlock_NonAutomaticSequence_Cued         = 1704' , ...
% 'StartBlock_NonAutomaticSequence_Uncued       = 1705' , ...
% 'StartBlock_AutomaticSequence_Dual_Cued       = 1706' , ...
% 'StartBlock_AutomaticSequence_Dual_Uncued     = 1707' , ...
% 'StartBlock_NonAutomaticSequence_Dual_Cued    = 1708' , ...
% 'StartBlock_NonAutomaticSequence_Dual_Uncued  = 1709' , ...
% 'EndBlock_AutomaticSequence_Cued              = 1710' , ...
% 'EndBlock_AutomaticSequence_Uncued            = 1711' , ...
% 'EndBlock_NonAutomaticSequence_Cued           = 1712' , ...
% 'EndBlock_NonAutomaticSequence_Uncued         = 1713' , ...
% 'EndBlock_AutomaticSequence_Dual_Cued         = 1714' , ...
% 'EndBlock_AutomaticSequence_Dual_Uncued       = 1715' , ...
% 'EndBlock_NonAutomaticSequence_Dual_Cued      = 1716' , ...
% 'EndBlock_NonAutomaticSequence_Dual_Uncued    = 1717' , ...
% 'Keypress = 1777' , ...
% 'CHECK = 1255' , ...       
% 'start = 1555' , ...        
% 'stop = 1500', ...
% 'Marker_StopMovie = 1799', ...
% 'StartMovie = 1798', ...
% 'Marker_Letter = 1797') 
       
% 
% 'CHECK = 1255','        
% 'start = 1555','       
% 'stop = 1500')        
% scatter(time, data{1}.time_series==1700,  '.')
% hold on
% scatter(time, data{1}.time_series==1701,  '.')
% legend('Start Cue','End Cue')
% 
% scatter(time,data{1}.time_series==1777,  '.')
% legend('Keypress')
%  
% scatter(time,data{1}.time_series==1702,  '.')
% legend('Start-AutomaticSequence_Cued')
%  
% scatter(time,data{1}.time_series==1703,  '.')
% legend('Start-AutomaticSequence_Uncued')
%  
% scatter(time,data{1}.time_series==1704,  '.')
% legend('Start-NonAutomaticSequence_Cued')
%  
% scatter(time,data{1}.time_series==1705,  '.')
% legend('Start - NonAutomaticSequence_Uncued')
%  
% scatter(time,data{1}.time_series==1706,  '.')
% legend('Start - AutomaticSequence_Dual_Cued')
%  
% scatter(time,data{1}.time_series==1707,  '.')
% legend('Start - AutomaticSequence_Dual_Unued')
%  
% scatter(time,data{1}.time_series==1708,  '.')
% legend('Start - NonAutomaticSequence_Dual_Cued')
%  
% scatter(time,data{1}.time_series==1709,  '.')
% legend('Start - NonAutomaticSequence_Dual_Uncued')
%  
% scatter(time,data{1}.time_series==1710,  '.')
% legend('End - AutomaticSequence_Cued')
%  
% scatter(time,data{1}.time_series==1711,  '.')
% legend('End - AutomaticSequence_Uncued')
%  
% scatter(time,data{1}.time_series==1712,  '.')
% legend('End - NonAutomaticSequence_Cued')
%  
% scatter(time,data{1}.time_series==1713,  '.')
% legend('End - NonAutomaticSequence_Uncued')
%  
% scatter(time,data{1}.time_series==1715,  '.')
% legend('End - AutomaticSequence_Dual_Cued')
%  
% scatter(time,data{1}.time_series==1716,  '.')
% legend('End - AutomaticSequence_Dual_Uncued')
%  
% scatter(time,data{1}.time_series==1717,  '.')
% legend('End - NonAutomaticSequence_Dual_Cued')
%  
% scatter(time,data{1}.time_series==1718,  '.')
% legend('End - NonAutomaticSequence_Dual_Uncued')
%  
% scatter(data{1}.time_series==1255, time, '.')
% legend('Check flip')
%  
% scatter(time,data{1}.time_series==1555,  '.')
% legend('Start - Checkerboard')
%  
% scatter(time,data{1}.time_series==1500,  '.')
% legend('End - Checkerboard')
%  
% scatter( time,data{1}.time_series==1798, '.')
% legend('Start Movie')
%  
% scatter( time,data{1}.time_series==1799, '.')
% legend('End Movie')
%  
% scatter(time, data{1}.time_series==1797, '.')
% legend('Letter')
% hold off
% 


function C = nma_rescale(A,new_min,new_max)
%Nasser M. Abbasi 011212
%NO ERROR CHECKING DONE ON INPUT. Rescale a matrix or a vector A
current_max = max(A(:))';
current_min = min(A(:))';
C =((A-current_min)*(new_max-new_min))/(current_max-current_min) + new_min';
end