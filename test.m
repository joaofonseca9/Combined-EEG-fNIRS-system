N_trials=10;
events_autodual=randCuedTrials(N_trials);
%% HELPER FUNCTIONS
function events=randCuedTrials(n)
    %Generate a vector where half is 0's and half is 1's
    %1=cued, 0=uncued
    isCued = zeros(n, 1);
    isCued(randperm(numel(isCued), round(n/2))) = 1;
    for j=1:n
        events.trial(j).cue=isCued(j);
    end
end

% To Play Back Sound
function [WAVstruct] = CreateWAVstruct(WAVfilename)
% This function creates a struct with the information from the wav-files.

    wav = WAVfilename;                                          
    WAVstruct = struct('wavedata',[],'fs',[],'nrChan',[]);      
    [WAVstruct.wavedata, WAVstruct.fs] = psychwavread(wav);     
    WAVstruct.wavedata = WAVstruct.wavedata';                   
    WAVstruct.nrChan = size(WAVstruct.wavedata,1);  
end