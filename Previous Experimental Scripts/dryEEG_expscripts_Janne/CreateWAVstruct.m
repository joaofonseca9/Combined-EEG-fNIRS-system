function [WAVstruct] = CreateWAVstruct (WAVfilename)
%% CreateWAVstruct
% > load a wav-file and create a struct containing the relevant information
% > WAVstruct
%   > wavedata = sound data 
%   > fs = sample frequency [Hz]
%   > nrChan = number of channels [#]


wav = WAVfilename;                                          % > load wav file
WAVstruct = struct('wavedata',[],'fs',[],'nrChan',[]);      % > create struct for wavedata
[WAVstruct.wavedata, WAVstruct.fs] = psychwavread(wav);     % > read wav file (fs, wavedata)
WAVstruct.wavedata = WAVstruct.wavedata';                   % > switch rows and columns
WAVstruct.nrChan = size(WAVstruct.wavedata,1);              % > indicate number of channels             
 
