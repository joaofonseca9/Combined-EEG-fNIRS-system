function [time,type]=StimuliIn(type,outlet)
%lib = lsl_loadlib();
%fs=0;
%id='sdfwerr32432';
%info = lsl_streaminfo(lib,'EOEC','Markers',1,fs,'cf_int32',id);
% info2=lsl_streaminfo(lib,'TimeAudiIso','Markers',1,fs,'cf_int32',id);
%outlet = lsl_outlet(info); % thing you push your data through/ This is the
%marker associated with the task
outlet.push_sample(type)
% outlet.push_sample(GetSecs)
time=GetSecs;
% outlet.delete()
end 