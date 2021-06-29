function [Ploc, Ptime, Pamp] = get_mainVEPcomponents (time, avgGFP, window)

% find main VEP component
Pwin = find(time>=window(1,1) & time<=window(1,2));     % time window in which peak should be found
[Pamp,idx_loc] = max(avgGFP(Pwin));                     % find max peak in average GFPvalues within time window
Ploc = Pwin(idx_loc); clear idx_loc;                    % sample number of peak
Ptime = time(Ploc);                                     % time (s) of  peak

end