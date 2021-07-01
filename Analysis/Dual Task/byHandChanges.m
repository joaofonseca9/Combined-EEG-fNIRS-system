function [events_dual, events_single] =...
    byHandChanges(sub, events_nonautodual, events_nonautosingle)
% For the subject passed as parameter, check if there are any by hand
% changes to be made from the lab notes and apply them.
  
    if sub == "02"
        events_nonautodual.trial(1) = [];
    end
    
    events_dual = events_nonautodual;
    events_single = events_nonautosingle;
              
end