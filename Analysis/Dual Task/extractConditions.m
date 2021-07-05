function [nirs_autocued, nirs_nonautocued, nirs_autouncued,...
    nirs_nonautouncued] = extractConditions(conditions, nirs_HbO2)
% Extract the different conditions to analyse
for con = 1:length(conditions)
    if con==1
        nirs_autocued = nirs_HbO2{con};
    elseif con==2
        nirs_nonautocued = nirs_HbO2{con};
    elseif con==3
        nirs_autouncued = nirs_HbO2{con};
    elseif con==4
        nirs_nonautouncued = nirs_HbO2{con};
    end
end
end