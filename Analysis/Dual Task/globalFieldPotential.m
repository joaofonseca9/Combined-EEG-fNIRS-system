function [GFPt] = globalFieldPotential(ERP)
% calculate difference between two channels, per time point 
% take the sum of the differences, per time point
% gives the GFPt over time

for iTime = 1:size(ERP,2)
    Udiff = zeros(size(ERP,1), size(ERP,1));
    for i = 1:size(ERP,1)
        for j = 1:size(ERP,1)
            Udiff(i,j) = (ERP(i,iTime)-ERP(j,iTime)).^2;
        end
    end
    GFPt(1,iTime) = sqrt( (1/(2*size(ERP,1))) * sum(Udiff,'all') );
end
end