function array_boundaries = removeExtraBoundaries(array_boundaries,...
    array_task)
% Eliminate boundary points which are not followed by a start of task.
    
for i=length(array_boundaries):-1:1
    if ~any(array_boundaries(i)+1==array_task(:))
        array_boundaries(i) = [];
    end
end
end