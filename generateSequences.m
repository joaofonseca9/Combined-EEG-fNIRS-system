numbers = [1, 2, 3, 4];
prevElement = 0;
count = zeros(1, 4);
seq = zeros(1, 12);

for i=1:12
    
    while true
        newElement = randsample(numbers, 1);
        if newElement ~= prevElement
            if count(newElement)~=3
                seq(i) = newElement;
                prevElement = newElement;
                count(newElement) = count(newElement)+1;
                break;
            end
        end
    end
end

seq

% 2 4 3 4 1 3 4 1 2 1 3 2
% 4 2 4 1 3 1 4 2 3 2 1 3

