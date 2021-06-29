function power_array_out = compensateRemovedChannels(power_array_in, EEG, list_channels, sub)
% Add NaN in the lines where channels were removed during pre-processing.

% Array to see which channels are missing.
list_present = zeros(30, 1);

% If there are less than 30 channels.
if size(power_array_in, 1)~=30
    
    % Initialize new power array.
    power_array_out = zeros(30, 1);
    power_array_out(1:size(power_array_in, 1), 1) = power_array_in;
    % Get the channels present in the signal.
    channels_present = EEG.chanlocs;
    % Go through the lists of channels supposed to be present and channels
    % actually present and mark them as 1 if present and as 0 if not
    % present.
    for i=1:size(list_channels, 1)
        for j=1:size(channels_present, 2)
            if convertCharsToStrings(channels_present(j).labels) == list_channels(i)
                list_present(i) = 1;
                break;
            end
        end
    end
    
    numMissing=0;
    % For the non-present channels, put NaN value in them and update the
    % rest of the values.
    for k=1:size(list_present, 1)
        if list_present(k)==0
            numMissing = numMissing+1;
            power_array_out(k) = NaN;
            for x=k+1:size(power_array_in, 1)
                power_array_out(x)=power_array_in(x-numMissing);
            end
        end
    end
    power_array_out(30:-1:30-numMissing+1, 1) =...
        power_array_in(size(power_array_in,1):-1:size(power_array_in, 1)-numMissing+1);
else
    power_array_out = power_array_in;
end
if sub == "64"
    power_array_out(1)=[];
end
end