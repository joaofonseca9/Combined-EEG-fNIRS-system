function correctSequence = checkCorrectSequence(typed_sequence, real_sequence)
% Check if the sequence performed was equal to the real sequence.

if all(strcmp(typed_sequence, real_sequence))
    correctSequence = true;
else
    correctSequence = false;
end

end