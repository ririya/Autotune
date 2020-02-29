function [freqs_adjusted,notes] = getNotes(freqs)
%takes in a vector of frequency values and normalizes them to music freqs
%outputs adjusted frequencies as well as the corresponding notes

global config;
%instantiate outputs
n = length(freqs);
freqs_adjusted = zeros(n,1);
notes = {};

%import pitch table
pitches = importdata(config.pitchChartPath);
note_freqs = pitches.data;
note_list = pitches.textdata;


%cycle through each input
for i = 1:n
    
    %check to see if the input is close to a freq value in the pitch table
    for j = 1:length(note_freqs)
        %if it's within a certain small range, adjust
        upper = note_freqs(j)*(2^(1/24));
        lower = note_freqs(j)*(2^(-1/24));
        if freqs(i) < upper && freqs(i) > lower
            freqs(i) = note_freqs(j);
        end
    end
    
    %reassign
    freqs_adjusted(i) = freqs(i);
    
    %get the note from the note table
    notes{i} = note_list{find(note_freqs==freqs_adjusted(i))};
end

end