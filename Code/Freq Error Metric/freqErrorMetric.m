function [percent_error,key] = freqErrorMetric(input_freqs, all_scales);
%takes in a vector of frequencies, and determines what key it's in. Then it
%determines how often the input frequencies stay in that key. Outputs a
%percent error, as a function of input vector length and error.



%generate a counter for each key, then increment that counter if a note is
%found in the input sequence that matches a note in the current key.
key_counter = zeros(12,1);

for i = 1:length(input_freqs)
    for j = 1:12
        for k = 1:length(all_scales{1,j})
            if input_freqs(i) == all_scales{1,j}(k)
                key_counter(j) = key_counter(j) + 1;
            end
        end
    end
    
    
    
end

%once the counter is complete, find the location of the maximum value.
% max_loc = find(key_counter == max(key_counter));
[maxValue, max_loc] = max(key_counter);

%the output key is the one where the counter had the highest value.
key = all_scales{2,max_loc};

%percent error is calculated by taking the amount of errors and dividing it
%by the total number of notes in the input sequence (aka the length)
percent_error = abs((max(key_counter)-length(input_freqs)))/length(input_freqs);

end