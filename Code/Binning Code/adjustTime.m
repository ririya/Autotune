function [adjusted_times] = adjustTime(input_times)
%takes in a vector of time values and normalizes them to music times
%outputs adjusted times as a vector

n = length(input_times);
%instantiate new time array
% new_times = [];
% bound = 0.03125;
step = 0.0625;
% step = 2*0.0625;
bound = step/2;
minTime = step;
maxTime = 2;

adjusted_times = input_times;

%fill new time array
% for i = 1:64
%     new_times(i) = 0.0625*i;
% end

new_times = minTime:step:maxTime;

for i = 1:n
    
    %check to see if the input is close to a freq value in the pitch table
    for j = 1:length(new_times)
        %if it's within a certain small range, adjust
        upper = new_times(j) + bound;
        lower = new_times(j) - bound;
        if input_times(i) < upper && input_times(i) > lower
%             input_times(i) = new_times(j);
            adjusted_times(i) = new_times(j);
            break;
%         else if input_times(i) < new_times(1)
%                 input_times(i) = new_times(1);
%                  break;
%             else if input_times(i) >= new_times(end)
%                     input_times(i) = new_times(end);
%                      break;
%                 end
%             end
            
            %reassign
%             adjusted_times(i) = input_times(i);
            
        end
        
    end
    
end

adjusted_times(adjusted_times>new_times(end))=new_times(end);
adjusted_times(adjusted_times<new_times(1))=new_times(1);

