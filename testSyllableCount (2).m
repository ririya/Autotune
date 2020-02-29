addpath('C:\ClassProjects\speech\FinalProject\voicebox');


speechDir = 'C:\ClassProjects\speech\FinalProject\Database\speech\';

resultsPath = 'C:\ClassProjects\speech\FinalProject\nSyllables10.mat';

fileList = dir([speechDir '*.wav']);

thresNanCount = 1;
maxThres = 200;

maxThres = 1;


errors = 10000*ones(maxThres,length(fileList));

minSumErrors = 10000;
minThres=1;

if exist(resultsPath)
   load( resultsPath);
end

correctNumberSyllables = [27 29 30 28 20 20];
correctNumberSyllables = repmat(correctNumberSyllables, 1, 3);

errors = [errors; 10000*ones(maxThres-size(errors,1),length(fileList))];
% minSumErrors = 10000;
% minThres=1;

for thresNanCount = 1:1:maxThres
    thresNanCount
    minThres
    minSumErrors
    if sum(abs(errors(thresNanCount,:))) ~= length(fileList)*10000;
        continue
    end
    
    nSyllables = zeros(1,length(fileList));
duration = {};



for i=1:length(fileList)
    
%     if nSyllables(i) > 0
%         continue
%     end
    

    fileName = [speechDir '\' fileList(i).name];
    
    [x,Fs] =audioread(fileName);
    x = x(:,1);    
    [fx,tt]=fxrapt(x,Fs,'u');

%     plot(fx)
    
% yinResult = yin(fileName);
% fx = yinResult.good';
% figure
% plot(fx)
    
    isNan = ~isnan(fx);
    
    while length(isNan)>0
         isNotNan = ~isNan;
        ind1stPitch = find(isNotNan,1,'first');   %looks for beginning of voiced sound
        
        if isempty(ind1stPitch)  %no more syllables          
            break;
        end
        
        isNan = isNan(ind1stPitch:end);
        
        indNanBegin = find(isNan);
        
        if isempty(indNanBegin)   % file ends during this syllable
           duration{i}.nSyllables(nSyllables(i)) =  length(isNan);
             nSyllables(i) = nSyllables(i)+1;
             break;
        else
            
            indLastNan = 0;
            
            noChanges = 1;
            
            for j=1:length(indNanBegin)
                
                ind1stNan = indNanBegin(j);                
                
%                 if indLastNan > ind1stNan
%                     continue
%                 end
                
                indNanEnd = find(~isNan(ind1stNan:end),1,'first');
                
                indLastNan = ind1stNan+indNanEnd-2;
                
                if ~isempty(indNanEnd)    %if found beginning of new syllable
                    
                    nanCount = indNanEnd - 1;
                    
                    if nanCount > thresNanCount   %if there are enough nans to make sure we are not breaking a syllable
                        
                        nSyllables(i) = nSyllables(i)+1;
                        
                         duration{i}.nSyllables(nSyllables(i)) = ind1stNan-1;
                        
                        isNan = isNan(indNanEnd:end);
                        
                        noChanges = 0;
                    end
                    
                else   %reached end of file (this is the last syllable)
                     nSyllables(i) = nSyllables(i)+1;      
                      duration{i}.nSyllables(nSyllables(i)) = ind1stNan-1;
                     isNan = [];
                     noChanges = 0
                     break;
                end
                
            end
            
            if (noChanges)
                  duration{i}.nSyllables(nSyllables(i)) =  length(isNan);
             nSyllables(i) = nSyllables(i)+1;
                          break
            end
            
        end
        
        
    end
    
%     save(resultsPath, 'nSyllables','duration');
    
        
end



errors(thresNanCount,:) = correctNumberSyllables - nSyllables;
sumErrors = sum(abs(errors),2);
[minSumErrors, minThres] = min(sumErrors);
save(resultsPath, 'errors', 'sumErrors','minSumErrors','minThres');

end
