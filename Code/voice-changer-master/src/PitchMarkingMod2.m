%--------------------------------------------------------------------------
% Do pitch marking based-on pitch contour using dynamic programming
%--------------------------------------------------------------------------
function [waveOut,pm,dur] = PitchMarkingMod2(x, p, fs,targetTimeInterval,targetPitch,s,sEnd)

global config;
global data;

% split voiced / unvoiced segments
% [u, v] = UVSplit(p);

%  pitch marking for voiced segments
pm = [];
ca = [];
first = 1;
waveOut = [];


waveOut = [waveOut x(1:s(1))'];

pitchNonZero = p(p>0);
avgPitch = mean(pitchNonZero);

for i = 1 : length(s)
    i
    range = s(i) : sEnd(i);
    xSyllable = x(range);
    pSylabble = p(range);
    
    pSylabble(pSylabble == 0) = avgPitch;
    pSylabble(isinf(pSylabble)) = avgPitch;
    pSylabble(isnan(pSylabble)) = avgPitch;
%     for j = 2:length(pSylabble)
%        if pSylabble(j) == 0
% %            pSylabble(j) = pSylabble(j)-1;
%            pSylabble(j) = 100;
%        end
%     end
    
%     [u, v] = UVSplit(pSylabble);
    
%     if u(1) < v(1)      
% %         waveOut = [waveOut x(s(i)+u(1):s(i)+v(1)-1)'];       
%            if i==1
%           timeScaleUnvoiced = 1;
%            else
%           timeScaleUnvoiced = modification(i-1,2);
%            end
% %           waveOut =  [waveOut UnvoicedMod(x(s(i)+u(1):s(i)+v(1)-1), fs, timeScaleUnvoiced)'];
%             
%     end
  
%     in = xSyllable(range);
%     prange = pSylabble(range);

in = xSyllable;
 prange = pSylabble ;
 
    ca = [];

    
    marks = [];
    diffPRange = diff(prange);
    pitchChangeInd = find(abs(diffPRange) > 0.1);
    
    k=1;
    
    marks = k;
    
    while k < length(prange)
        
         marks = [marks k];
%         
        currentF0 = prange(k);
        currentPeriod = round((1/currentF0)*fs);
%         
        k = marks(end) + currentPeriod;
         
    end
    
    %marks(length(marks)) is the end of pitch cycle (relative index, for abs index sum range(1))
    %marks(1) is the beginning of pitch cycle   ((relative index), for abs index sum range(1))


pitchScale =  targetPitch(i)./prange;

%     timeScale = modification(i,2);
timeScale = targetTimeInterval(i)/(length(range)/fs);
   
   timeScale = timeScale*ones(length(pitchScale),1);

      synWav =  PSOLA2(in, fs, marks,timeScale,pitchScale);
     dur(i) = length(synWav)/fs;
    waveOut = [waveOut synWav];    
       
%    first = s(i) + range(1)+ marks(end)+1;  %beginning of unvoiced region
%    

% last = sEnd(i);
% %    if i < length(s)
% %    last = s(i+1);
% % last = sEnd(i);
% %    else
% %        last = length(x);
% %    end
% 
%        
%        
%   ra = first:last;  %sets unvoiced region
%     %     waveOut = [waveOut UnvoicedMod(x(ra), fs, config.timeScale)'];
%     waveOut = [waveOut x(ra)'];   %no need to timescale the unvoiced part

%     waveOut =  [waveOut UnvoicedMod(x(ra), fs, modification(i,2))'];
    
%             waveOut = [waveOut PSOLA(in, fs, marks)];
%              soundsc(synWav,fs)
%            soundsc(waveOut,fs)

end

% for i = 1 : length(v(:,1))
%     i
%     range = (v(i, 1) : v(i, 2));
%     in = x(range);
%     prange = p(range);
%     [marks1, cans] = VoicedSegmentMarking(in, prange, fs);
%     
%     pm = [pm  (marks1 + range(1))];
%     ca = [ca;  (cans + range(1))];
%     
%     %     if find(marks == 0)
%     
%     marks = [];
%     %     marks(1) = 1;
%     diffPRange = diff(prange);
%     pitchChangeInd = find(abs(diffPRange) > 0.1);
%     
%     k = 1;
%     
%     while i < length(prange)
%         
%         marks = [marks i];
%         
%         currentF0 = prange(i);
%         currentPeriod = round((1/currentF0)*fs);
%         
%         i = marks(end) + currentPeriod;
%         
%     end
%     
%     %marks(length(marks)) is the end of pitch cycle (relative index, for abs index sum range(1))
%     %marks(1) is the beginning of pitch cycle   ((relative index), for abs index sum range(1))
%     ra = first:marks(1)+range(1)-1;  %sets unvoiced region
%     first = marks(length(marks))+range(1)+1;  %updates beginning of unvoiced region as last mark
%     %     waveOut = [waveOut UnvoicedMod(x(ra), fs, config.timeScale)'];
%     waveOut = [waveOut x(ra)'];   %no need to timescale the unvoiced part
%  
%    pitchScale = pitchScaleAll(range);
%    timeScale = timeScaleAll(range);
%    
%    timeScale = config.timeScale*ones(length(pitchScale),1);
% 
%        synWav =  PSOLA2(in, fs, marks,timeScale,pitchScale);
%     waveOut = [waveOut synWav];
% 
% %             waveOut = [waveOut PSOLA(in, fs, marks)];
% %              soundsc(synWav,fs)
% %            soundsc(waveOut,fs)
%    
% end

% data.waveOut = waveOut;
% data.pitchMarks = pm;
% data.candidates = ca;

