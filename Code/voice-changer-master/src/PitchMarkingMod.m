%--------------------------------------------------------------------------
% Do pitch marking based-on pitch contour using dynamic programming
%--------------------------------------------------------------------------
function [waveOut,pm,dur] = PitchMarkingMod(x, p, fs)

global config;
global data;

% split voiced / unvoiced segments
[u, v] = UVSplit(p);

%  pitch marking for voiced segments
pm = [];
ca = [];
first = 1;
waveOut = [];

%   s(1) = 1.177*fs;
%     s(2) = 1.557*fs;
%     s(3) = 2.104*fs;
%     s(4) = 2.474*fs;
%     s(5) = 2.993*fs;
%     s(6) = 3.308*fs;
%     s(7) = 4.466*fs;
%     s(8) = 4.698*fs;
%     s(9) = 5.319 *fs;
%     s(10) = 5.514 *fs;
%     s(11) = 5.727*fs;
%     s(12) = 6.422*fs;
%     s(13)= 6.514*fs;
%     s(14) = 7.469*fs;
%     s(15) = 7.802*fs;
%     s(16) = 8.312*fs;
%     s(17) = 8.877*fs;
%     s(18) = 9.109*fs;
%     s(19) = 9.359*fs;
%     
%     s(20) = 10.314 *fs;
%     s(21)= 10.582*fs;
%     s(22) = 10.777*fs;
%      s(23) = 11.240*fs;
%      
%       s(24) = 11.917*fs;
%        s(25) = 12.158 *fs;
%         s(26) = 12.445*fs;
%         s(27) = 12.936*fs;
% %          s14 = 1*fs;
% 
%     sSeconds = s/fs;
    
%     sSecondsDiff = diff(sSeconds);

xls = xlsread('C:\ClassProjects\speech\FinalProject\Database\speech\cam1.xlsx');

sSeconds=xls(:,2);
sEndSeconds = xls(:,3);

s = sSeconds*fs;
sEnd = sEndSeconds*fs;


    sSecondsDiff = sEndSeconds - sSeconds;
    BPM = 60;
    musicalTime = 0.25:0.25:4;
    musicalTimeSeconds = musicalTime*60/BPM;
    
    distances = (repmat(sSecondsDiff,1, length(musicalTimeSeconds)) - repmat(musicalTimeSeconds,length(sSecondsDiff),1)).^2;
    
    for i=1:length(sSecondsDiff)
       [minValue ind] = min(distances(i,:)); 
        minDist(i) = ind;
%         targetDist(1,i) = musicalTimeSeconds(ind);
        targetDist(1,i) = 0.5;
        targetDist(2,i) = sSecondsDiff(i);
        timeScaleSyllables(i) = targetDist(1,i)/targetDist(2,i);
    end
    targetDist(1,end+1) = 0.5;
    timeScaleSyllables(end+1) = 1;

    s = round(s);
    
    pitchScaleAll = ones(length(p),1);
%     timeScaleAll = ones(length(p),1);
    
    Abase = 110;
    
    notes = [Abase Abase*2^(2/12) Abase*2^(4/12) Abase*2^(5/12) Abase*2^(7/12) Abase*2^(9/12) Abase*2^(11/12) Abase*2^(11/12)]; 
    
    possibleTimeScales = [0.25 0.5 1 2 1 0.5 0.25]*2;
    
%     possibleTimeScales = circshift(possibleTimeScales,-1);
    
%      possibleTimeScales = [0.5 1 2];
    
     possibleTimeScales = 1;
    
    ascend = 1;
    
    noteInd = 1;
    timeInd = 1;

   


   for i = 1:length(s)

    note = notes(noteInd); 
    
%     modification(i,:) = [note possibleTimeScales(timeInd)];
      modification(i,:) = [note timeScaleSyllables(i)];  
    
    pitchScaleAll(s(i):sEnd(i)) = note./p(s(i):sEnd(i));
%      timeScaleAll(s(i):s(i+1)) = possibleTimeScales(timeInd);
%      timeInd = timeInd +1;
     if(timeInd > length(possibleTimeScales))
         timeInd = 1;
     end
     
    if (ascend)
    noteInd = noteInd+1;
    else
       noteInd = noteInd-1; 
    end
        if(noteInd > 7)
            ascend = 0;
            noteInd = 6; 
        end
    
     if(noteInd < 1)
        ascend = 1;
        noteInd = 1; 
    end
   
        
   end
         
         pitchScaleAll(isnan(pitchScaleAll)) = 1;
    pitchScaleAll(isinf(pitchScaleAll)) = 1;

    
waveOut = [waveOut x(1:s(1))'];

for i = 1 : length(s)
    i
    range = s(i) : sEnd(i);
    xSyllable = x(range);
    pSylabble = p(range);
    
    [u, v] = UVSplit(pSylabble);
    
    curSyn = [];
    
    if u(1) < v(1)      
%         waveOut = [waveOut x(s(i)+u(1):s(i)+v(1)-1)'];       
           if i==1
          timeScaleUnvoiced = 1;
           else
          timeScaleUnvoiced = modification(i-1,2);
           end
           
           unvoicedSyn = UnvoicedMod(x(s(i)+u(1):s(i)+v(1)-1), fs, timeScaleUnvoiced)';   %doesnt seem to work
%           waveOut =  [waveOut unvoicedSyn];          
          curSyn = [curSyn unvoicedSyn];
            
    end
    
    nVoicedRegions = length(v(:,1));
    
    lengthVoiced = sum(diff(v,1));
    durVoiced = lengthVoiced/fs;
    
    for j = 1 : nVoicedRegions
    
    range = (v(j, 1) : v(j, 2));
    in = xSyllable(range);
    prange = pSylabble(range);
  
%     [marks1, cans] = VoicedSegmentMarking(in, prange, fs);
%     
%     pm = [pm  (marks1 + range(1))];
%     ca = [ca;  (cans + range(1))];

    ca = [];
    
    %     if find(marks == 0)
    
    marks = [];
    %     marks(1) = 1;
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

pitchScale =  modification(i,1)./prange;

    timeScale = modification(i,2);

%     timeScale = 0.9*targetDist(1,i)/durVoiced;

%     timeScale = targetDist(1,i)/(length(range)/fs);
%  timeScale = 1;
   
   timeScale = timeScale*ones(length(pitchScale),1);

      synWav =  PSOLA2(in, fs, marks,timeScale,pitchScale);      
      curSyn = [curSyn synWav];
      
%     waveOut = [waveOut synWav];
    
       
   first = s(i) + range(1)+ marks(end)+1;  %beginning of unvoiced region
   if (j < nVoicedRegions)
       last =  s(i) + v(j+1, 1)-1;     %end of unvoiced region
   else       
       last = sEnd(i);          %end of unvoiced region
%        if i < length(s)
%        last = s(i+1);
% %         last = sEnd(i);
%        else
%            last = length(x);
%        end
   end
       
       
  ra = first:last;  %sets unvoiced region
%     waveOut = [waveOut x(ra)'];   %no timescale on unvoiced part
    unvoicedSyn = UnvoicedMod(x(ra), fs, modification(i,2))';    %doesnt seem to work
  
    curSyn = [curSyn unvoicedSyn];
%             waveOut = [waveOut PSOLA(in, fs, marks)];
%              soundsc(synWav,fs)
%            soundsc(waveOut,fs)
   
    end
    
    dur(i) = length(curSyn)/fs;
    waveOut =  [waveOut curSyn];
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

data.waveOut = waveOut;
data.pitchMarks = pm;
data.candidates = ca;

