clear;
directory = 'C:/Users/Matt Jackson/Documents/MATLAB/SHS/Project/Speech';
fnames = dir(directory);
numf = length(fnames);

for n = 3: numf %Had to start at 3 for some reason as there are 2 hidden files called '...'
    fn = fnames(n).name;
    
    cmp1 = strcmp(fn,'cam1.wav') || strcmp(fn,'matt1.wav') || strcmp(fn,'rafael1.wav');
    cmp2 = strcmp(fn,'cam1.wav') || strcmp(fn,'matt2.wav') || strcmp(fn,'rafael2.wav');
    cmp3 = strcmp(fn,'cam3.wav') || strcmp(fn,'matt3.wav') || strcmp(fn,'rafael3.wav');
    cmp4 = strcmp(fn,'cam4.wav') || strcmp(fn,'matt4.wav') || strcmp(fn,'rafael4.wav');
    cmp5 = strcmp(fn,'cam5.wav') || strcmp(fn,'matt5.wav') || strcmp(fn,'rafael5.wav');
    cmp6 = strcmp(fn,'cam6.wav') || strcmp(fn,'matt6.wav') || strcmp(fn,'rafael6.wav');
    
    if cmp1 == 1
        numWords = 27;
    elseif cmp2 ==1
        numWords = 29;
    elseif cmp3 ==1
        numWords = 30;
    elseif cmp4 ==1
        numWords = 28;
    elseif cmp5 ==1
        numWords = 20;
    else
        numWords = 20;
    end;
    
    [voice, Fs] = audioread(fn);
    sample = voice(:,1);
    Energy = (abs(sample)).^2;
    E = Energy/sum(Energy);
    filt = ones(1, 200);
    Smoothed = 1000*conv(E, filt);
    %figure
    %plot(Smoothed)
    bin = Smoothed>.2;
    [ pairs ] = bitparser( bin );
    pairstemp = pairs;
    for p = 1:length(pairstemp)-1
        for m = 1:length(pairstemp)-1
            if p+m <= length(pairstemp)
                if (abs(pairstemp(p+m,1)-pairstemp(p,2)) <  2000 || abs(pairs(p+m,1)-pairs(p,2))<2000)&& pairs(p+1,1)-pairs(p,2)<2000 && pairstemp(p,1)~=0
                    pairstemp(p,2) = pairstemp(p+m,2);  
                    pairstemp(p+m,:)=0;
                end
            end
        end
    end
    tosspairstemp = (pairstemp(:,1)==0);
    pairstemp(tosspairstemp,:) = [];
    lens = pairstemp(:,2)-pairstemp(:,1);
    
    
    % while length(pairs) > numWords
    startsAndEnds = pairstemp;
    StartsAndEnds = [];
    formant1= [];
    durations = [];
    f0s = [];
    for p = 1: length(pairstemp)
        len = pairstemp(p,2)-pairstemp(p,1);
        %     lens = pairs(:,2)-pairs(:,1);
        %     [mins, index] = min(lens);
        currentFrame = sample((pairstemp(p,1)): (pairstemp(p,2)));
        if length(currentFrame) > 2000 %100
            zer = find(currentFrame==0);
            currentFrame(zer) = 1*10^-15;
            
            A = lpc(currentFrame,10);
            e = filter(A,1,currentFrame);
            B = ones(1,11);
            [h,w] = freqz(B,A);
            h = abs(h);
            f = w*Fs/(2*pi);
            
            [pks,locs] = findpeaks(h);
            j = locs(1);
            f1 = f(j);
%            if f1 <150 || f1>1000
%                startsAndEnds(p,:) = 0;
%            else
                [f0_time,f0,SHR,f0_candidates]=shrp(currentFrame,Fs,[50 550],1000*length(currentFrame)/Fs);
%                 if f0 > 0 && f0<1000
                    formant1=[formant1; f1];
                    f0s = [f0s; f0];
                    duration = pairstemp(p,2)-pairstemp(p,1);
                    durations = [durations; duration];
                    StartsAndEnds = [StartsAndEnds; pairstemp(p,:)];
%                 else
%                     startsAndEnds(p,:) = 0;
%                 end
   
%             end
        end
    end
    distance = zeros((length(StartsAndEnds)-1),1);
    for l = 1:(length(StartsAndEnds)-1)
        distance(l) = StartsAndEnds(l+1,1)-StartsAndEnds(l,2);
    end
    newdurations = durations;
    ftemp = formant1;
    for fval = 1:length(ftemp)-1
        for m = 1:length(ftemp)-1
            if fval+m <= length(ftemp)
                if ftemp(fval) == ftemp(fval+m) && ftemp(fval) == ftemp(fval+1) && newdurations(fval)~=0
                    newdurations(fval) = newdurations(fval)+newdurations(fval+m);
                    newdurations(fval+m)=0;
                end
            end
        end
    end
%     while length(durations) > numWords
%        loc = find(durations==min(durations));
%        durations(loc) = [];
%        f0s(loc) = [];
%     end
    loc = find(newdurations==0);
    newdurations(loc)=[];
    StartsAndEnds(loc,:)=[];
    f0s(loc) = [];
    Mat = [StartsAndEnds, newdurations, f0s];
    Cs{n} = Mat;
end
