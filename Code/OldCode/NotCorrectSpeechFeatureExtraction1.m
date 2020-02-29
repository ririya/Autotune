function [ Cs, pairs, formant1, formant2 ] = SpeechFeatureExtraction1(directory)
%Note: This is not robust yet as we had to use a shortcut to know monosyllabic word number
%Input is directory containing speech files
%Output cells that contain the start and end of a syllable, the formant1,
%and the formant 2 of the syllables, along with those raw values

% fnames = dir('C:/Users/Matt Jackson/Documents/MATLAB/SHS/Project/Speech');
fnames = dir(directory);
numf = length(fnames);

for n = 3: numf %Had to start at 3 for some reason as there are hidden files called '...'
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
bin = Smoothed>.1;
[ pairs ] = bitparser( bin );
while length(pairs) > numWords
lens = pairs(:,2)-pairs(:,1);
[mins, index] = min(lens);
pairs(index,:) = [];
end
lens = pairs(:,2)-pairs(:,1);

formant1 = [];
formant2 = [];
for i = 1:length(lens)
   currentFrame = sample((pairs(i,1)): (pairs(i,2)));
   zer = find(currentFrame==0);
   currentFrame(zer) = 1*10^-15;
   
    A = lpc(currentFrame,10);
   B =ones(1,11);
   [h,w] = freqz(B,A);
   h = abs(h);
   f = w*Fs/(2*pi);
   
   [pks,locs] = findpeaks(h);
   j = locs(1);
   k = locs(2);
   
   formant1 = [formant1; f(j)];
   formant2 = [formant2; f(k)];

end
Mat = [pairs, formant1, formant2];
Cs{n} = Mat;
end;
end

