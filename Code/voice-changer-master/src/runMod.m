% addpath('C:\ClassProjects\speech\FinalProject\voicebox');

addpath('D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\Code\voicebox');

%--------------------------------------------------------------------------
% main script to do pitch and time scale modification of speech signal
%--------------------------------------------------------------------------
% config contain all parameter of this program
global config;
% config.pitchScale           = 1/(2^(12/12));	%pitch scale ratio
config.pitchScale           = (2^(4/12));	%pitch scale ratio
 config.timeScale            = 0.5;	%time scale ratio
% config.timeScale            = 1;	%time scale ratio
config.resamplingScale      = 1;		%resampling ratio to do formant shifting
config.reconstruct          = 0;		%if true do low-band spectrum reconstruction
config.displayPitchMarks    = 0;		%if true display pitch mark results
config.playWavOut           = 1;		%if true send output waveform to speaker
config.cutOffFreq           = 900;	%cut of frequency for lowpass filter
% config.fileIn               = '..\waves\m2.wav';		%input file full path and name

% config.fileIn               = '..\waves\m2.wav';		%input file full path and name

% config.fileIn = 'C:\ClassProjects\speech\FinalProject\Database\speech\cam1.wav';

rootDir = 'D:\ASUDropbox\Dropbox (ASU)\ClassProjects\';

config.fileIn = [rootDir '\speech\FinalProject\Database\speech\cam1.wav'];
% config.fileIn = 'C:\ClassProjects\speech\FinalProject\vocoder\clar.wav';

% config.fileOut              = '..\waves\syn.wav';		%output file full path and name

config.fileOut   = [rootDir '\speech\FinalProject\Database\speech\modified_cam1.wav'];		%output file full path and name

directoryInput = [rootDir '\speech\FinalProject\Database\speech\'];
directoryOutput = [rootDir '\speech\FinalProject\Database\results\'];

if ~exist(directoryOutput)
    mkdir(directoryOutput)
end

filelist = dir([directoryInput '*.wav']);

%data contain analysis results
global data;
data.waveOut = [];		%waveform after do pitch and time scale modification
data.pitchMarks = [];	%pitch marks of input signal
data.Candidates = [];	%pitch marks candidates

% for file=1:length(filelist)
    
for file=1:1

    filename = [directoryInput filelist(file).name];
    
    outputFile = [directoryOutput strrep(filelist(file).name, '.wav', '_mod.wav')];

% [WaveIn, fs] = audioread(config.fileIn);	%read input signal from file
[WaveIn, fs] = audioread(filename);	%read input signal from file
WaveIn = WaveIn - mean(WaveIn); 				%normalize input wave

WaveIn = WaveIn(:,1);

[LowPass] = LowPassFilter(WaveIn, fs, config.cutOffFreq); %low-pass filter for pre-processing

% [PitchContour,FramePitch] = PitchEstimation(LowPass, fs);	

% % pitch contour estimation

WaveLength = length(LowPass);

WaveLength = length(WaveIn);

FramePitch=fxrapt(LowPass,fs,'u');

FrameRate = round(fs * 0.010);
 
PitchContour = zeros(WaveLength, 1);

lenPitch = length(FramePitch);

% Using median filter for Pos-processing
FramePitch = medfilt1(FramePitch, 5);

% calculate pitch contour

for i = 1 : WaveLength
    
    indx = floor(i / FrameRate) + 1;
    
    if indx <= lenPitch
    
        PitchContour(i) = FramePitch(indx);
    
    end
end
 PitchContour(isnan(PitchContour)) = 0;
 
 
 xls = xlsread([rootDir '\speech\FinalProject\Database\speech\cam1.xlsx']);

sSeconds=xls(:,2);
sEndSeconds = xls(:,3);

s = sSeconds*fs;
sEnd = sEndSeconds*fs;

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
%          s14 = 1*fs;

%     sSeconds = s/fs;
%     sSecondsDiff = diff(sSeconds);


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
    
    pitchScaleAll = ones(length(PitchContour),1);
%     timeScaleAll = ones(length(p),1);
    
    Abase = 110;
    
    notes = [Abase Abase*2^(2/12) Abase*2^(4/12) Abase*2^(5/12) Abase*2^(7/12) Abase*2^(9/12) Abase*2^(11/12) Abase*2^(11/12)]; 
          
    ascend = 1;
    
    noteInd = 1;
    timeInd = 1;


   for i = 1:length(s)

    note = notes(noteInd);

      modification(i,:) = [note timeScaleSyllables(i)];  
    
    pitchScaleAll(s(i):sEnd(i)) = note./PitchContour(s(i):sEnd(i));
     
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
 

[waveOut,~,dur]= PitchMarkingMod2(WaveIn, PitchContour, fs,targetDist(1,:),modification(:,1),s,sEnd);										%do pitch marking and PSOLA

% soundsc(WaveIn, fs)
% soundsc(data.waveOut, fs)



% audiowrite(config.fileOut, data.waveOut, fs);								%write output result to file
audiowrite(outputFile, waveOut, fs);	

end

% if config.playWavOut
%     wavplay(data.waveOut, fs);
% end
% 
% if config.displayPitchMarks
%     PlotPitchMarks(WaveIn, data.candidates, data.pitchMarks, PitchContour); %show the pitch marks
% end