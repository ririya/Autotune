% rootDir = 'D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\';
rootDir = 'E:\Dropbox (ASU)\Dropbox (ASU)\ClassProjects\speech\FinalProject\';

addpath([rootDir '\Code\Binning Code']);
addpath([rootDir '\Code\MIDIFeatureExtraction']);
addpath([rootDir '\Code\MIDIFeatureExtraction\miditoolbox']);
addpath([rootDir '\Code\MIDIFeatureExtraction\miditoolbox\private']);
addpath([rootDir '\Code\voicebox']);
addpath([rootDir '\Code\voice-changer-master\src']);
addpath([rootDir '\Code\voice-changer-master\src\private']);

%params for psola pitch and time modification
global config;
config.pitchScale           = (2^(4/12));	%pitch scale ratio
config.timeScale            = 0.5;	%time scale ratio
% config.timeScale            = 1;	%time scale ratio
config.resamplingScale      = 1;		%resampling ratio to do formant shifting
config.reconstruct          = 0;		%if true do low-band spectrum reconstruction
config.displayPitchMarks    = 0;		%if true display pitch mark results
config.playWavOut           = 1;		%if true send output waveform to speaker
config.cutOffFreq           = 900;	%cut of frequency for lowpass filter
config.BPM = 60;
 config.maxInitialPitchChange = 3;    
    config.maxTimeScale = 4;
    config.maxSeqLen = 10;
    config.initialPitch = 110;
    config.initialDur = 0.5;
    
speechDirectory = [rootDir '\Database\speech'];
fileList = dir([speechDirectory '\' '*.wav']);

speechDurDirectory = [rootDir '\Code\SpeechFeatureExtraction'];

outputDirectory = [rootDir  '\Database\results'];

% midiDirectory = 'D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\Database\jazz';
% 
% modelPath = 'D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\jazzModel.mat';

midiDirectory = [rootDir '\Database\pop'];

modelPath = [rootDir  '\popModel.mat'];

if ~exist(modelPath)
    
mFeatures = MusicFeatures( midiDirectory );

feats_adjustedTime = mFeatures;
for i = 1:length(mFeatures)
    feats_adjustedTime{i}(:,1) = adjustTime(mFeatures{i}(:,1));
end

[dictionary, vmap] =  createMusicModel(feats_adjustedTime,config.maxSeqLen);

save(modelPath,'dictionary','vmap');

else 
    
    load(modelPath);
    
end


for f = 1:1
% for f = 1:length(fileList)
    
    filename = [speechDirectory '\' fileList(f).name];
    
    outputFile = [outputDirectory '\' fileList(f).name];
    
    [WaveIn, fs] = audioread(filename);
    WaveIn = WaveIn(:,1);
    
    speechDurFile = [speechDirectory '\' strrep(fileList(f).name, '.wav', '.xlsx')];
    
    startEnd = xlsread(speechDurFile);
    startEnd = startEnd(:,2:3);
    nSyllables = size(startEnd,1);  
    
    startSamples = round(startEnd(:,1)*fs);
    endSamples = round(startEnd(:,2)*fs);
    
    durations = diff(startEnd,1,2);    

    %obtain initial time interval
    musicalIntervals= durations.*(config.BPM/60); %change to musical interval
    adjusted_MusicalIntervals = adjustTime(musicalIntervals); %Adjust to musical scale
%     iTimeChange = adjusted_MusicalIntervals(1);

    %obtain initial musical interval
    maxInitialTimeChange = durations(1)*config.maxTimeScale;
    minInitialTimeChange = durations(1)/config.maxTimeScale;
%     iPitchChange = -maxInitialChange-1 + unidrnd(2*maxInitialChange+1,1000);
%     unique(iPitchChange)
%     iPitchChange = -maxInitialChange-1 + unidrnd(2*maxInitialChange+1);
    
    indAllowedPitch = find(abs(vmap(2,:))<=config.maxInitialPitchChange);
    indAllowedMaxTime = find(vmap(1,:)<=maxInitialTimeChange);
    indAllowedMinTime = find( vmap(1,:)>=minInitialTimeChange);
    indAllowed = intersect(indAllowedPitch,indAllowedMaxTime);
    indAllowed = intersect(indAllowed,indAllowedMinTime);
    
    if isempty(indAllowed)
        iPitchChange = config.initialPitch;
        iTimeChange =  config.initialDur;
    else      
        R =unidrnd(length(indAllowed));
        iPitchChange = vmap(2,indAllowed(R));
        iTimeChange = vmap(1,indAllowed(R));
    end
       
    Xreal = [iTimeChange;iPitchChange]; %initial input
    for k = 2:nSyllables-1
        [Y,maxProb,bestSeqLen] = getNextNoteFreqAndTime(Xreal, dictionary,vmap,config.maxSeqLen,config.maxTimeScale,durations(k));
        Xreal = [Xreal, Y];
    end
    
    
    
 %obtain initial pitch    
PitchContour = getPitchContour(WaveIn,fs,config.cutOffFreq);
          
for s = 1:nSyllables
    
avgPitchSyllables(s) = mean(PitchContour(startSamples(s):endSamples(s)));
    
end
    
avgPitch1stSyllable = avgPitchSyllables(1);
%       avgPitch1stSyllable = mean(avgPitchSyllables);   %another option
    [MusicFreq, LetterNote]= adjustFreq(avgPitch1stSyllable);
% endSamples = 
MusicTime   = iTimeChange;

for p = 1:length(Xreal)
   MusicFreq =  [MusicFreq , MusicFreq(end)*2^(Xreal(2,p)/12)];
   MusicTime = [MusicTime, Xreal(1,p)*60/config.BPM];
end

[~, LetterNote]= adjustFreq(MusicFreq);


[waveOut,~,dur]= PitchMarkingMod2(WaveIn, PitchContour, fs,MusicTime,MusicFreq,startSamples,endSamples);			
    
audiowrite(outputFile, waveOut, fs);	
    
end




