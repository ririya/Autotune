function [Output] = OutputChanges(BPM)
%Input is BPM (usually 60). Maybe want input of a filename and start/end of
%syllables
%Output is the sequence of freqencies in row 1 and BPM in row 2. A cell for
%each file for a given type of music. May want to add both types of music
%in here if we are going to get it all over with
clear all;


%We may want to do 1 speech file at a time do to how long this takes to run

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

%Not sure yet how to automize the start ad end time etraction from the
%excel file or if it should just be an input....

%%%%Input
%filename = 'cam5.wav';
[Speech, Fs] = audioread(fn);
Speech = Speech(:,1);

%Dont know if I should automize this or just add as an input
sylStart = [1.692,1.871,2.1, 2.289,2.605,2.867,3.366,3.591,4.031,4.27,5.29,5.612,6.046,6.555,7.348,7.6,7.9,8.24,8.6,8.768]; 
sylEnd =[1.849,2.047,2.256,2.599,2.835,3.312,3.564,3.993,4.25,4.711,5.574,6.019,6.523,7.107,7.584,7.842,8.238,8.554,8.715,8.98];



%%%%%From Start to Start
durations = sylEnd-sylStart; %find the length of each syllable in time
MusicalTime= durations.*(BPM/60); %change to BPM
adjusted_times = adjustTime(MusicalTime); %Adjust to musical scale
BPMTime = adjusted_times(1); %Initial duration
iTimeChange = adjusted_times(2)-adjusted_times(1); %Initial duration Change

%%%%%%Pitches
 [fx,tt]=fxrapt(Speech,Fs,'u'); %Solve for f0 values
 exists = fx(~isnan(fx)); %take the values that are not NaN
 avg = sum(exists)/length(exists); %find the average f0 value
[MusicFreq, LetterNote]= adjustFreq(avg); %Adjust to musical scale
%Randomly generate second note
iPitchChange = round(10*(rand(1)-.5)); %Intiail pitch change between -5 and 5 for pitch change
SecondNote= MusicFreq*2^(iPitchChange/12); %For checkin the random note

%%%%%Model
%Create dictionary and vmap
maxSeqLen = numWords-1; %How long the sequence can be max????????
midiDirectory = 'C:/Users/Matt Jackson/Documents/MATLAB/SHS/Project/Music/pop';
mFeatures = MusicFeatures( midiDirectory ); %Get the music features
feats_adjustedTime = mFeatures;
for i = 1:length(mFeatures)
    feats_adjustedTime{i}(:,1) = adjustTime(mFeatures{i}(:,1));
end
[dictionary, vmap] =  createMusicModel(feats_adjustedTime,maxSeqLen); %create model based off of Midi

%Generate Sequences
Xreal = [iPitchChange; iTimeChange]; %initial input
for k = 1:numWords-2 %only want numwords -1 changes total so need -2 here
Y = getNextNoteFreqAndTime(Xreal, dictionary,vmap);
Xreal = [Xreal, Y];
end

MusicFreq = [MusicFreq];
MusicTime = [BPMTime*60/BPM];
for p = 1:length(Xreal)
   MusicFreq =  [MusicFreq , MusicFreq(p)*2^(Xreal(1,p)/12)];
   MusicTime = [MusicTime, MusicTime(p)+Xreal(2,p)*60/BPM];
end

Values = [MusicFreq;MusicTime];
Cell{n} = Values;
Output = Cell(~cellfun(@isempty, Cell));

end

