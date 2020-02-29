function [mFeatures] = MusicFeatures( directory )
%This function will read the database of MIDI files and output a mat file
%Input is directory (string)
%Output is 2 columns. Midi length and Midi pitch step

% fnames = dir('C:/Users/Matt Jackson/Documents/MATLAB/SHS/Project/Music');
fnames = dir(directory);
numf = length(fnames);

for i = 3: numf %Had to start at 3 for some reason as there are hidden files called '...'
    fn = [directory '\' fnames(i).name];

%From MIDI ToolBox
mstr = mdlMidiToMStr(fn); %Read MIDI file and store in a Matlab structure
nmat = mdlMStrToNMat(mstr);
%   1. Onset (beats) (calc'ed solely from ticks-per-quarter-note)
%   2. Duration beats
%   3. Channel
%   4. MIDI pitch
%   5. mIdI velocity
%   6. Onset (sec)
%   7. Duration (sec)
nmat = nmat(setdiff(1:size(nmat,1),find(nmat(:,2)==-1)),:); % remove notes with duration -1; PT 17.5.2016

features = nmat(:,1:4);
features(:,2) = []; %Remove the Duration
features(:,2) = [];%Remove the Channel which is now in column 2
der = diff(features); % getting the distance between starts of pitches (assuming no pauses) and MIDI pitch changes
%Octaves start at C = [0,12, 24, 36...120] 
%From C--> B is 11 steps
%Stops at MIDI 127 at G
c{i} = der;
mFeatures = c(~cellfun(@isempty, c));
end
end

