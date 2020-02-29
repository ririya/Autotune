addpath('D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\Code\Binning Code');
addpath('D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\Code\MIDIFeatureExtraction');
addpath('D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\Code\MIDIFeatureExtraction\miditoolbox');
addpath('D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\Code\MIDIFeatureExtraction\miditoolbox\private');
addpath('D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\Code\voicebox');

maxSeqLen = 10;

midiDirectory = 'D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\Database\jazz';

modelPath = 'D:\ASUDropbox\Dropbox (ASU)\ClassProjects\speech\FinalProject\jazzModel.mat';

mFeatures = MusicFeatures( midiDirectory );

feats_adjustedTime = mFeatures;
for i = 1:length(mFeatures)
    feats_adjustedTime{i}(:,1) = adjustTime(mFeatures{i}(:,1));
end
% mFeaturesAdjusted = adjustTime(mFeatures);

[dictionary, vmap] =  createMusicModel(feats_adjustedTime,maxSeqLen);

save(modelPath,'dictionary','vmap');


