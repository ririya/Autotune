%dir = 'C:\Users\Cameron\Documents\Academic\SHS 598\Database\pop';
%feats = MusicFeatures(dir);
feats_adjustedTime = feats;
for i = 1:length(feats)
    feats_adjustedTime{i}(:,1) = adjustTime(feats{i}(:,1));
end