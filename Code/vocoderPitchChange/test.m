% [d,sr]=audioread('sf1_cln.wav'); 

addpath('C:\ClassProjects\speech\FinalProject\voice-changer-master\src');
addpath('C:\ClassProjects\speech\FinalProject\voice-changer-master\src\private');

sr 
sr = 16000;
% 1024 samples is about 60 ms at 16kHz, a good window 
% y=pvoc(d,.75,1024); 
% 
% % Compare original and resynthesis 
% sound(d,16000) 
% sound(y,16000)

% [d,sr]=audioread('clar.wav'); 

 [d,sr]=audioread('C:\ClassProjects\speech\FinalProject\Database\speech\cam1.wav'); 
 
%  PSOLA(d, sr, anaPm, timeScale, pitchScale)
 
ratio = 2.^(2/12);
 
e = pvoc(d, ratio); 
f = resample(e,10^4,round(ratio*10^4)); % NB: 0.8 = 4/5 
% soundsc(d(1:length(f))+f,sr)
% soundsc(d,sr)
soundsc(f,sr)


% y=pvoc(d,.75,1024); 

% Compare original and resynthesis 
% sound(d,16000) 
% sound(y,16000)
