addpath('C:\ClassProjects\speech\FinalProject\voicebox')

% yinpath = 'C:\ClassProjects\speech\FinalProject\yin';
% 
% addDirAndSubDirToPath(yinpath);

close all
clear
% filename = 'C:\ClassProjects\speech\sine300.wav';
filename = 'C:\ClassProjects\speech\Q7speech\Track 50.wav';
% filename = 'C:\ClassProjects\speech\Q7speech\Track 49.wav';

MaxFreq = 600; %800 Hz
MinFreq = 50; %50 Hz
MaxPeriod =  1/MinFreq;

N = 3;

yinResult = yin(filename)
% yin(filename)


figure

% yinPitch = yinResult.f0;
% yinPitch(isnan(yinPitch)) = 0;
% 
% plot(yinResult.timeAxis, yinPitch)

% figure
% hold on
% yinPitch = yinResult.best;
% indBad = find(isnan(yinPitch));
% % yinPitch(indBad) = 0;
% 
% plot(yinResult.timeAxis, yinPitch)
% % plot(yinResult.timeAxis(indBad), yinPitch(indBad) ,'.r')
% % hold off
% figure

yinPitch = yinResult.good;
% yinPitch(isnan(yinPitch)) = 0;

plot(yinResult.timeAxis, yinPitch)

[x,Fs] =audioread(filename);

x = x(1:10*Fs,1);

lenX = length(x);

[fx,tt]=fxrapt(x,Fs,'u');

lenFx = length(fx);

ratio = lenX/ (lenFx+1)

ratio/ Fs

% fx(isnan(fx)) = 0;

minfx = 0;

maxfx = 400;

time = 1:lenFx;

time = time*0.01;


figure
plot(fx)
ylim([minfx, maxfx])
hold on



soundsc(x,Fs)
for i=1:lenFx
h1 = line([i i], [minfx maxfx]); 

pause(0.0092) 

delete(h1);
refresh;

end







    
