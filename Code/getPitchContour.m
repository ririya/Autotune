function PitchContour = getPitchContour(WaveIn,fs,cutOffFreq)

	%read input signal from file
WaveIn = WaveIn - mean(WaveIn); 				%normalize input wave

WaveIn = WaveIn(:,1);

[LowPass] = LowPassFilter(WaveIn, fs,cutOffFreq); %low-pass filter for pre-processing

WaveLength = length(LowPass);

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