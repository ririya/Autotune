function [scale] = generateScale(key)
%generates a list of frequencies belonging to a scale specificed with an
%input key (i.e. C,C#,D,etc.)

%initialize scale vector
scale = zeros(49,1);

%choose a starting frequency based on input key
switch key
    case 'C'
        init_note = 16.35;
    case 'C#'
        init_note = 17.32;
    case 'D'
        init_note = 18.35;
    case 'D#'
        init_note = 19.45;
    case 'E'
        init_note = 20.6;
    case 'F'
        init_note = 21.83;
    case 'F#'
        init_note = 23.12;
    case 'G'
        init_note = 24.5;
    case 'G#'
        init_note = 25.96;
    case 'A'
        init_note = 27.5;
    case 'A#'
        init_note  = 29.14;
    case 'B'
        init_note = 30.87;
end

%call init_note n, just to condense code
n = init_note;

%generate initial frequencies
init_freqs = [n n*(2^(2/12)) n*(2^(4/12)) n*(2^(5/12)) n*(2^(7/12)) n*(2^(9/12)) n*(2^(11/12))];

%build scale vectors by making additional notes as multiples of previous
for i = 1:7
    if i == 1
        scale(1:7) = init_freqs;
    else
        scale((7*i)-6:(7*i)) = 2^(i-1).*scale(1:7);
    end
end

%adjust frequencies for accuracy
scale = adjustFreq(scale);

end

