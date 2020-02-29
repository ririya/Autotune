function all_scales = generateAllScales()

%generate all possible scales using the generateScale function, then load
%into one cell, all_scales.
all_scales = {};

C_scale   = generateScale('C');  all_scales{1}  = C_scale;    
Csh_scale = generateScale('C#'); all_scales{2}  = Csh_scale;  
D_scale   = generateScale('D');  all_scales{3}  = D_scale;    
Dsh_scale = generateScale('D#'); all_scales{4}  = Dsh_scale;  
E_scale   = generateScale('E');  all_scales{5}  = E_scale;    
F_scale   = generateScale('F');  all_scales{6}  = F_scale;    
Fsh_scale = generateScale('F#'); all_scales{7}  = Fsh_scale;  
G_scale   = generateScale('G');  all_scales{8}  = G_scale;    
Gsh_scale = generateScale('G#'); all_scales{9}  = Gsh_scale;  
A_scale   = generateScale('A');  all_scales{10} = A_scale;    
Ash_scale = generateScale('A#'); all_scales{11} = Ash_scale;  
B_scale   = generateScale('B');  all_scales{12} = B_scale;    

all_scales{2,1} = 'C';
all_scales{2,2} = 'C#';
all_scales{2,3} = 'D';
all_scales{2,4} = 'D#';
all_scales{2,5} = 'E';
all_scales{2,6} = 'F';
all_scales{2,7} = 'F#';
all_scales{2,8} = 'G';
all_scales{2,9} = 'G#';
all_scales{2,10} = 'A';
all_scales{2,11} = 'A#';
all_scales{2,12} = 'B';

end