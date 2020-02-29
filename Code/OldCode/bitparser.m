function [ pairs ] = bitparser( bitvector )
%bitparser Summary of this function goes here
%   Given a vector of 1s and 0s return a vector of starting and ending
%   indicies
    row = bitvector(:); % make the vector into a row vector
    
    change01 = find(((row(1:end-1) == 0).*(row(2:end) == 1)) == 1)+1;
    change10 = find(((row(1:end-1) == 1).* (row(2:end) == 0)) == 1);
    
    if bitvector(1)==1
       change01 = [1; change01];
    end
    if bitvector(end) == 1
       change10 = [change10; length(bitvector)]; 
    end
    
    pairs = [change01,change10];
    

end

