clear

N = 10000;

seqLen = 2;

X = 1:seqLen;
for i=1:N
     X(end+1) = 0;
for j = 1:seqLen
   X(end) = X(end) + X(end-j);
end
  if abs(X(end)) > 2
       X(end) = -X(end); 
   end


end

X1 = circshift(X,-1);

X = [X; X1];
X = X(:,1:end-1);  %circshift makes the last sample be different

maxSeqLen = 5;
[dictionary, vmap] =  createMusicModel(X,maxSeqLen);

Y = X(:,1:2);   
       
for i=length(Y):length(X) - 1
    
    Ynext = getNextNoteFreqAndTime(Y, dictionary,vmap);
    
    Y = [Y Ynext];
    
end


  for i=1:size(vmap,2)
        indx = find(all(bsxfun(@eq, X', vmap(:,i)'), 2));
        indy = find(all(bsxfun(@eq, Y', vmap(:,i)'), 2));
        Xmap(indx) = i;
        Ymap(indy) = i;

  end
 
difference = Xmap - Ymap;
errors = sum(difference)  

  



