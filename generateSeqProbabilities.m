clear

N = 10000;

seqLen = 3;

X = 1:seqLen;
for i=1:N
     X(end+1) = 0;
for j = 1:seqLen
%    X(end+1) = X(end) + X(end-1) + X(end-2);
   X(end) = X(end) + X(end-j);
end
  if abs(X(end)) > 2
       X(end) = -X(end); 
   end

%   X(end+1) = X(end) + X(end-1);
%    if abs(X(end)) > 2
%        X(end) = -X(end); 
%    end

%  X(end+1) = X(end) +1;
%  if (X(end) > 100)
%      X(end) =1;
%  end
end

% X = repmat(X,3,1);
X1 = circshift(X,-1);
X2 = circshift(X,-2);

X = [X; X1; X2];

vmap = unique(X', 'rows')';
vocabSize = size(vmap,2)
% 
for i=1:vocabSize
        indx = find(all(bsxfun(@eq, X', vmap(:,i)'), 2));
%         [~,indx]=ismember(vmap(:,i)',Xlong','rows')
        Xmap(indx) = i;

end

xreal = vmap(Xmap);

X = Xmap;
T = X (:,seqLen+1:end);


X = X(:,1:length(T));

% seqLen2 = 2;

maxSeqLen = 5;

dictionary = {};

tstart = tic;
for seqLen2 = 1:maxSeqLen
    
for i=1:length(X)-seqLen2
   currSeq = X(i-1+1:i-1+seqLen2);  
   
   indSeq = 0;
   foundInd = 0;
   
   if seqLen2>(length(dictionary))
       dictionary{seqLen2}.Sequences = {}; 
   end
   
   if isempty(dictionary{seqLen2}.Sequences)
       indSeq =1;
       dictionary{seqLen2}.Sequences{indSeq}.Seq = currSeq;
   else
     
       for k = 1:length(dictionary{seqLen2}.Sequences)
           if (dictionary{seqLen2}.Sequences{k}.Seq == currSeq)
               indSeq = k;
               foundInd = 1;
               break
           end                 
       end  
       if indSeq==0
           indSeq = length( dictionary{seqLen2}.Sequences) +1;
           dictionary{seqLen2}.Sequences{indSeq}.Seq = currSeq;   
           
       end
   end
   
   if foundInd == 0
   
   for j=1:vocabSize
%        [currSeq vmap(j)]
%     indNumber = findstr(X, [currSeq  vmap(j)]);
    indNumber = findstr(X, [currSeq  j]);
    
   dictionary{seqLen2}.Sequences{indSeq}.count(j) = length(indNumber);
    
   end
   count = dictionary{seqLen2}.Sequences{indSeq}.count;
   dictionary{seqLen2}.Sequences{indSeq}.prob = count./sum(count);
     
   end
    
end

end

toc(tstart)
a= 1;

Y = X(1:2);

       
 for i=length(Y):length(X) - 1
     
     maxProb = 0;
     
     for seqLen2 = 1:maxSeqLen
                  
         if length(Y) < seqLen2
            continue; 
         end
         
          currSeq = Y(i-seqLen2+1:i);  
         
           for k = 1:length(dictionary{seqLen2}.Sequences)
               
               if (dictionary{seqLen2}.Sequences{k}.Seq == currSeq) %consider changing to smallest error
                   indSeq = k;
                   break;
               end 
               
           end
           
                [currMax indMax] = max(dictionary{seqLen2}.Sequences{indSeq}.prob);
                
                if (currMax > maxProb)
                  maxProb = currMax;  
%                   next = vmap(indMax);
                  next = indMax;
                end
           
     end
           
        Y = [Y next];
         
     end
     
  



