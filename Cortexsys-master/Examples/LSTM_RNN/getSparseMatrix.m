function [X, vmap] = getSparseMatrix(Xlong, seqLen, offset,vmap)

Xlong = Xlong(:,1+offset:end);
N = size(Xlong,2);
    
    m = floor(N/seqLen);
    
    Xlong = Xlong(1:m*seqLen);
    
    dimX = size(Xlong,1);    

 Npad = size(Xlong,2);
 
 if isempty(vmap)

vmap = unique(Xlong', 'rows')';
vocabSize = size(vmap,2);
indRand= randperm(vocabSize);
vmap = vmap(indRand);
 end

vocabSize = size(vmap,2);



  Xmap = zeros(1,Npad);

  % Map all ASCII values to the reduced set
    for i=1:vocabSize
        indx = find(all(bsxfun(@eq, Xlong', vmap(:,i)'), 2));
%         [~,indx]=ismember(vmap(:,i)',Xlong','rows')
        Xmap(indx) = i;
    end
    
     Xmap = reshape(Xmap,seqLen,[]);
%  
%      Xmap= hankel(Xmap, 1:seqLen);
%      Xmap = Xmap';%      

     xreal = vmap(Xmap);

  Xca = cell(seqLen,1);
    I = speye(vocabSize); 

   I = speye(vocabSize); 

    for t=1:seqLen
        Xt = Xmap(t,:);
        Xca{t} = I(Xt(:),:)';
        
    end
   
%     Xfull1 = full(Xca{1});
%     Xfull2 = full(Xca{2});
    
     X = Xca;
    