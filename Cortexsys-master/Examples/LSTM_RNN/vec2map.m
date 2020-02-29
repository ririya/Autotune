function V = vec2map(X, vmap)  
    vocabSize = size(vmap,2);

%     vec = double(str);
    % Map all ASCII values to the reduced set
%     for i=1:vocabSize
%         vec(vec==vmap(i)) = i;
%     end

Xmap = zeros(1,size(X,2));


for i=1:size(X,2)
    indx = find(ismember(vmap',X(:,i)','rows')); 
     Xmap(i) = indx;    
end

%    for i=1:vocabSize
%         indx = find(all(bsxfun(@eq, X', vmap(:,i)'), 2));
% %         [~,indx]=ismember(vmap(:,i)',Xlong','rows')
%         Xmap(indx) = i;
%     end
    
    V = speye(vocabSize);
    V = V(Xmap,:)';
end