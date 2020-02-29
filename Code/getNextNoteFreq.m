function [Y,maxProb,bestSeqLen] = getNextNoteFreq(Xreal, dictionary,vmap,maxSeqLen)

if maxSeqLen > length(dictionary)
    maxSeqLen = length(dictionary);
end

% for i=1:size(vmap,2)
%     indx = find(all(bsxfun(@eq, Xreal', vmap(:,i)'), 2));
%     X(indx) = i;
% end

for i=1:size(Xreal,2)
    indx = find(all(bsxfun(@eq, vmap', Xreal(:,i)'), 2));
    X(i) = indx;
end

maxProb = 0;

lastIndex = length(X);

for seqLen = 1:maxSeqLen
    
    if  lastIndex < seqLen
        continue;
    end
    
    currSeq = X(lastIndex-seqLen+1:lastIndex);
    
    error = repmat(currSeq,size(dictionary{seqLen}.SequenceList,1),1) - dictionary{seqLen}.SequenceList;
    indSeq = find(sum(abs(error),2) == 0);    %consider changing to smallest error
    
    if isempty(indSeq)
        break;
    end
    
    %     for k = 1:length(dictionary{seqLen2}.Sequences)
    %
    %         if (dictionary{seqLen2}.Sequences{k}.Seq == currSeq) %consider changing to smallest error
    %             indSeq = k;
    %             break;
    %         end
    %
    %     end
    
    
    [currMax indMax] = max(dictionary{seqLen}.Sequences{indSeq}.prob);
    
    if (currMax >= maxProb)   %if prob is equal, favor longer sequences
        maxProb = currMax;
        Y = vmap(:,indMax);
        bestSeqLen = seqLen;
    end
    
    
end

end
