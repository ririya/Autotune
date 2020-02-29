function [Y,maxProb,bestSeqLen] = getNextNoteFreqAndTime(Xreal, dictionary,vmap,maxSeqLen,maxTimeScale,originalDuration)


% Y = NaN;

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
    

%     [currMax indMax] = max(dictionary{seqLen2}.Sequences{indSeq}.prob);

%   if (currMax > maxProb)
%         maxProb = currMax;
%         Y = vmap(:,indMax);
%     end

    [sortedProbs sortedInd] = sort(dictionary{seqLen}.Sequences{indSeq}.prob,2,'descend');

    indNon0 = find(sortedProbs>0);
    sortedProbs = sortedProbs(indNon0);
    sortedInd = sortedInd(indNon0);
    
    for k=1:length(sortedProbs)
        
        if sortedProbs(k) >= maxProb  %if probabilities are the same, favor longer sequences
            
              %
%             if ( sortedProbs(k) == maxProb) && (unidrnd(N)==1)  %if probabilities are the same, 1/N% chance of keeping the same note 
%                 continue
%             else
                
                candidateProb = sortedProbs(k);
                candidateInd = sortedInd(k);
                candidateY = vmap(:,candidateInd);
                
                timeScale = candidateY(1)/originalDuration;
                
                if (timeScale >= 1/maxTimeScale) &&  (timeScale <= maxTimeScale)
                    Y = candidateY;
                    maxProb = candidateProb;
                    bestSeqLen = seqLen;
                    break;                  %found Y for this seqLen proceed to next seqLen
                end
            
        else
            break; % maxProbability is not higher than what was already selected so no need to keep going
        end
        
    end
    
end

% if isnan(Y)
%     error=1;
% end

end