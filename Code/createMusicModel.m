function [dictionary, vmap] =  createMusicModel(X,maxSeqLen)

isCell = iscell(X);

if isCell
    
    Xcell = X;
    X = [];
    for c=1:length(Xcell)
        X = [X Xcell{c}'];
    end
    
else
    Xcell{1} = X;
end

vmap = unique(X', 'rows')';
vocabSize = size(vmap,2);

tstart = tic;

dictionary = {};

for c=1:length(Xcell)
    
    c
    
    X = Xcell{c}';
    
    if size(X,1) > size(X,2)
        X = X';
    end
    
    %
    for i=1:vocabSize
        indx = find(all(bsxfun(@eq, X', vmap(:,i)'), 2));
        Xmap(indx) = i;
        
    end
    
%     xreal = vmap(:,Xmap);
    
    X = Xmap;
    

    
    
    for seqLen2 = 1:maxSeqLen
        
        c
        seqLen2
        
        for i=1:length(X)-seqLen2
            currSeq = X(i-1+1:i-1+seqLen2);
            
            indSeq = 0;
            foundInd = 0;
            
            if seqLen2>(length(dictionary))
                dictionary{seqLen2}.Sequences = {};
            end
            
            if isempty(dictionary{seqLen2}.Sequences)   % if dictionary is empty for this sequence length, add first entry (no entry found)
                indSeq =1;              
                dictionary{seqLen2}.Sequences{indSeq}.Seq = currSeq;
                dictionary{seqLen2}.SequenceList = currSeq;
                dictionary{seqLen2}.Sequences{indSeq}.count = zeros(1,vocabSize);
            else  % look for entry in the dictionary, 
                
                error = repmat(currSeq,size(dictionary{seqLen2}.SequenceList,1),1) - dictionary{seqLen2}.SequenceList;
                indMatch = find(sum(abs(error),2) == 0);
                
                if ~isempty(indMatch)                  % if there is a match, probabilities should not be updated
                    indSeq =   indMatch;
                    foundInd = 1;
                end
                
%                 for k = 1:length(dictionary{seqLen2}.Sequences)
%                     if (dictionary{seqLen2}.Sequences{k}.Seq == currSeq)
%                         indSeq = k;
%                         foundInd = 1;
%                         break
%                     end
%                 end
                if indSeq==0        %
      
                    indSeq = length( dictionary{seqLen2}.Sequences) +1;
                    dictionary{seqLen2}.SequenceList = [dictionary{seqLen2}.SequenceList; currSeq];
                    dictionary{seqLen2}.Sequences{indSeq}.Seq = currSeq;
                    dictionary{seqLen2}.Sequences{indSeq}.count = zeros(1,vocabSize);
                    
                end
            end
            
            if foundInd == 0     % if first time this sequence was found, update probabilities
                
                for j=1:vocabSize
                    %        [currSeq vmap(j)]
                    %     indNumber = findstr(X, [currSeq  vmap(j)]);
                    indNumber = strfind(X, [currSeq  j]);
                    
                    dictionary{seqLen2}.Sequences{indSeq}.count(j) = dictionary{seqLen2}.Sequences{indSeq}.count(j) + length(indNumber);
                    
                end
                count = dictionary{seqLen2}.Sequences{indSeq}.count;
                dictionary{seqLen2}.Sequences{indSeq}.prob = count./sum(count);
                
            end
            
        end
        
    end
    
end

toc(tstart)






