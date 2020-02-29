% [X,T] = simpleseries_dataset;

clear

N = 10000;

seqLen = 2;

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

% X = [X; X1; X2];

% vmap = unique(X', 'rows')';
% vocabSize = size(vmap,2)
% 
% for i=1:vocabSize
%         indx = find(all(bsxfun(@eq, X', vmap(:,i)'), 2));
% %         [~,indx]=ismember(vmap(:,i)',Xlong','rows')
%         Xmap(indx) = i;
% end
% 
% 
% 
% X = Xmap;
T = X (:,seqLen+1:end);


X = X(:,1:length(T));

% na=1;
% nb=2;
% % sys=arimax(X,[na nb]);
% % sys=armax(X',[na nb]);
% % y=predict(sys,X',2);


T = num2cell(T);


X = num2cell(X );


% net = layrecnet(1:seqLen,10);
% [Xs,Xi,Ai,Ts] = preparets(net,X,T);
% % net = train(net,Xs,Ts,Xi,Ai);
% [net, tr] = traingdm(net,Xs,Ts,Xi,Ai);
% view(net)
% Y = net(Xs,Xi,Ai);

% net = narnet(1:seqLen,10);
net = narnet(1:2,10);
[Xs,Xi,Ai,Ts] = preparets(net,{},{},X);
% [Xs,Xi,Ai,Ts] = preparets(net,X,T);
net = train(net,Xs,Ts,Xi,Ai);
view(net)
Y = net(Xs,Xi);


perf = perform(net,Ts,Y)



yvec =  cell2mat(Y);

% yvec = yvec(1:end-2);
%  xvec = cell2mat(X);
% xvec = xvec(5:end);
 xvec = cell2mat(Xs);

yvec = round(yvec);
numErrors = sum(abs(sign(xvec - yvec)))

percErrors = numErrors/length(xvec)
% 
% netc = closeloop(net);
% [xc,xic,aic,tc] = preparets(netc,{},{},T);
% yc = netc(xc,xic,aic);
% perfc = perform(net,tc,yc)
% 
% ycvec =  cell2mat(yc);
% tcvec = cell2mat(tc);
% 
% ycvec = round(ycvec);
% numErrors = sum(abs(sign(tcvec - ycvec)))


Ypredict = net(Xs(1:2),Xi)

nets = removedelay(net);
[xs,xis,ais,ts] = preparets(nets,{},{},T);
ys = nets(xs,xis,ais);
stepAheadPerformance = perform(net,ts,ys)

Ypredict2 = nets(Xs,[1 2],ais)



