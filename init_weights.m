function subp=init_weights(NSubProblems, objDim)
% init_weights function initialize a pupulation of subproblems structure
% with the generated decomposition weight and the neighbourhood
% relationship.W = [];
if NSubProblems ~= 1
    iter = 0;
        [W,~] = UniformPoint(NSubProblems,objDim);
        W = W((sum((W<1e-5),2)==0|sum((W<1e-5),2)==objDim-1),:);
        while size(W,1) < NSubProblems
            [W1,~] = UniformPoint(NSubProblems + iter ,objDim);
%             W1 = W1((sum((W1<1e-5),2)==0),:);
            W = [W;W1];
            W = unique(W,'rows');
            iter = iter + 1;
        end
        I = randperm(size(W,1),NSubProblems);
        W = W(I,:);
    %% Generalized Decomposition
    W(W<1e-5) = 0.001;
    W = 1./W;
    W = W./sum(W,2);
    %%end of Generalized Decomposition 
else
    W = ones(objDim,1)*1/objDim;
end
subp=[];
for i=1:NSubProblems
    p=struct('weight',[]);
    p.weight=W(i,:);
    subp=[subp p];
    
end


%Set up the neighbourhood. (NR)
% leng=length(subp);
% distanceMatrix=zeros(leng, leng);
% for i=1:leng
%     for j=i+1:leng
%         A=subp(i).weight;B=subp(j).weight;
%         distanceMatrix(i,j)=(A-B)*(A-B)';
%         distanceMatrix(j,i)=distanceMatrix(i,j);
%     end
%     [s,sindex]=sort(distanceMatrix(i,:));
%     subp(i).neighbour=sindex(1:niche)';
% end
% 
% end