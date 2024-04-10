function EP = updateEP(EP,Offsprings,nEP)
% Update the external population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Lucas Farias

    %% Select the non-dominated solutions
    EP = [EP,Offsprings];
    F = [EP(:).Fit]';
    [~,ia,~] = unique(F,'rows');
    EP = EP(ia);
    F = [EP(:).Fit]';
    EP = EP(NDSort(F,1)==1);
    N = length(EP);
    M = length(EP(1).Fit);
    F = [EP(:).Fit]';
	
    %% Delete the overcrowded solutions
    Dis = pdist2(F,F);
    Dis(logical(eye(length(Dis)))) = inf;
    Del = false(1,N);
    while sum(Del) < N-nEP
        Remain = find(~Del);
        subDis = sort(Dis(Remain,Remain),2);
        [~,worst] = min(prod(subDis(:,1:min(M,length(Remain))),2));
        Del(Remain(worst)) = true;
    end   
    
    EP = EP(~Del);
	
end