function [Newindividuals] = ScoreAssigment(SubProblems,Newindividuals,idealpoint,Gen_ratio)
NSubProblems = length(SubProblems);
%% del unfeasble solutions
for indIndex = 1 : length(Newindividuals)
    if sum(isnan(Newindividuals(indIndex).Fit)) > 0
        Newindividuals(indIndex).Fit = inf(size(Newindividuals(indIndex).Fit));
    end
end
%%
if NSubProblems == 1
    for index_individual = 1 : length(Newindividuals)
        Newindividuals(index_individual).SubProbelmsScore = teNormNav(SubProblems.weight, Newindividuals(index_individual).Fit);
    end
else
    %% Calculate the Score for each subproblem
    for i = 1 : NSubProblems
        for index_individual = 1 : length(Newindividuals)
             Newindividuals(index_individual).SubProbelmsScore(i) = te(SubProblems(i).weight, [Newindividuals(index_individual).Fit]',idealpoint);
%             if Gen_ratio > 1
%                 Newindividuals(index_individual).SubProbelmsScore(i) = teNormNav(SubProblems(i).weight, Newindividuals(index_individual).Trainfitness./idealpoint);
%             else
              %  Newindividuals(index_individual).SubProbelmsScore(i) = teNorm(SubProblems(i).weight, [Newindividuals(index_individual).Fit]'./idealpoint);
%             end
        end
    end
end
end



function gte = te(Lmbda, F, Z)
gte = max((Lmbda+0.00001*(Lmbda==0)).*abs(F - Z));
end

function gte = teNorm(Lmbda, F)
gte = max((Lmbda+0.00001*(Lmbda==0)).*abs(F - ones(size(F))));
end
function gte = teNormNav(Lmbda, F)
panelty = 1/((10*sum((Lmbda>0).*(F==0))+1));
Fc = panelty.*F;
gte = max((Lmbda+0.00001*(Lmbda==0)).*abs(Fc - ones(size(Fc))));
end