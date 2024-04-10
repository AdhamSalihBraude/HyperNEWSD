function [NEWpopulation, UpdatedTR] = updatePopulation(Newindividuals,E,PS,TopologyTracker,NSubProblems)
%% Saving the best Inds
NEWpopulation = [];

UpdatedTR = TopologyTracker;
EidxCoped = [];
%apply E-set

Etemp = [E,Newindividuals];
for index_individual = 1 : length(Etemp)
    EtempScore(index_individual,:) = Etemp(index_individual).SubProbelmsScore;
end
[V,IndexBestScore] = min(EtempScore,[],1);
IndexBestScore = unique(IndexBestScore);
% Copy E
NEWpopulation = [NEWpopulation,Etemp(IndexBestScore)];

%apply Protection
for index_individual = 1 : length(NEWpopulation)
    NewTopology(index_individual) = NEWpopulation(index_individual).TID;
end
for i = 1 : length(PS)
    protectedTid = PS(i).TID;
    if ~any(NewTopology == protectedTid)
        NEWpopulation = [NEWpopulation,PS(i)];
    end
end
% add inds until PopSize = length(pop)
while length(NEWpopulation) < length(Newindividuals)
    NEWpopulation = [NEWpopulation, Newindividuals(randi(length(Newindividuals)))];
end

% update TR
for index_topology=1:size(UpdatedTR,2)
    UpdatedTR(index_topology).number_individuals=sum([NEWpopulation(:).TID]==index_topology);
end

end