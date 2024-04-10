function PS = protectedSet(population,TopologyTracker,idealpoint)
PS = [];
for Topologyindex = 1 : length(TopologyTracker)
    if TopologyTracker(Topologyindex).number_individuals ~= 0 && TopologyTracker(Topologyindex).protected == 1
        TopologyInds = population([population.TID] == TopologyTracker(Topologyindex).ID);
        TopologyFitI = mean(abs([TopologyInds(1).Fit]'-idealpoint));
        I = 1;
        for i = 2 : length(TopologyInds)
            TopologyFit = mean(abs([TopologyInds(i).Fit]'-idealpoint));
            if TopologyFitI > TopologyFit
                I = i;
                TopologyFitI = TopologyFit;
            end
        end
        PS = [PS,TopologyInds(I)];
    end
end