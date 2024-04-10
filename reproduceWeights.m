function Newindividuals  = reproduceWeights(crossover ,mutation, GenCounter,mating_pool,PopSize)
%% ver1 adaptive mutation 

Ind_index = randperm(length(mating_pool));
Newindividuals = [];
parentIndex = 1;
while length(Newindividuals) < PopSize
    if parentIndex + 1 < length(mating_pool)
        parent1 = mating_pool(Ind_index(parentIndex));
        parent2 = mating_pool(Ind_index(parentIndex + 1));
        parentIndex = parentIndex + 2;
    else
        Index1 = randi(length(mating_pool));
        Index2 = randi(length(mating_pool));
        while Index1 == Index2
            Index2 = randi(length(mating_pool));
        end
        parent1 = mating_pool(Index1);
        parent2 = mating_pool(Index2);
    end
    if parent1.TID == parent2.TID
        [child1,child2] = garep(parent1,parent2,mutation,crossover);
    else
        child1 = esrep(parent1,mutation,GenCounter);
        child2 = esrep(parent2,mutation,GenCounter);
    end
    if length(Newindividuals) + 1 <  PopSize
            Newindividuals = [Newindividuals,child1,child2];
        else
            if rand() < 0.5
                Newindividuals = [Newindividuals,child1];
            else
                Newindividuals = [Newindividuals,child2];
            end
    end
end
if length(Newindividuals) > PopSize
    Newindividuals = Newindividuals(1:PopSize);
end