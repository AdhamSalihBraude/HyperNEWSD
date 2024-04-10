function [Newindividuals,updated_innovation_record,updated_TopologyTracker] = reproduceTopologies(crossover ,mutation, TopologyTracker ,innovation_record, GenCounter,mating_pool,PopSize,MaxTopologeis,ActivFunc)
yl = -1*mutation.weight_range;
yu = mutation.weight_range;
Ind_index = randperm(length(mating_pool));
Newindividuals = [];
parentIndex = 1;
numberoftopologies = sum([TopologyTracker(:).number_individuals]~=0);
while length(Newindividuals) < PopSize
    if numberoftopologies < MaxTopologeis
        
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
        if rand() < crossover.percentage
            %inherit nodes from both parents
            new_individual.nodegenes=[];
            matrix_node_lineup=[[parent1.nodegenes(1,:);1:size(parent1.nodegenes,2);zeros(1,size(parent1.nodegenes,2))],[parent2.nodegenes(1,:);zeros(1,size(parent2.nodegenes,2));1:size(parent2.nodegenes,2)]];
            [discard,sort_node_vec]=sort(matrix_node_lineup(1,:));
            matrix_node_lineup=matrix_node_lineup(:,sort_node_vec);
            node_number=0;
            for index_node_sort=1:size(matrix_node_lineup,2)
                if node_number~=matrix_node_lineup(1,index_node_sort)
                    if matrix_node_lineup(2,index_node_sort)>0
                        new_individual.nodegenes=[new_individual.nodegenes,parent1.nodegenes(:,matrix_node_lineup(2,index_node_sort))];
                    else
                        new_individual.nodegenes=[new_individual.nodegenes,parent2.nodegenes(:,matrix_node_lineup(3,index_node_sort))];
                    end
                    node_number=matrix_node_lineup(1,index_node_sort);
                end
            end
            %Crossover connection genes
            %first do lineup of connection genes
            matrix_lineup=[[parent1.connectiongenes(1,:);1:size(parent1.connectiongenes,2);zeros(1,size(parent1.connectiongenes,2))],[parent2.connectiongenes(1,:);zeros(1,size(parent2.connectiongenes,2));1:size(parent2.connectiongenes,2)]];
            [discard,sort_vec]=sort(matrix_lineup(1,:));
            matrix_lineup=matrix_lineup(:,sort_vec);
            final_matrix_lineup=[];
            innovation_number=0;
            for index_sort=1:size(matrix_lineup,2)
                if innovation_number~=matrix_lineup(1,index_sort)
                    final_matrix_lineup=[final_matrix_lineup,matrix_lineup(:,index_sort)];
                    innovation_number=matrix_lineup(1,index_sort);
                else
                    final_matrix_lineup(2:3,size(final_matrix_lineup,2))=final_matrix_lineup(2:3,size(final_matrix_lineup,2))+matrix_lineup(2:3,index_sort);
                end
            end
            % O.K. Connection Genes are lined up, start with crossover
            new_individual.connectiongenes=[];
            
            for index_lineup=1:size(final_matrix_lineup,2)
                if (final_matrix_lineup(2,index_lineup)>0) & (final_matrix_lineup(3,index_lineup)>0) %check for matching genes, do crossover
                    if rand<0.5 %random crossover for matching genes
                        new_individual.connectiongenes=[new_individual.connectiongenes,parent1.connectiongenes(:,final_matrix_lineup(2,index_lineup))];
                    else
                        new_individual.connectiongenes=[new_individual.connectiongenes,parent2.connectiongenes(:,final_matrix_lineup(3,index_lineup))];
                    end
                    if rand>crossover.probability_multipoint %weight averaging for offspring, otherwise the randomly inherited weights are left undisturbed
                        new_individual.connectiongenes(4,size(new_individual.connectiongenes,2))=(parent1.connectiongenes(4,final_matrix_lineup(2,index_lineup))+parent2.connectiongenes(4,final_matrix_lineup(3,index_lineup)))/2;
                    end
                end
                parent1_flag=sum(final_matrix_lineup(2,index_lineup+1:size(final_matrix_lineup,2)));  % test if there exist further connection genes from index_lineup+1 to end of final_matrix_lineup for parent1 (to detect excess)
                parent2_flag=sum(final_matrix_lineup(3,index_lineup+1:size(final_matrix_lineup,2)));  % test if there exist further connection genes from index_lineup+1 to end of final_matrix_lineup for parent1 (to detect excess)
                % Two cases to check (excess is taken care of in the disjoint gene checks)
                if (final_matrix_lineup(2,index_lineup)>0) && (final_matrix_lineup(3,index_lineup)==0) %Disjoint parent1
                    %                 if parent1.rank_local < parent2.rank_local || (parent1.rank_local == parent2.rank_local && parent1.distance_local > parent2.distance_local)
                    new_individual.connectiongenes=[new_individual.connectiongenes,parent1.connectiongenes(:,final_matrix_lineup(2,index_lineup))];
                    %                 end
                end
                if (final_matrix_lineup(2,index_lineup)==0) && (final_matrix_lineup(3,index_lineup)>0) %Disjoint parent2
                    %                 if parent1.rank_local > parent2.rank_local || (parent1.rank_local == parent2.rank_local && parent1.distance_local < parent2.distance_local)
                    new_individual.connectiongenes=[new_individual.connectiongenes,parent2.connectiongenes(:,final_matrix_lineup(3,index_lineup))];
                    %                 end
                end
            end
            new_individual.SubProbelmsScore =[]; %has no impact on algorithm, only required for assignment to new population
            new_individual.TID= parent1.TID; %will be species hint for speciation
            new_individual.Fit = [];
            new_individual.NC = [];
        else % no crossover, just copy a individual of the species and mutate in subsequent steps
            new_individual=mating_pool(randperm(length(mating_pool),1));
        end
        
        % Hidden nodes culling (remove any hidden nodes where there is no corresponding connection gene in the new individual)
        connected_nodes=[];
        for index_node_culling=1:size(new_individual.nodegenes,2)
            node_connected_flag=sum(new_individual.connectiongenes(2,:)==new_individual.nodegenes(1,index_node_culling))+sum(new_individual.connectiongenes(3,:)==new_individual.nodegenes(1,index_node_culling));
            if (node_connected_flag>0) || (new_individual.nodegenes(2,index_node_culling)~=3)
                connected_nodes=[connected_nodes,new_individual.nodegenes(:,index_node_culling)];
            end
        end
        new_individual.nodegenes=connected_nodes;
        
        % Disabled Genes Mutation
        %run through all connection genes in a new_individual, find disabled connection genes, enable again with crossover.probability_gene_reenabled probability
        for index_connection_gene=1:size(new_individual.connectiongenes,2)
            if (new_individual.connectiongenes(5,index_connection_gene)==0)&& (rand<mutation.probability_gene_reenabled)
                new_individual.connectiongenes(5,index_connection_gene)=1;
            end
        end
        
        % Weight Mutation
        %run through all connection genes in a new_individual, decide on mutating or not
        for index_connection_gene=1:size(new_individual.connectiongenes,2)
            if rand<mutation.probability_mutate_weight %*index_connection_gene/size(new_individual.connectiongenes,2) %linearly biased towards higher probability of mutation at end of connection genes
                new_individual.connectiongenes(4,index_connection_gene)= new_individual.connectiongenes(4,index_connection_gene)+mutation.weight_range*(rand-0.5);
            end
            % weight capping
            new_individual.connectiongenes(4,index_connection_gene)=new_individual.connectiongenes(4,index_connection_gene)*(abs(new_individual.connectiongenes(4,index_connection_gene))<=mutation.weight_cap)+(sign(new_individual.connectiongenes(4,index_connection_gene))*mutation.weight_cap)*(abs(new_individual.connectiongenes(4,index_connection_gene))>mutation.weight_cap);
        end
        
        % IMPORTANT: The checks for duplicate innovations in the following two types of mutation can only check in the current GenCounter
        % Add Connection Mutation
        flag_recurrency_enabled=rand<mutation.probability_recurrency;
        vector_possible_connect_from_nodes=new_individual.nodegenes(1,:); %connections can run from every node
        vector_possible_connect_to_nodes=new_individual.nodegenes(1,find((new_individual.nodegenes(2,:)==2)+(new_individual.nodegenes(2,:)==3))); %connections can only run into hidden and output nodes
        number_possible_connection=length(vector_possible_connect_from_nodes)*length(vector_possible_connect_to_nodes)-size(new_individual.connectiongenes,2);
        
        flag1=(rand<mutation.probability_add_node);
        
        if (rand<mutation.probability_add_connection) && (number_possible_connection>0) && (flag1==0) %check if new connections can be added to genes (if there are any possible connections which are not already existing in genes of new individual)
            % First build matrix containing all possible new connection for nodegene of new individual
            
            new_connection_matrix=[];
            for index_connect_from=1:length(vector_possible_connect_from_nodes)
                for index_connect_to=1:length(vector_possible_connect_to_nodes)
                    possible_connection=[vector_possible_connect_from_nodes(index_connect_from);vector_possible_connect_to_nodes(index_connect_to)];
                    if sum((new_individual.connectiongenes(2,:)==possible_connection(1)).*(new_individual.connectiongenes(3,:)==possible_connection(2)))==0 % Check if proposed connection is not already contained in gene
                        new_connection_matrix=[new_connection_matrix,possible_connection];
                    end
                end
            end
            % Shuffle possible new connections randomly
            [discard,shuffle]=sort(rand(1,size(new_connection_matrix,2)));
            new_connection_matrix=new_connection_matrix(:,shuffle);
            
            index_new_connection=0;
            flag_connection_ok=0;
            % check if connection is o.k. (meaning either non-recurrent or recurrent and flag_recurrency_enabled set to 1) if not connection is found which is o.k.,no connection will be added to connection genes of new individual
            while (flag_connection_ok==0) && (index_new_connection<size(new_connection_matrix,2))
                index_new_connection=index_new_connection+1;
                new_connection=new_connection_matrix(:,index_new_connection);
                
                % test new connection if it is recurrent (i.e. at least one of the possibles path starting from connect_to node in the network leads back to the connect_from node
                flag_recurrent=0;
                if new_connection(1)==new_connection(2) %trivial recurrency
                    flag_recurrent=1;
                end
                nodes_current_level=new_connection(2);
                depth=0;
                while flag_recurrent==0 && depth<size(new_individual.connectiongenes,2) && ~isempty(nodes_current_level)
                    depth=depth+1;
                    nodes_next_level=[];
                    for index_check=1:size(nodes_current_level);
                        nodes_next_level=[nodes_next_level,new_individual.connectiongenes(3,find(new_individual.connectiongenes(2,:)==nodes_current_level(index_check)))];
                    end
                    if sum(nodes_next_level(:)==new_connection(1))>0
                        flag_recurrent=1;
                    end
                    nodes_current_level=nodes_next_level;
                end
                if flag_recurrent==0
                    flag_connection_ok=1;
                elseif flag_recurrency_enabled
                    flag_connection_ok=1;
                end
            end
            
            % Now we test if it is a true innovation (i.e. hasn't already happened in this GenCounter) we can only do this if a valid new connection has been found
            if flag_connection_ok
                index_already_happened=find((innovation_record(5,:)==GenCounter).*(innovation_record(2,:)==new_connection(1)).*(innovation_record(3,:)==new_connection(2))); %set flag signifying new innovation (connection not contained in innovation_record of this GenCounter)
                new_innovation=not(sum(index_already_happened));
                if new_innovation==1 % O.K. is new innovation
                    new_connection=[max(innovation_record(1,:))+1;new_connection]; %Update the new connection with its innovation number
                    % Update connection_genes
                    new_individual.connectiongenes=[new_individual.connectiongenes,[new_connection;rand*2-1;1]];
                    % Update innovation_record
                    innovation_record=[innovation_record,[new_connection;0;GenCounter]];
                else % connection gene already exists in innovation_record of this GenCounter
                    % Update connection_genes
                    new_individual.connectiongenes=[new_individual.connectiongenes,[innovation_record(1:3,index_already_happened);rand*2-1;1]];
                end
            end
        end
        % Add (Insert) Node Mutation
        new_innovation = 0;
        if flag1==1
            max_old_innovation_number=max((innovation_record(5,:)<GenCounter).*innovation_record(1,:)); %highest innovation number from last GenCounter (to ensure that only connections from from last GenCounter or older are chosen for add node mutation, otherwise a new connection added in the last mutation might instantly be disabled)
            vector_possible_connections=[new_individual.connectiongenes(2:3,find((new_individual.connectiongenes(5,:)==1) & (new_individual.connectiongenes(1,:)<=max_old_innovation_number)));find((new_individual.connectiongenes(5,:)==1) & (new_individual.connectiongenes(1,:)<=max_old_innovation_number))];  %compute vector of connections into which a new node could be inserted and their positions in the connection_gene matrix. This vector is composed of all nondisabled connections which stem at least from the last GenCounter or older
            insert_node_connection=vector_possible_connections(:,round(rand*size(vector_possible_connections,2)+0.5));
            new_innovation=1; %set provisionally to 1, will be checked
            exist_innovation=find((innovation_record(5,:)==GenCounter).*(innovation_record(4,:)>0).*(innovation_record(2,:)==insert_node_connection(1))); %Beginning of check innovation record to test for real innovation. exist_innovation contains vector of index of elements in innovation record which fulfil three things: current GenCounter, add node mutation and same connect from as current innovation
            if sum(exist_innovation)>0 %if these are fulfilled, we have to test for connect_to node to see if innovation really is the same
                for index_check=1:length(exist_innovation)
                    if innovation_record(3,exist_innovation(index_check)+1)==insert_node_connection(2)
                        new_innovation=0;
                        index_already_existent_this_GenCounter=exist_innovation(index_check);
                        
                    end
                end
            end
            if new_innovation==1 %O.K. is true innovation for current GenCounter
                % Update node_genes
                new_node_number=max(innovation_record(4,:))+1;
                new_individual.nodegenes=[new_individual.nodegenes,[new_node_number;3;0;0;randi([1,ActivFunc])]];
                % Update connection_genes
                new_individual.connectiongenes(5,insert_node_connection(3))=0; %disable old connection gene
                new_connections=[[max(innovation_record(1,:))+1;insert_node_connection(1);new_node_number;1;1],[max(innovation_record(1,:))+2;new_node_number;insert_node_connection(2);new_individual.connectiongenes(4,insert_node_connection(3));1]];
                new_individual.connectiongenes=[new_individual.connectiongenes,new_connections]; %extend connection_genes by the two new connections
                % Update innovation_record
                innovation_record=[innovation_record,[new_connections(1:3,:);new_node_number,0;GenCounter,GenCounter]];
            else %no new innovation, has already happened at least once in this GenCounter
                % Update node_genes
                node_number=innovation_record(4,index_already_existent_this_GenCounter);
                new_individual.nodegenes=[new_individual.nodegenes,[node_number;3;0;0;randi([1,ActivFunc])]];
                % Update connection_genes
                new_individual.connectiongenes(5,insert_node_connection(3))=0; %disable old connection gene
                new_connections=[innovation_record(1:3,index_already_existent_this_GenCounter:index_already_existent_this_GenCounter+1);1,new_individual.connectiongenes(4,insert_node_connection(3));1,1];
                length_con_gen=size(new_individual.connectiongenes,2); %length of the connection genes of current new_individual
                if new_individual.connectiongenes(1,length_con_gen)>new_connections(1,2) % check if there was an add_connection_mutation to current new_individual which has a higher innovation number than current add_node_mutation
                    new_individual.connectiongenes=[new_individual.connectiongenes(:,1:length_con_gen-1),new_connections,new_individual.connectiongenes(:,length_con_gen)];
                    
                    
                else
                    new_individual.connectiongenes=[new_individual.connectiongenes,new_connections];
                end
            end
            if rand() < mutation.mutateActivFunc
                new_individual.nodegenes(5,:) = randi([1,ActivFunc],size(new_individual.nodegenes(5,:)));
            end
        end
        %% TopologyTracker
        % Loop through comparison vector
        topology_assigned=0;
        index_topology_ref=0;
        while topology_assigned==0 & index_topology_ref<size(TopologyTracker,2)
            %extract reference_individual from reference population
            index_topology_ref=index_topology_ref+1;
            %run through both connection genes, compute disjoint, excess, and average weight difference
            max_num_genes=max([size(new_individual.connectiongenes,2),size(TopologyTracker(index_topology_ref).connections,2)]);
            max_num_innovation=max([new_individual.connectiongenes(1,:),TopologyTracker(index_topology_ref).connections(1,:)]);
            vector_innovation_new=[zeros(1,max(new_individual.connectiongenes(1,:))),ones(1,max_num_innovation-max(new_individual.connectiongenes(1,:)))];
            vector_innovation_new(new_individual.connectiongenes(1,:))=2;
            vector_weight_new=zeros(1,max_num_innovation);
            vector_weight_new(new_individual.connectiongenes(1,:))=new_individual.connectiongenes(4,:);
            vector_innovation_ref=[4*ones(1,max(TopologyTracker(index_topology_ref).connections(1,:))),8*ones(1,max_num_innovation-max(TopologyTracker(index_topology_ref).connections(1,:)))];
            vector_innovation_ref(TopologyTracker(index_topology_ref).connections(1,:))=16;
            vector_weight_ref=zeros(1,max_num_innovation);
            vector_weight_ref(TopologyTracker(index_topology_ref).connections(1,:))=TopologyTracker(index_topology_ref).connections(4,:);
            vector_lineup=vector_innovation_new+vector_innovation_ref;
            excess=sum(vector_lineup==10)+sum(vector_lineup==17);
            disjoint=sum(vector_lineup==6)+sum(vector_lineup==16);
            vector_matching=find(vector_lineup==18);
            average_weight_difference=sum(abs(vector_weight_new(vector_matching)-vector_weight_ref(vector_matching)))/length(vector_matching);
            max_num_genes=1;
            idxNotInputBias = TopologyTracker(index_topology_ref).nodes(1,:) ~= 1 & TopologyTracker(index_topology_ref).nodes(1,:) ~= 4;
            TopologyrefActiv = TopologyTracker(index_topology_ref).nodes(5,idxNotInputBias);
            idxNotInputBias = new_individual.nodegenes(1,:) ~= 1 & new_individual.nodegenes(1,:) ~= 4;
            newIndActiv = new_individual.nodegenes(5,idxNotInputBias);
            notSamActive = 1;
            if length(newIndActiv) == length(TopologyrefActiv)
                notSamActive = all(newIndActiv==TopologyrefActiv);
            end
            distance=excess+disjoint+notSamActive;
            if distance == 0
                % assign individual to same species as current reference individual
                
                new_individual.TID = TopologyTracker(index_topology_ref).ID;
                topology_assigned=1; %set flag indicating new_individual has been assigned to Topology
                if TopologyTracker(index_topology_ref).number_individuals == 0
                    new_individual = [];
                end
            end
        end
        % not compatible with any? well, then create new topology
        if topology_assigned==0
            new_Topology_ID=size(TopologyTracker,2)+1;
            % assign individual to new species
            new_individual.TID=new_Topology_ID;
            % update species_record
            TopologyTracker(new_Topology_ID).ID=new_Topology_ID;
            TopologyTracker(new_Topology_ID).nodes = new_individual.nodegenes;
            TopologyTracker(new_Topology_ID).connections = new_individual.connectiongenes;
            TopologyTracker(new_Topology_ID).gen = GenCounter;
            TopologyTracker(new_Topology_ID).protected = 1;
            numberoftopologies = numberoftopologies + 1;
        end
        Newindividuals = [Newindividuals,new_individual];
        clear new_individual
    else
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
end
if length(Newindividuals) > PopSize
    Newindividuals = Newindividuals(1:PopSize);
end
% final update of species_record (can only be done now since old population sizes were needed during reproduction cycle)
for index_topology=1:size(TopologyTracker,2)
    TopologyTracker(index_topology).number_individuals=sum([Newindividuals(:).TID]==index_topology);
end
%assign updated species_record to output
updated_TopologyTracker = TopologyTracker;
%assign updated innovation_record to output
updated_innovation_record=innovation_record;

