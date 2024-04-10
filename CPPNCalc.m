%% NetCalc

function output = CPPNCalc(Ind,Input)
ActivationFunc(1).F = @(x) sin(x);
ActivationFunc(2).F = @(x) tanh(x);
ActivationFunc(3).F = @(x) gaussmf(x,[1,0]);
ActivationFunc(4).F = @(x) cos(x);
ActivationFunc(5).F = @(x) 1/(1+exp(-x));
ActivationFunc(6).F = @(x) (x + abs(x))/2;
ActivationFunc(7).F = @(x) x;
nodegenes = Ind.nodegenes;
connectiongenes = Ind.connectiongenes;
OutputSize = sum(nodegenes(2,:)==2);
number_input_nodes = length(Input);
no_change_threshold=1e-3; %threshold to judge if state of a node has changed significantly since last iteration
number_nodes=size(nodegenes,2);
number_connections=size(connectiongenes,2);
% set node input steps for first timestep
nodegenes(3,number_input_nodes+1:number_nodes) = 0; %set all node input states to zero
nodegenes(3,number_input_nodes+1) = 1; %bias node input state set to 1
nodegenes(3,1:number_input_nodes) = Input; %node input states of the two input nodes are consecutively set to the XOR input pattern
%set node output states for first timestep (depending on input states)
nodegenes(4,1:number_input_nodes) = nodegenes(3,1:number_input_nodes);
for idx = number_input_nodes+2:number_nodes
    F = ActivationFunc(nodegenes(5,idx)).F;
    nodegenes(4,idx) = F(nodegenes(3,idx));
end
% nodegenes(4,number_input_nodes+1:number_nodes) = F(nodegenes(3,number_input_nodes+1:number_nodes));
no_change_count = 0;
index_loop = 0;
while (no_change_count<number_nodes) && index_loop<3*number_connections
    index_loop = index_loop + 1;
    vector_node_state = nodegenes(4,:);
    for index_connections = 1 : number_connections
        if connectiongenes(5,index_connections)==1 %Check if Connection is enabled
            %read relevant contents of connection gene (ID of Node where connection starts, ID of Node where it ends, and connection weight)
            ID_connection_from_node = connectiongenes(2,index_connections);
            ID_connection_to_node = connectiongenes(3,index_connections);
            connection_weight = connectiongenes(4,index_connections);
            %map node ID's (as extracted from single connection genes above) to index of corresponding node in node genes matrix
            index_connection_from_node = nodegenes(1,:)==ID_connection_from_node;
            index_connection_to_node = find(nodegenes(1,:)==ID_connection_to_node);
            nodegenes(3,index_connection_to_node) = nodegenes(3,index_connection_to_node)+nodegenes(4,index_connection_from_node)*connection_weight; %take output state of connection_from node, multiply with weight, add to input state of connection_to node
        end
    end
    %pass on node input states to outputs for next timestep
    for idx = number_input_nodes+2:number_nodes
        F = ActivationFunc(nodegenes(5,idx)).F;
        nodegenes(4,idx) = F(nodegenes(3,idx));
    end
    %Re-initialize node input states for next timestep
    nodegenes(3,number_input_nodes+1:number_nodes) = 0; %set all output and hidden node input states to zero
    no_change_count = sum(abs(nodegenes(4,:)-vector_node_state)<no_change_threshold); %check for alle nodes where the node output state has changed by less than no_change_threshold since last iteration through all the connection genes
end
if index_loop >= 3*number_connections
    output = nan(size(nodegenes(4,number_input_nodes+1 : number_input_nodes+1 + OutputSize-1)));
else
    outputidx = nodegenes(2,:)==2;
    output = nodegenes(4,outputidx);
end
end