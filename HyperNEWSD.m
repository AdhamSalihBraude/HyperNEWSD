%% HyperNEWSD Scource Code
%% Code for :
%% Adham Salih, and Amiram Moshaiov. "Neuro-evolution-based generic missile guidance law for many-scenarios." Applied Soft Computing 152 (2024): 111210.
%------------------------------- References --------------------------------
% Kenneth O. Stanley, and Risto Miikkulainen. "Evolving neural networks
% through augmenting topologies." Evolutionary Computation 10, no. 2 (2002):
% 99-127. % Coding by Christian Mayr (matlab_neat@web.de)
%---------------------------------------------------------------------------
% Qingfu Zhang, and Hui Li. "MOEA/D: A multiobjective evolutionary
% algorithm based on decomposition." IEEE Transactions on evolutionary
% computation 11.6 (2007): 712-731.
%--------------------------------------------------------------------------
%"Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

clc;clear;clf;close all;
warning('off','all')
RunIdx = 1;
TotalRuns = 31;
% CPPN parameters (TWEANN)

InputLength = 3; % 4 - [from, to, layer, weight or alpha]; 3- [from, to, layer]
outSize = 1; % the weight value

% NC parameters
LayerRes = [5,5,5,1]; % input, Sheets, output
Zres = length(LayerRes) - 1;


ActivFunc = 5; %sin(x);tanh(x);gaussmf(x,[1,0]);cos(x);1/(1+exp(-x));(x + abs(x))/2;x;
% addpath('CommonFunc')
%% Senarios
SenNum = [1:6];
%X = [Xm;Ym;Zm;phiM;gammaM;Vm,Xt;Yt;Zt;phiT;gammaT;integral(am^2);integral(Wt)]
d = 0.005;
pltflg = 0;
tspan = 0:0.01:50;


Senario(1).X0 = [0 0 3000 45*pi/180 -45*pi/180 400 9000 3000 3000 180*pi/180 0 200 0 0];
Senario(1).TarA = 1;
Senario(1).d = d;
Senario(1).tspan = tspan;

Senario(2).X0 = [0 0 3000 45*pi/180 -45*pi/180 400 9000 3000 3000 180*pi/180 0 200 0 0];
Senario(2).TarA = 2;
Senario(2).d = d;
Senario(2).tspan = tspan;

Senario(3).X0 = [0 0 3000 45*pi/180 -45*pi/180 400 9000 3000 3000 180*pi/180 0 200 0 0];
Senario(3).TarA = 3;
Senario(3).d = d;
Senario(3).tspan = tspan;

Senario(4).X0 = [0 0 3000 -75*pi/180 -70*pi/180 400 3000 0 3000 0*pi/180 0 200 0 0];
Senario(4).TarA = 1;
Senario(4).d = 0.05;
Senario(4).tspan = tspan;

Senario(5).X0 = [0 0 3000 -75*pi/180 -70*pi/180 400 3000 0 3000 0*pi/180 0 200 0 0];
Senario(5).TarA = 2;
Senario(5).d = d;
Senario(5).tspan = tspan;

Senario(6).X0 = [0 0 3000 -75*pi/180 -70*pi/180 400 3000 0 3000 0*pi/180 0 200 0 0];
Senario(6).TarA = 3;
Senario(6).d = d;
Senario(6).tspan = tspan;
Senario = Senario(SenNum);
[Senario, J] = calcSatFit(Senario,pltflg); % calculating PPN performances

ObjNum = length(Senario)*2; % the number of objectives to be optimized
%% Run parameters
PopSize = 25*ObjNum/2; % Population Size
MaxGen = 20*(ObjNum/2); % Maximum number of generations
NSubProblems = floor(PopSize*3/4); % number of sub-problems
nEP = floor(PopSize); % number of members in the non-dominated archive
MaxTopologies = floor(PopSize*1/4); % maximum number of evolved topologies
protectGen = floor(MaxGen*1/5); % Number of generations to proteced a novice topology


GenCounter = 1;
gen0 = GenCounter;
%% Run Loop
while RunIdx < TotalRuns + 1
    % if it is a new run
    if gen0 == 1
        %crossover
        crossover.percentage=0.5; %percentage governs the way in which new population will be composed from old population.  exception: species with just one individual can only use mutation
        crossover.probability_interspecies=0.2 ; %if crossover has been selected, this probability governs the intra/interspecies parent composition being used for the
        crossover.probability_multipoint=0.6; %standard-crossover in which matching connection genes are inherited randomly from both parents. In the (1-crossover.probability_multipoint) cases, weights of the new connection genes are the mean of the corresponding parent genes
        crossover.SBXprop = [0.8,0.1]; %first is the probability second is the distibution

        %mutation
        mutation.probability_add_node=0.3;
        mutation.probability_add_connection=0.3;
        mutation.probability_recurrency=0.0; %if we are in add_connection_mutation, this governs if a recurrent connection is allowed. Note: this will only activate if the random connection is a recurrent one, otherwise the connection is simply accepted. If no possible non-recurrent connections exist for the current node genes, then for e.g. a probability of 0.1, 9 times out of 10 no connection is added.
        mutation.probability_mutate_weight=0.6;
        mutation.weight_cap=5; % weights will be restricted from -mutation.weight_cap to mutation.weight_cap
        mutation.weight_range=2; % random distribution with width mutation.weight_range, centered on 0. mutation range of 5 will give random distribution from -2.5 to 2.5
        mutation.probability_gene_reenabled=0.25; % Probability of a connection gene being reenabled in offspring if it was inherited disabled
        mutation.mutateActivFunc = 0.2; % Probability to change the activation function of a node
        mutation.pol = [0.2,15];
        mutation.adaptiverate = [mutation.weight_range,0.2]; %step = a1*exp(-a2*gen)

        ComputerID = getenv('computername');
        gen = 1;
        rng('default')
        warning('off','all')


        %% parameters for the CPPN
        number_input_nodes = InputLength;
        number_output_nodes = outSize;
        vector_connected_input_nodes=1:number_input_nodes; % Try other initial topologies (not fully connected)

        %% First Generation
        idealpoint = inf*ones(1,ObjNum);
        SubProblems = init_weights(NSubProblems, ObjNum);
        NSubProblems = size(SubProblems,2);
        rng('shuffle')
        [population,innovation_record] = init_individual(PopSize,vector_connected_input_nodes,number_input_nodes,number_output_nodes,ObjNum,ActivFunc);
        population = CreatNCs(population, LayerRes, Zres);
        %topology traking
        TopologyTracker = [];
        TopologyTracker(1).ID = 1;
        TopologyTracker(1).nodes = population(1).nodegenes; %node ID; type
        TopologyTracker(1).connections = population(1).connectiongenes;% connction from , connction to
        TopologyTracker(1).gen = 0;
        TopologyTracker(1).number_individuals = length(population);
        TopologyTracker(1).protected = 1;
        population = evalua(population,Senario,pltflg);
        idealpoint = update_idealpoint(idealpoint,population);
        EP = [];
        EP = updateEP(EP,population,nEP);
        population = ScoreAssigment(SubProblems,population,idealpoint,GenCounter/MaxGen);
        SolutionsToKeep = keepsolutions(population,Senario);
        PS = [];
        Data(1).EP = EP;
        Data(1).population = population;
        Data(1).IT = 1;
        gen0 = 2;
    end
    %% Gen Loop

    for GenCounter = gen0 : MaxGen
        Ntopologies = 1;
        IT = [];
        for Topologyindex = 1 : length(TopologyTracker)
            if TopologyTracker(Topologyindex).number_individuals ~= 0
                Ntopologies = Ntopologies + 1;
                IT = [IT,TopologyTracker(Topologyindex).ID];
            end
            if (GenCounter - TopologyTracker(Topologyindex).gen < protectGen && TopologyTracker(Topologyindex).gen ~= 0)
                TopologyTracker(Topologyindex).protected = 1;
            else
                TopologyTracker(Topologyindex).protected = 0;
            end
        end
        E = eliteSet(population,NSubProblems);
        SI = selection(population,PopSize-length([E,PS]));
        matingPool = [SI,PS,E];
        if Ntopologies < MaxTopologies
            [Newindividuals,innovation_record,TopologyTracker] = reproduceTopologies(crossover ,mutation, TopologyTracker ,innovation_record, GenCounter,matingPool,PopSize,MaxTopologies,ActivFunc);
        else
            Newindividuals  = reproduceWeights(crossover, mutation, GenCounter, matingPool, PopSize);
        end
        Newindividuals = CreatNCs(Newindividuals, LayerRes, Zres);
        Newindividuals = evalua(Newindividuals,Senario,pltflg);
        idealpoint = update_idealpoint(idealpoint,Newindividuals);
        SolutionsToKeep = keepsolutions([SolutionsToKeep,Newindividuals],Senario);
        Newindividuals = ScoreAssigment(SubProblems,Newindividuals,idealpoint,GenCounter/MaxGen);
        PS = protectedSet([Newindividuals,population],TopologyTracker,idealpoint);
        if d > 0
            E = evalua(E,Senario,pltflg);
        end
        PS = ScoreAssigment(SubProblems,PS,idealpoint,GenCounter/MaxGen);
        E = ScoreAssigment(SubProblems,E,idealpoint,GenCounter/MaxGen);
        [population,TopologyTracker] = updatePopulation(Newindividuals,E,PS,TopologyTracker,NSubProblems);
        EP = updateEP(EP,Newindividuals,nEP);

        if rem(GenCounter,10)==0
            if d > 0
                EP = evalua(EP,Senario,pltflg);
            end
            if rem(GenCounter,20)==0
                disp([GenCounter, idealpoint,length(SolutionsToKeep)])
            end
            Data(floor(GenCounter/10)+1).EP = EP; % reminder: start idx form 0
            Data(floor(GenCounter/10)+1).population = population;
            Data(floor(GenCounter/10)+1).IT = IT;
        end
        save('tempRun')
    end
    % save
    save(['Test',num2str(NSubProblems),'MT',num2str(MaxTopologies),'Pg',num2str(protectGen),'obj',num2str(ObjNum),'d',num2str(d>0),'Run',num2str(RunIdx)], 'EP','SolutionsToKeep', 'population','Senario','TopologyTracker','Data','MaxTopologies','protectGen','LayerRes');
    fprintf('\nRun %d/%d\n',RunIdx,TotalRuns)
    RunIdx = RunIdx + 1;
    gen0 = 1;





end
warning('on','all')

