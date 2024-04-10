function population_withNC = CreatNCs(population, LayerRes, Zres)
population_withNC = population;
for i = 1 : length(population)
    Zs = linspace(-1,1,Zres);
    NC = [];
    for L = 2 : length(LayerRes)
        z = Zs(L-1);
        inputX = linspace(-1,1,LayerRes(L-1));
        OutputX = linspace(-1,1,LayerRes(L));
        for x1idx = 1 : LayerRes(L-1)
            for x2idx = 1 : LayerRes(L)
                
                if sum(population(i).nodegenes(1,:)==1)==4
                    NC(i).Layer(L).W(x1idx,x2idx) = CPPNCalc(population(i),[inputX(x1idx),OutputX(x2idx),z,1]);
                    NC(i).Layer(L).alpha(x2idx) =  CPPNCalc(population(i),[inputX(x1idx),OutputX(x2idx),z,2]);
                else
                    NC(i).Layer(L).W(x1idx,x2idx) = CPPNCalc(population(i),[inputX(x1idx),OutputX(x2idx),z]);
                    NC(i).Layer(L).alpha(x2idx) =  0.5;
                end
            end
        end
    end
    NC(i).LayerRes = LayerRes;
    population_withNC(i).NC = NC(i);
end
end