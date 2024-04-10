function offspring = esrep(parent,mutation,GenCounter)
yl = -1*mutation.weight_cap;
yu = mutation.weight_cap;
offspring = parent;
vecsize = size(parent.connectiongenes,2);
pmutat = mutation.pol(1);
sigma0 = mutation.adaptiverate(1);
rate = mutation.adaptiverate(2);
    for j =  1:vecsize
        if rand(1) <= pmutat
            y = parent.connectiongenes(4,j);
            rang = -1 + 2*rand(1);
            sigma = sigma0*exp(-rate*GenCounter/200);
            neW = y + rang*sigma;
            if neW > yu
               neW = yu;
            elseif neW < yl
                neW = yl;
            end
            offspring.connectiongenes(4,j) = neW;
            
        end
    end
end
            