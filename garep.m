function [offspring1,offspring2] = garep(parent1,parent2,mutation,crossover)
offspring1 = parent1;
offspring2 = parent2;
yl = -1*mutation.weight_cap;
yu = mutation.weight_cap;
vecsize = size(parent1.connectiongenes,2);
pcross = crossover.SBXprop(1);
mu = crossover.SBXprop(2);
if rand(1) <= pcross
    for j =  1:vecsize
        par1 = parent1.connectiongenes(4,j);
        par2 = parent2.connectiongenes(4,j);
        
        rnd = rand(1);
        if rnd <= 0.5
            if abs(par1 - par2) > 0.000001
                if par2 > par1
                    y2 = par2;
                    y1 = par1;
                else
                    y2 = par1;
                    y1 = par2;
                end
                if (y1 - yl) > (yu - y2)
                    beta = 1 + (2*(yu - y2)/(y2 - y1));
                else
                    beta = 1 + (2*(y1 - yl)/(y2 - y1));
                end
                expp = mu + 1;
                beta = 1/beta;
                alpha = 2 - beta^expp;
                rnd = rand(1);
                if rnd <= 1/alpha
                    alpha = alpha*rnd;
                    expp = 1/(mu + 1);
                    betaq = alpha^expp;
                else
                    alpha = alpha*rnd;
                    alpha = 1/(2 - alpha);
                    expp = 1/(mu + 1);
                    betaq = alpha^expp;
                    
                end
                child1 = 0.5*((y1 + y2) - betaq*(y2 - y1));
                child2 = 0.5*((y1 + y2) + betaq*(y2 - y1));
            else
                betaq = 1;
                y1 = par1;
                y2 = par2;
                child1 = 0.5*((y1 + y2) - betaq*(y2 - y1));
                child2 = 0.5*((y1 + y2) + betaq*(y2 - y1));
            end
            if child1 < yl
                child1 = yl;
            end
            if child2 < yl
                child2 = yl;
            end
            if child1 > yu
                child1 = yu;
            end
            if child2 > yu
                child2 = yu;
            end
        else
            child1 = par1;
            child2 = par2;
        end
        if ~isreal(child1) || ~isreal(child2) 
            disp('error in ga rep sbx');
        end
        offspring1.connectiongenes(4,j) = child1;
        offspring2.connectiongenes(4,j) = child2;
    end
else
end



pmutat = mutation.pol(1);
mu = mutation.pol(2);
for index = 1:2
    if index == 1
        parent = offspring1;
    else
        parent = offspring2;
    end
    for j =  1:vecsize
        rnd = rand(1);
        if rnd <= pmutat
            y = parent.connectiongenes(4,j);
            if y > yl
                if (y - yl) < (yu - y)
                    delta = (y - yl)/(yu - yl);
                else
                    delta = (yu - y)/(yu - yl);
                end
                rnd = rand(1);
                indi = 1/(mu + 1);
                
                if rnd <= 0.5
                    u = 1 - delta;
                    val = 2*rnd + (1-2*rnd)*(u^(mu+1));
                    deltaq = val^indi - 1;
                else
                    u = 1 - delta;
                    val = 2*(1 - rnd) + 2*(rnd - 0.5)*(u^(mu + 1));
                    deltaq = 1 - val^indi;
                end
                y = y + deltaq * (yu - yl);
                if y > yu
                    y = yu;
                end
                if y < yl
                    y = yl;
                end
                if index == 1
                    offspring1.connectiongenes(4,j) = y;
                else
                    offspring2.connectiongenes(4,j) = y;
                end
            else
                if index == 1
                    offspring1.connectiongenes(4,j) = rand(1)*(yu - yl) + yl;
                else
                    offspring2.connectiongenes(4,j) = rand(1)*(yu - yl) + yl;
                end
            end
            if ~isreal(offspring1.connectiongenes(4,j)) 
            disp('error in ga rep polmutation');
            end
        end
    end
end
end

