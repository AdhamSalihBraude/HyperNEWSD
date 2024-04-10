function population_with_Fit = evalua(population,Senario,plotF)
population_with_Fit = population;
N = 3;
Wt = 1;%2*10^5;
g = 9.81;
rth = 10;
fig = 0;
Pen = 1;

parfor i = 1 : length(population_with_Fit)
    
    Opt    = odeset('Events', @myEvent);
    J = nan(length(Senario),1);
    NC = population_with_Fit(i).NC;
    warning('off','all')
    for j = 1 : length(Senario)
        X0 = Senario(j).X0;
        TarA = Senario(j).TarA;
        d = Senario(j).d;
        tspan = Senario(j).tspan;
        [tout, X] = ode45(@(t, x) PEsimulationNC(t, x, N, TarA, NC, d, Wt, g, rth), tspan, X0', Opt);
        r = sqrt((X(:,1)-X(:,7)).^2+(X(:,2)-X(:,8)).^2+(X(:,3)-X(:,9)).^2);
        rmin = min(r);
        if rmin < rth
            rFinal = rmin;
            timeidx = find(r < rth, 1, "first");
        else
            rFinal = r(end);
            rmin = rFinal;
            timeidx = length(tout);
        end
        r0 = sqrt((X(1,1)-X(1,7))^2+(X(1,2)-X(1,8))^2+(X(1,3)-X(1,9))^2);

        h = min(X(:,3));
        J1 = X(timeidx,13);
        J2 = X(timeidx,14);
%         J(2*j-1,1) = rmin/rth + (J2)*(1 + (2+(rmin/rth)^2)*(rmin > rth || h < 0) );%/Senario(j).SatFit;
%         J(2*j,1) = rmin/rth + (J1)*(1 + (2+(rmin/rth)^2)*(rmin > rth || h < 0) )/Senario(j).SatFit;
        J(2*j-1,1) =  (J2)*(1 + (2+(rmin/rth)^2)*(rmin > rth || h < 0) );%/Senario(j).SatFit;
        J(2*j,1) = (J1)*(1 + (2+(rmin/rth)^2)*(rmin > rth || h < 0) )/Senario(j).SatFit;

        %       if plotF
        %           fig = fig + 1;
        %           figure(fig)
        %           plot3(X(1:timeidx,1),X(1:timeidx,2),X(1:timeidx,3),X(1:timeidx,7),X(1:timeidx,8),X(1:timeidx,9),'k',X(1,1),X(1,2),X(1,3),'ob')
        %           title(['Dist to Target  ',num2str(rmin)])
        %           grid on
        %           fig = fig + 1;
        %           figure(fig)
        %           plot(tout,X(:,13))
        %           title(['J = ',num2str(J(j,1))])
        %           grid on
        %       end
    end
    population_with_Fit(i).Fit = J;
end
end
