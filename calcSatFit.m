function [Updated_Senario,J] = calcSatFit(Senario,plotF)
Updated_Senario = Senario;
Nspan = 2:0.1:6;
Wt = 10^4;
g = 9.81;
rth = 10;
fig = 0;
Pen = 10e20;

J = nan(length(Nspan),length(Senario));
for j = 1 : length(Senario)
        Opt    = odeset('Events', @myEvent);
        X0 = Senario(j).X0;
        TarA = Senario(j).TarA;
        d = 0;
        tspan = Senario(j).tspan;
    
    for i = 1 : length(Nspan)
        N = Nspan(i);

        [tout, X] = ode45(@(t, x) PEsimulationNC(t, x, N, TarA, [], d, Wt, g, rth), tspan, X0', Opt);
        r = sqrt((X(:,1)-X(:,7)).^2+(X(:,2)-X(:,8)).^2+(X(:,3)-X(:,9)).^2);
        rmin = min(r);
        if rmin < rth
            rFinal = rmin;
            timeidx = find(r < rth, 1, "first");
            T(i,j) = tout(timeidx);
        else
            rFinal = r(end);
            timeidx = length(tout);
            T(i,j) = tspan(end);
        end
        r0 = sqrt((X(1,1)-X(1,7))^2+(X(1,2)-X(1,8))^2+(X(1,3)-X(1,9))^2);

        h = min(X(:,3));
        J1 = X(timeidx,13);
        J2 = X(timeidx,14);
        J(i,j) = rmin/rth + (J1 + J2)*(1 + (2+(rmin/rth)^2)*(rmin > rth || h < 0) );
% 
%         if plotF
%             fig = fig + 1;
%             figure(fig)
%             plot3(X(:,1),X(:,2),X(:,3),X(:,7),X(:,8),X(:,9),'k')
%             grid on
%             fig = fig + 1;
%             figure(fig)
%             plot(tout,X(:,13))
%             title(['J = ',num2str(J(j,i))])
%             grid on
%         end
    end
    [SatFit, Nidx] = min(J(:,j));
    Updated_Senario(j).SatFit = SatFit;
    Updated_Senario(j).N = Nspan(Nidx);
    
end
end
