%% This program analyze the sensetivity of the PPN guidance law

clc
clear
Wt = 10^4;
g = 9.81;
rth = 10;
fig = 0;
Pen = 1;
plotF = 0;
SenNum = [4:9] ;
%X = [Xm;Ym;Zm;phiM;gammaM;Vm,Xt;Yt;Zt;phiT;gammaT;integral(am^2);integral(Wt)]
NtoSave = [];
d = 0.005;
pltflg = 0;
tspan = 0:0.01:50;
Senario(1).X0 = [0 0 3000 30*pi/180 -125*pi/180  400 3000 0 3000 0 0 200 0 0];
Senario(1).TarA = 1;
Senario(1).d = d;
Senario(1).tspan = tspan;


Senario(2).X0 = [0 0 3000 30*pi/180 -125*pi/180  400 3000 0 3000 0 0 200 0 0];
Senario(2).TarA = 2;
Senario(2).d = d;
Senario(2).tspan = tspan;

Senario(3).X0 = [0 0 3000 30*pi/180 -125*pi/180  400 3000 0 3000 0 0 200 0 0];
Senario(3).TarA = 3;
Senario(3).d = d;
Senario(3).tspan = tspan;

Senario(4).X0 = [0 0 3000 45*pi/180 -45*pi/180 400 9000 3000 3000 180*pi/180 0 200 0 0];
Senario(4).TarA = 1;
Senario(4).d = d;
Senario(4).tspan = tspan;

Senario(5).X0 = [0 0 3000 45*pi/180 -45*pi/180 400 9000 3000 3000 180*pi/180 0 200 0 0];
Senario(5).TarA = 2;
Senario(5).d = d;
Senario(5).tspan = tspan;

Senario(6).X0 = [0 0 3000 45*pi/180 -45*pi/180 400 9000 3000 3000 180*pi/180 0 200 0 0];
Senario(6).TarA = 3;
Senario(6).d = d;
Senario(6).tspan = tspan;

Senario(7).X0 = [0 0 3000 -75*pi/180 -70*pi/180 400 3000 0 3000 0*pi/180 0 200 0 0];
Senario(7).TarA = 1;
Senario(7).d = 0.05;
Senario(7).tspan = tspan;

Senario(8).X0 = [0 0 3000 -75*pi/180 -70*pi/180 400 3000 0 3000 0*pi/180 0 200 0 0];
Senario(8).TarA = 2;
Senario(8).d = d;
Senario(8).tspan = tspan;

Senario(9).X0 = [0 0 3000 -75*pi/180 -70*pi/180 400 3000 0 3000 0*pi/180 0 200 0 0];
Senario(9).TarA = 3;
Senario(9).d = d;
Senario(9).tspan = tspan;
Senario = Senario(SenNum);

Nspan = [2:0.01:6];

for N = Nspan
    disp(N)
    for j = 1 : length(Senario)
        X0 = Senario(j).X0;
        TarA = Senario(j).TarA;
        d = Senario(j).d;
        tspan = Senario(j).tspan;
        [tout, X] = ode45(@(t, x) PEsimulationNC(t, x, N, TarA, [], d, Wt, g, rth), tspan, X0', Opt);
        r = sqrt((X(:,1)-X(:,7)).^2+(X(:,2)-X(:,8)).^2+(X(:,3)-X(:,9)).^2);
        rmin(j) = min(r);
        if rmin(j) < rth
            timeidx = find(r < rth, 1, "first");
            T(j) = tout(timeidx);
            J2(j) = X(timeidx,13);
        else
            T(j) = -100;
            J2(j) = -100;
        end
    end
    disp(rmin)
    if all(rmin<=rth)
        disp('V')
        Ind.N = N;
        Ind.T = T';
        Ind.Ja = J2';
        NtoSave = [NtoSave,Ind];
    end
end


if isempty(NtoSave)
    disp('is empty')
else
    figure(1);
    parallelcoords([NtoSave(:).T]')
    grid on
    title('Time')
    figure(2)
    parallelcoords([NtoSave(:).Ja]')
    grid on
    title('Energy')
end
save('PPNResultswithoutD','NtoSave')