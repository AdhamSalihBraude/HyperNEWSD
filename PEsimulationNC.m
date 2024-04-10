function [dx] = PEsimulationNC(t, x, N, Tar, NC, d, Wt, g, rth)
%x = [xM;yM;zM;phiM;gammaM;VM,xT;yT;zT;phiT;gammaT;J1;J2]
global Nsave A
xm = x(1);
ym = x(2);
zm = x(3);
xt = x(7);
yt = x(8);
zt = x(9);
Dx = xt - xm;
Dy = yt - ym;
Dz = zt - zm;
r = sqrt(Dx^2+Dy^2+Dz^2);

if r > rth
    h = zm;
    phiM = x(4);

    gammaM = x(5);

    Vm = x(6);

    phiT = x(10);

    gammaT = x(11);
    Vt = x(12);
    VmG = [Vm*cos(gammaM)*cos(phiM);Vm*cos(gammaM)*sin(phiM);Vm*sin(gammaM)];
    VtG = Vt*[cos(gammaT)*cos(phiT);cos(gammaT)*sin(phiT);sin(gammaT)];
    if Tar == 2
        % bank-to-bank maneuver
        ayt = -1.52*g*(-1)^floor(t/6);
        apt = 1.52*g*(-1)^floor(t/6);
    elseif Tar == 3
        % step maneuvering
        ayt = -0.39*g;
        apt = 0.39*g;
    elseif Tar == 4
        % Sin-Cos
        apt = 1.52*g*cos(0.3*t);
        ayt = 1.52*g*sin(0.3*t);
    else
        % non maneuvering
        ayt = 0*g;
        apt = 0*g;
    end



    %% Control calc

    gammaL = atan2(Dz,(sqrt(Dy^2+Dx^2)));
    % r = sqrt(Dx^2+Dy^2+Dz^2);
    phiL = atan2(Dy,Dx);
    R1 = [cos(phiL),sin(phiL),0;-sin(phiL),cos(phiL),0;0,0,1];
    R2 = [cos(gammaL),0,sin(gammaL);0,1,0;-sin(gammaL),0,cos(gammaL)];
    VmL = R2*R1*VmG;
    VtL = R2*R1*VtG;
    Ldoty = (VmL(3)-VtL(3))/r;
    Ldotz = (VtL(2)-VmL(2))/r;
    if ~isempty(NC)
        Input = [[phiL; gammaL;phiM;gammaM]/(pi);1];
        for i = 2 : length(NC.LayerRes)
            W = NC.Layer(i).W;
            alpha = NC.Layer(i).alpha;
            Sig = alpha'.*W'*Input;
            Input = 1./(1+exp(-Sig));
        end
        N = 1.5 + (10-1.5)*Input;
        Nsave = [Nsave;[t,Input]];
    else
        Input = [0];
    end
    
    amL2 = [VmL(3)*Ldoty-VmL(2)*Ldotz;VmL(1)*Ldotz;-VmL(1)*Ldoty]*N;
    %% Netcalc

    RLG1 = [cos(-phiL),sin(-phiL),0;-sin(-phiL),cos(-phiL),0;0,0,1];
    RLG2 = [cos(-gammaL),0,sin(-gammaL);0,1,0;-sin(-gammaL),0,cos(-gammaL)];
    amG = RLG1*RLG2*amL2;
    RGm1 = [cos(phiM),sin(phiM),0;-sin(phiM),cos(phiM),0;0,0,1];
    RGm2 = [cos(gammaM),0,sin(gammaM);0,1,0;-sin(gammaM),0,cos(gammaM)];
    amM = RGm2*RGm1*amG;
    
    amax = 15*g;
    aym = amM(2);
    apm = amM(3);

    amax = 15*g/sqrt(2);
    aym = min([aym,amax]);
    aym = max([aym,-amax]);
    apm = min([apm,amax]);
    apm = max([apm,-amax]);
    
     am = sqrt(aym^2 + apm^2);
     
    % if ( am > amax)% && ~isempty(NC))
    %     aym = aym*(amax-0.001)/am;
    %     apm = apm*(amax-0.001)/am;
    %     am = sqrt(aym^2 + apm^2);
    % end
A = [A;[t,aym/g,apm/g]];
    %% updating the model parameters

    m = (165-4*t)*(t<=10) + (125-(t-10))*(t>10 && t <=40 ) + (95)*(t>40);
    Lambda = 12000*(t<=10) + 2000*(t>10 && t <=40 );



    R = 288;
    T = (288.16-0.0065*h)*(h<11000) + (216.66-0.0065*h)*(h>11000);

    M = Vm/sqrt(1.4*R*T);

    K = 0.2*(M<=1.15) + (0.2 + 0.246*(M-1.15))*(M>1.15);

    Cd0 = 0.02*(M<0.93) + (0.02+0.2*(M-0.93))*(M>=0.93 && M < 1.03) +...
        (0.04+0.06*(M-1.03))*(M>=1.03 && M < 1.1) + (0.0442+0.007*(M-1.03))*(M>=1.1);

    ru = 1.5579-1.058e-4*max([h,20000])+3.725e-9*(max([h,20000]))^2 - 6e-14*(max([h,20000]))^3;
    Q = 0.5*ru*Vm^2;
    if Q > 0
        s = 1;
        D0 = Cd0*Q*s;
        Di = K*m^2*sum((am).^2)^0.5/(Q*s);
    else
        D0 = 0;
        Di = 0;
    end
    D = D0 + Di;

    dis = d*randn(1,2);
    %% updating the model variables

    dx(1,1) = Vm*cos(gammaM)*cos(phiM);
    dx(1,2) = Vm*cos(gammaM)*sin(phiM);
    dx(1,3) = Vm*sin(gammaM);
    dx(1,4) = aym/(Vm*cos(gammaM)) + dis(1);
    dx(1,5) = (apm-g*cos(gammaM))/Vm - dis(2);
    dx(1,6) = (Lambda-D)/m-g*sin(gammaM);
    dx(1,7) = Vt*cos(gammaT)*cos(phiT);
    dx(1,8) = Vt*cos(gammaT)*sin(phiT);
    dx(1,9) = Vt*sin(gammaT);
    dx(1,10) = ayt/(Vt*cos(gammaT));
    dx(1,11) = apt/Vt;
    dx(1,12) = 0;
    dx(1,13) = am;
    dx(1,14) = Wt;
    dx = dx';
else
    dx = zeros(size(x));
end

end

