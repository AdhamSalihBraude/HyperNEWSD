function [value, isterminal, direction] = myEvent(T, x)
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

value      = (r < 10||zm<0);
isterminal = 1;   % Stop the integration
direction  = 0;
end