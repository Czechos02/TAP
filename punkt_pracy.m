function [x0, u0, p] = punkt_pracy()

p.q = 1e6;  p.q_c = 1e6;
p.c_p = 1;  p.c_pc = 1;
p.k0 = 1e10; p.E_R = 8330.1;
p.h = 130e6; p.a = 1.678e6; p.b = 0.5;
p.V = 1; p.F = 1; p.F_in = 1;

u0 = [2, 15, 323, 365];

K = (p.a*u0(2)^(p.b+1)) / (u0(2) + (p.a*u0(2)^p.b)/(2*p.q_c*p.c_pc));

eq_fun = @(y) [
    p.F_in*u0(1)/p.V - p.F*y(1)/p.V - p.k0*exp(-p.E_R/y(2))*y(1);
    p.F_in*u0(3)/p.V - p.F*y(2)/p.V + (p.h*p.k0/(p.q*p.c_p))*exp(-p.E_R/y(2))*y(1) ...
    - K*(y(2)-u0(4))/(p.V*p.q*p.c_p)
];

opts = optimoptions('fsolve', 'Display','off', 'TolFun',1e-15, 'TolX',1e-15);
x0 = fsolve(eq_fun, [0.26, 393.9], opts);

end
