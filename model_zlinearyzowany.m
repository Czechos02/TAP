clear; clc; close all;

%  ==========================================
% PARAMETRY 
% ===========================================
q = 10^6;
q_c = 10^6;
c_p = 1;
c_pc = 1;
k0 = 10^10;
E_R = 8330.1;
h = 130*10^6;
a = 1.678*10^6;
b = 0.5;

% punkt pracy
F = 1;
F_in = 1;
V = 1;

C_A_0 = 0.26;
T_0 = 393.9;

p.p11 = F_in/V;
p.p12 = -F/V;
p.p13 = -k0;
p.p14 = -E_R;
p.p21 = F_in/V;
p.p22 = -F/V;
p.p23 = h*k0/(q*c_p);
p.p24 = -E_R;
p.p25 = -V*q*c_p/a;
p.p26 = -V/2;
p.b = 0.5;
p.C_A_0 = C_A_0;
p.T_0 = T_0;
p.F_C0 = 15;
p.T_Cin0 = 365;

p.A = [p.p12 + p.p13*exp(p.p14/p.T_0), (-p.p14/p.T_0^2)*p.p13*p.C_A_0*exp(p.p14/p.T_0);
    p.p23*exp(p.p24/p.T_0), p.p22+ (-p.p24/p.T_0^2)*p.p23*p.C_A_0*exp(p.p24/p.T_0)+p.F_C0/(p.p25*p.F_C0^(1-p.b)+p.p26)];

p.B = [p.p11, 0, 0, 0;
    0, (p.T_0-p.T_Cin0)*(p.p25*p.F_C0^(1-b)+p.p26-(1-b)*p.F_C0*p.p25*p.F_C0^(-b))/((p.p25*p.F_C0^(1-p.b)+p.p26)^2), p.p21,(-p.F_C0)/(p.p25*p.F_C0^(1-b)+p.p26)];

% ==============================================
% WARUNKI POCZĄTKOWE
% ==============================================

% y(1) -> C_A
% y(2) -> T
% u(1) -> C_Ain
% u(2) -> F_C
% u(3) -> T_in
% u(4) -> T_Cin

y0   = [0.26 393.9];
u0   = [2 15 323 365]; 
sim_end = 20;
sim_time = [0, sim_end];      

% ==============================================
%  SYMULACJA
% ==============================================

u1_step = 0.5;
u2_step = 0;
z1_step = 0;
z2_step = 0;
step_time = 10;

u1 = @(t) step_function(t, u0(1), u0(1)+u1_step, 10);
u2 = @(t) step_function(t, u0(2), u0(2)+u2_step, 10);
z1 = @(t) step_function(t, u0(3), u0(3)+z1_step, 10);
z2 = @(t) step_function(t, u0(4), u0(4)+z2_step, 10);

u1_lin = @(t) step_function(t, 0, u1_step, 10);
u2_lin = @(t) step_function(t, 0, u2_step, 10);
z1_lin = @(t) step_function(t, 0, z1_step, 10);
z2_lin = @(t) step_function(t, 0, z2_step, 10);

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.01);
[t, Y] = ode15s(@(t, y) model_nl(t, y, u1, u2, z1, z2, p), sim_time, y0, opts);
[t_lin, Y_lin] = ode15s(@(t, y) model_lin(t, y, u1_lin, u2_lin, z1_lin, z2_lin, p), sim_time, [0 0], opts);

Y_lin = Y_lin + y0;

y1 = Y(:,1);             % np. y1 = x1
y2 = Y(:,2);             % np. y2 = x2

y1_lin = Y_lin(:,1);             % np. y1 = x1
y2_lin = Y_lin(:,2);             % np. y2 = x2

U1 = arrayfun(u1_lin, t);
U2 = arrayfun(u2_lin, t);
Z1 = arrayfun(z1_lin, t);
Z2 = arrayfun(z2_lin, t);

% ==============================================
% SYMULACJA TRANSMITANCJI MIMO
% ==============================================
s = tf('s');

G_u = tf(ss(p.A,p.B,eye(2),0));   % transmitancja dla wejść sterujących

t_tr = 0:0.01:sim_end;

U1_tr = arrayfun(u1_lin, t_tr);
U2_tr = arrayfun(u2_lin, t_tr);
Z1_tr = arrayfun(z1_lin, t_tr);
Z2_tr = arrayfun(z2_lin, t_tr);
% odpowiedź na wymuszenia

Y_u = lsim(G_u,[U1_tr;U2_tr;Z1_tr;Z2_tr],t_tr);

% wyraz wolny + punkt pracy
Y_tr = Y_u + y0;

% ==============================================
% WYKRESY
% ==============================================

figure('Name','Symulacja nieliniowa 2×2','NumberTitle','off', ...
       'Position',[100 100 900 700]);

subplot(3,2,1);
plot(t, y1, 'b', 'LineWidth', 1.5);
hold on
plot(t_lin, y1_lin, 'r', 'LineWidth', 1.5);
plot(t_tr, Y_tr(:,1),'g--','LineWidth',1.5);
hold off
xlabel('t [s]'); ylabel('y_1'); title('Wyjście y_1');
grid on;

subplot(3,2,2);
plot(t, y2, 'c', 'LineWidth', 1.5);
hold on
plot(t_lin, y2_lin, 'r', 'LineWidth', 1.5);
plot(t_tr, Y_tr(:,2),'g--','LineWidth',1.5);
hold off
xlabel('t [s]'); ylabel('y_2'); title('Wyjście y_2');
grid on;

subplot(3,2,3);
plot(t, U1, 'y', 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('u_1'); title('Sterowanie u_1(t)');
grid on;

subplot(3,2,4);
plot(t, U2, 'y', 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('u_2'); title('Sterowanie u_2(t)');
grid on;

subplot(3,2,5);
plot(t, Z1, 'w', 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('z_1'); title('Zakłócenie z_1(t)');
grid on;

subplot(3,2,6);
plot(t, Z2, 'w', 'LineWidth', 1.5);
xlabel('t [s]'); ylabel('z_2'); title('Zakłócenie z_2(t)');
grid on;

sgtitle('Symulacja modelu nieliniowego 2 wejścia / 2 wyjścia');

% ==============================================
%  FUNKCJA
% ==============================================
function dy_dt = model_nl(t, y, u1, u2, z1, z2, p)

C_Ain = u1(t);
F_C = u2(t);
T_in = z1(t);
T_Cin = z2(t);

dy_dt = [
    p.p11 * C_Ain ...
    + p.p12 * y(1) ...
    + p.p13*exp(p.p14/y(2))*y(1);

    p.p21*T_in ...
    + p.p22*y(2) ...
    + p.p23*exp(p.p24/y(2))*y(1) ...
    + (y(2)-T_Cin)*F_C/(p.p25*((F_C)^(1-p.b))+p.p26); 
];

end

function dy_dt = model_lin(t, y, u1, u2, z1, z2, p)

C_Ain = u1(t);
F_C = u2(t);
T_in = z1(t);
T_Cin = z2(t);
C_A = y(1);
T = y(2);

dy_dt = p.A*[C_A; T]+p.B*[C_Ain;F_C;T_in;T_Cin];

end


function out = step_function(t, before, after, step_time)
    if t < step_time
        out = before; 
    else
        out = after; 
    end
end
