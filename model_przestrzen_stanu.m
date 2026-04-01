clear; clc; close all;

% ==========================================
% PARAMETRY 
% ==========================================
p.q = 10^6;
p.q_c = 10^6;
p.c_p = 1;
p.c_pc = 1;
p.k0 = 10^10;
p.E_R = 8330.1;
p.h = 130*10^6;
p.a = 1.678*10^6;
p.b = 0.5;

% punkt pracy
p.F = 1;
p.F_in = 1;
p.V = 1;

% ==============================================
% WARUNKI POCZĄTKOWE
% ==============================================

y0   = [0.26 393.9]; % zmienne stanu (Ca, T) w pp
u0   = [2 15 323 365]; % sterowania (Cain, Fc) i zaklocenia (Tin, Tcin) w pp  

%% Model LTI %% ciagly DO SPRAWDZENIA

A = jacobian_A(y0,u0,p);
B = jacobian_B(y0,u0,p);
C = eye(2);
D = zeros(2,4);
t = linspace(0,20,1000);

du = zeros(length(t),4);
step_time = 10;
du(:,1) = (t >= step_time) * 1;     % skok CAin
du(:,2) = (t >= step_time) * 1;     % skok Fc
du(:,3) = (t >= step_time) * 1;     % skok Tin
du(:,4) = (t >= step_time) * 1;     % skok Tcin

sys_c = ss(A,B,C,D);

%% Model LTI %% dyskretny DO SPRAWDZENIA
Ts = 0.1;
sys_d = c2d(sys_c, Ts, 'zoh');
t_dysk = (0:Ts:20)';

du_dysk = zeros(length(t_dysk),4);
step_time = 10;
du_dysk(t_dysk>=step_time,1) = 1;
du_dysk(t_dysk>=step_time,2) = 1;
du_dysk(t_dysk>=step_time,3) = 1;
du_dysk(t_dysk>=step_time,4) = 1;

% symulacja (BEZ PĘTLI)
[y_d, t_d, x_d] = lsim(sys_d, du_dysk, t_dysk);

y1_d = y_d(:,1) + y0(1);
y2_d = y_d(:,2) + y0(2);

% Transmitancja dyskretna
Gz = tf(sys_d)

figure;
subplot(3,2,1);
stairs(t_d, y1_d, 'b','LineWidth',1.5);
title('C_A (dyskretny)');
grid on;

subplot(3,2,2);
stairs(t_d, y2_d, 'c','LineWidth',1.5);
title('T (dyskretny)');
grid on;

subplot(3,2,3);
plot(t_d, du_dysk(:,1), 'y', 'LineWidth', 1.5);
xlabel('T_dysk'); ylabel('u_1'); title('u_1');
grid on;

subplot(3,2,4);
plot(t_d, du_dysk(:,2), 'y', 'LineWidth', 1.5);
xlabel('T_dysk'); ylabel('u_2'); title('u_2');
grid on;

subplot(3,2,5);
plot(t_d, du_dysk(:,3), 'g', 'LineWidth', 1.5);
xlabel('T_dysk'); ylabel('z_1'); title('z_1');
grid on;

subplot(3,2,6);
plot(t_d, du_dysk(:,4), 'g', 'LineWidth', 1.5);
xlabel('T_dysk'); ylabel('z_2'); title('z_2');
grid on;

%% symulacja ciagla LTI %%

[y_lin, t, x_lin] = lsim(sys_c, du, t);
y1_lin = y_lin(:,1) + y0(1);
y2_lin = y_lin(:,2) + y0(2);

U1 = du(:,1) + u0(1);
U2 = du(:,2) + u0(2);
Z1 = du(:,3) + u0(3);
Z2 = du(:,4) + u0(4);

figure('Name','Model liniowy vs wejścia','NumberTitle','off', ...
       'Position',[100 100 900 700]);

subplot(3,2,1);
plot(t, y1_lin, 'b', 'LineWidth', 1.5);
xlabel('t'); ylabel('C_A'); title('C_A (model liniowy)');
grid on;

subplot(3,2,2);
plot(t, y2_lin, 'c', 'LineWidth', 1.5);
xlabel('t'); ylabel('T'); title('T (model liniowy)');
grid on;

subplot(3,2,3);
plot(t, U1, 'y', 'LineWidth', 1.5);
xlabel('t'); ylabel('C_{Ain}'); title('u_1');
grid on;

subplot(3,2,4);
plot(t, U2, 'y', 'LineWidth', 1.5);
xlabel('t'); ylabel('F_c'); title('u_2');
grid on;

subplot(3,2,5);
plot(t, Z1, 'g', 'LineWidth', 1.5);
xlabel('t'); ylabel('T_{in}'); title('z_1');
grid on;

subplot(3,2,6);
plot(t, Z2, 'g', 'LineWidth', 1.5);
xlabel('t'); ylabel('T_{cin}'); title('z_2');
grid on;

sgtitle('Symulacja modelu liniowego (LTI)');


% ==============================================
% JACOBIAN Macierzy A oraz B
% ==============================================
function A = jacobian_A(x,u,p)

x1 = x(1); % C_A
x2 = x(2); % T

u2 = u(2);

exp_term = exp(-p.E_R/x2);

K = (p.a*u2^(p.b+1)) / (u2 + (p.a*u2^p.b)/(2*p.q*p.c_pc));

A = zeros(2,2);

A(1,1) = -p.F/p.V - p.k0*exp_term;
A(1,2) = -p.k0*x1*exp_term*(p.E_R/x2^2);
A(2,1) = (p.h*p.k0/(p.q*p.c_p))*exp_term;
A(2,2) = -p.F/p.V ...
    + (p.h*p.k0/(p.q*p.c_p))*x1*exp_term*(p.E_R/x2^2) ...
    - K/(p.V*p.q*p.c_p);

end

function B = jacobian_B(x,u,p)

x1 = x(1); % C_A
x2 = x(2); % T

u1 = u(1);
u2 = u(2);

z1 = u(3);
z2 = u(4);

exp_term = exp(-p.E_R/x2);

% K (to samo co w modelu)
K = (p.a*u2^(p.b+1)) / (u2 + (p.a*u2^p.b)/(2*p.q*p.c_pc));

% Pochodna dK/du2
N  = p.a * u2^(p.b+1);
D  = u2 + (p.a*u2^p.b)/(2*p.q*p.c_pc);

dN = p.a * (p.b+1) * u2^p.b;
dD = 1 + (p.a*p.b*u2^(p.b-1))/(2*p.q*p.c_pc);
dK = (dN*D - N*dD) / (D^2);

% Macierz B (2x4: [u1 u2 u3 u4]) (u3 = z1, u4 = z2)
B = zeros(2,4);

% Pochodne z pierwszego rownania
B(1,1) = p.F_in / p.V;   % df1/du1
B(1,2) = 0;              % df1/du2
B(1,3) = 0;              % df1/dz1
B(1,4) = 0;              % df1/dz2

% Pochodne z drugiego rownania
B(2,1) = 0; % df2/du1
B(2,2) = -(x2 - z2)/(p.V*p.q*p.c_p) * dK; % df2/du2
B(2,3) = p.F_in / p.V; % df2/du3
B(2,4) = K / (p.V*p.q*p.c_p); % df2/du4
end