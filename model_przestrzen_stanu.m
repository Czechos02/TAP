clear; clc; close all;
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end

%% Punkt pracy
[y0, u0, p] = punkt_pracy();

%% Model liniowy ciagly
A = jacobian_A(y0,u0,p);
B = jacobian_B(y0,u0,p);
C_mat = eye(2);
D_mat = zeros(2,4);

sys_c = ss(A, B, C_mat, D_mat);

fprintf('========================================\n');
fprintf('  MODEL CIAGLY\n');
fprintf('========================================\n');
fprintf('A =\n'); disp(A);
fprintf('B =\n'); disp(B);

eigA = eig(A);
fprintf('Wartosci wlasne A:\n');
fprintf('  l1 = %.4f + %.4fi\n', real(eigA(1)), imag(eigA(1)));
fprintf('  l2 = %.4f + %.4fi\n\n', real(eigA(2)), imag(eigA(2)));

fprintf('Transmitancje ciagle G(s):\n');
G_c = tf(sys_c);
G_c

%% Odpowiedzi skokowe - przestrzen stanow
t_sim = linspace(0, 20, 5000);
in_names  = {'C_{Ain}', 'F_C', 'T_{in}', 'T_{Cin}'};
colors_in = {'b', [0 0.6 0], 'r', [0.6 0 0.6]};

figure('Name','1a - Przestrzen stanow','NumberTitle','off', ...
    'Position',[50 400 1000 500]);

for in_i = 1:4
    du_i = zeros(length(t_sim), 4);
    du_i(:, in_i) = 1;
    [y_i, ~] = lsim(sys_c, du_i, t_sim);

    subplot(1,2,1); hold on;
    plot(t_sim, y_i(:,1), 'Color', colors_in{in_i}, 'LineWidth', 1.5);

    subplot(1,2,2); hold on;
    plot(t_sim, y_i(:,2), 'Color', colors_in{in_i}, 'LineWidth', 1.5);
end

leg_ss = cellfun(@(s) ['skok ', s, ' = +1'], in_names, 'UniformOutput', false);

subplot(1,2,1);
xlabel('t [min]'); ylabel('\Delta C_A'); title('Wyjscie C_A');
legend(leg_ss, 'Location','best', 'FontSize', 8); grid on;

subplot(1,2,2);
xlabel('t [min]'); ylabel('\Delta T'); title('Wyjscie T');
legend(leg_ss, 'Location','best', 'FontSize', 8); grid on;

sgtitle('Model liniowy - przestrzen stanow (punkt 1a)');
exportgraphics(gcf, 'wykresy/1a_przestrzen_stanow.pdf', 'ContentType', 'vector');

%% Odpowiedzi skokowe - transmitancje
figure('Name','1a - Transmitancje','NumberTitle','off', ...
    'Position',[50 50 1000 500]);

u_step = ones(length(t_sim), 1);

for in_i = 1:4
    for out_i = 1:2
        subplot(1,2,out_i); hold on;
        y_tf = lsim(G_c(out_i, in_i), u_step, t_sim);
        plot(t_sim, y_tf, 'Color', colors_in{in_i}, 'LineWidth', 1.5);
    end
end

subplot(1,2,1);
xlabel('t [min]'); ylabel('\Delta C_A'); title('Wyjscie C_A');
legend(leg_ss, 'Location','best', 'FontSize', 8); grid on;

subplot(1,2,2);
xlabel('t [min]'); ylabel('\Delta T'); title('Wyjscie T');
legend(leg_ss, 'Location','best', 'FontSize', 8); grid on;

sgtitle('Model liniowy - transmitancje G(s) (punkt 1a)');
exportgraphics(gcf, 'wykresy/1a_transmitancje.pdf', 'ContentType', 'vector');

%% Jakosc dyskretyzacji - rozne Tp
Tp_test = [0.01, 0.05, 0.1, 0.5];
t_fine = linspace(0, 20, 5000);
colors_tp = {'r--', 'g-.', 'm-', 'k--'};

% Skok CAin = +1
du_step = zeros(length(t_fine), 4);
du_step(:, 1) = 1;
[y_c_step, ~] = lsim(sys_c, du_step, t_fine);

figure('Name','Jakosc dyskretyzacji - skok CAin','NumberTitle','off', ...
    'Position',[50 50 1000 500]);

for i = 1:length(Tp_test)
    Tp = Tp_test(i);
    sys_d_i = c2d(sys_c, Tp, 'zoh');
    t_d_i = (0:Tp:20)';
    du_d_i = zeros(length(t_d_i), 4);
    du_d_i(:, 1) = 1;
    [y_d_i, ~] = lsim(sys_d_i, du_d_i, t_d_i);

    subplot(1,2,1); hold on;
    stairs(t_d_i, y_d_i(:,1), colors_tp{i}, 'LineWidth', 1.2);
    subplot(1,2,2); hold on;
    stairs(t_d_i, y_d_i(:,2), colors_tp{i}, 'LineWidth', 1.2);
end

subplot(1,2,1);
plot(t_fine, y_c_step(:,1), 'b-', 'LineWidth', 2.5);
xlabel('t [min]'); ylabel('\Delta C_A'); title('C_A');
legend('Tp=0.01','Tp=0.05','Tp=0.1','Tp=0.5','ciagly','Location','best');
grid on;

subplot(1,2,2);
plot(t_fine, y_c_step(:,2), 'b-', 'LineWidth', 2.5);
xlabel('t [min]'); ylabel('\Delta T'); title('T');
legend('Tp=0.01','Tp=0.05','Tp=0.1','Tp=0.5','ciagly','Location','best');
grid on;

sgtitle('Jakosc dyskretyzacji (ZOH) - skok C_{Ain} = +1');
exportgraphics(gcf, 'wykresy/1d_dyskretyzacja_CAin.pdf', 'ContentType', 'vector');

% Skok FC = +1
du_step2 = zeros(length(t_fine), 4);
du_step2(:, 2) = 1;
[y_c_step2, ~] = lsim(sys_c, du_step2, t_fine);

figure('Name','Jakosc dyskretyzacji - skok FC','NumberTitle','off', ...
    'Position',[100 100 1000 500]);

for i = 1:length(Tp_test)
    Tp = Tp_test(i);
    sys_d_i = c2d(sys_c, Tp, 'zoh');
    t_d_i = (0:Tp:20)';
    du_d_i = zeros(length(t_d_i), 4);
    du_d_i(:, 2) = 1;
    [y_d_i, ~] = lsim(sys_d_i, du_d_i, t_d_i);

    subplot(1,2,1); hold on;
    stairs(t_d_i, y_d_i(:,1), colors_tp{i}, 'LineWidth', 1.2);
    subplot(1,2,2); hold on;
    stairs(t_d_i, y_d_i(:,2), colors_tp{i}, 'LineWidth', 1.2);
end

subplot(1,2,1);
plot(t_fine, y_c_step2(:,1), 'b-', 'LineWidth', 2.5);
xlabel('t [min]'); ylabel('\Delta C_A'); title('C_A');
legend('Tp=0.01','Tp=0.05','Tp=0.1','Tp=0.5','ciagly','Location','best');
grid on;

subplot(1,2,2);
plot(t_fine, y_c_step2(:,2), 'b-', 'LineWidth', 2.5);
xlabel('t [min]'); ylabel('\Delta T'); title('T');
legend('Tp=0.01','Tp=0.05','Tp=0.1','Tp=0.5','ciagly','Location','best');
grid on;

sgtitle('Jakosc dyskretyzacji (ZOH) - skok F_C = +1');
exportgraphics(gcf, 'wykresy/1d_dyskretyzacja_FC.pdf', 'ContentType', 'vector');

%% Model dyskretny Tp = 0.1
Ts = 0.1;
sys_d = c2d(sys_c, Ts, 'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sys_d);

fprintf('========================================\n');
fprintf('  MODEL DYSKRETNY (Tp = %.2f min)\n', Ts);
fprintf('========================================\n');
fprintf('Ad =\n'); disp(Ad);
fprintf('Bd =\n'); disp(Bd);

eig_d = eig(Ad);
fprintf('Wartosci wlasne Ad:\n');
fprintf('  z1 = %.6f + %.6fi  (|z1| = %.6f)\n', real(eig_d(1)), imag(eig_d(1)), abs(eig_d(1)));
fprintf('  z2 = %.6f + %.6fi  (|z2| = %.6f)\n', real(eig_d(2)), imag(eig_d(2)), abs(eig_d(2)));

fprintf('\nTransmitancje dyskretne:\n');
Gz = tf(sys_d)

%% Symulacja modelu dyskretnego
t_dysk = (0:Ts:20)';
du_dysk = zeros(length(t_dysk),4);
step_time = 10;
du_dysk(t_dysk>=step_time,1) = 1;

[y_d, t_d] = lsim(sys_d, du_dysk, t_dysk);

figure('Name','Model dyskretny','NumberTitle','off', ...
    'Position',[100 100 900 400]);

subplot(1,2,1);
stairs(t_d, y_d(:,1) + y0(1), 'b', 'LineWidth', 1.5);
xlabel('t [min]'); ylabel('C_A [kmol/m^3]'); title('C_A (dyskretny)');
grid on;

subplot(1,2,2);
stairs(t_d, y_d(:,2) + y0(2), 'r', 'LineWidth', 1.5);
xlabel('t [min]'); ylabel('T [K]'); title('T (dyskretny)');
grid on;

sgtitle(sprintf('Symulacja modelu dyskretnego (Tp = %.2f min, lsim)', Ts));
exportgraphics(gcf, 'wykresy/1d_model_dyskretny.pdf', 'ContentType', 'vector');

%%

function A = jacobian_A(x,u,p)
x1 = x(1);
x2 = x(2);
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
x1 = x(1);
x2 = x(2);
u1 = u(1);
u2 = u(2);
z1 = u(3);
z2 = u(4);

exp_term = exp(-p.E_R/x2);
K = (p.a*u2^(p.b+1)) / (u2 + (p.a*u2^p.b)/(2*p.q*p.c_pc));

N  = p.a * u2^(p.b+1);
D  = u2 + (p.a*u2^p.b)/(2*p.q*p.c_pc);
dN = p.a * (p.b+1) * u2^p.b;
dD = 1 + (p.a*p.b*u2^(p.b-1))/(2*p.q*p.c_pc);
dK = (dN*D - N*dD) / (D^2);

B = zeros(2,4);
B(1,1) = p.F_in / p.V;
B(1,2) = 0;
B(1,3) = 0;
B(1,4) = 0;
B(2,1) = 0;
B(2,2) = -(x2 - z2)/(p.V*p.q*p.c_p) * dK;
B(2,3) = p.F_in / p.V;
B(2,4) = K / (p.V*p.q*p.c_p);
end
