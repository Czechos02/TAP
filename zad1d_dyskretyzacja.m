clear; clc; close all;

% ==========================================================
% PUNKT 1d: Modele liniowe dyskretne
%   - Rownania stanu i transmitancje dyskretne
%   - Porownanie jakosci dyskretyzacji (rozne Tp)
%   - Implementacja modelu dyskretnego w petli
% ==========================================================

%% PARAMETRY
p.q = 1e6;  p.q_c = 1e6;
p.c_p = 1;  p.c_pc = 1;
p.k0 = 1e10; p.E_R = 8330.1;
p.h = 130e6; p.a = 1.678e6; p.b = 0.5;
p.V = 1; p.F = 1; p.F_in = 1;

%% PUNKT PRACY
x0 = [0.26, 393.9];
u0 = [2, 15, 323, 365];

%% MODEL CIAGLY (linearyzacja)
A = jacobian_A(x0, u0, p);
B = jacobian_B(x0, u0, p);
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

%% ============================================================
%  DYSKRETYZACJA - rozne okresy probkowania
% =============================================================

Tp_test = [0.01, 0.05, 0.1, 0.5];

for i = 1:length(Tp_test)
    Tp = Tp_test(i);
    sys_d = c2d(sys_c, Tp, 'zoh');
    [Ad, Bd, Cd, Dd] = ssdata(sys_d);

    fprintf('\n======= MODEL DYSKRETNY Tp = %.2f min =======\n', Tp);
    fprintf('Ad =\n'); disp(Ad);
    fprintf('Bd =\n'); disp(Bd);

    eig_d = eig(Ad);
    fprintf('Wartosci wlasne Ad:\n');
    fprintf('  z1 = %.6f + %.6fi  (|z1| = %.6f)\n', ...
        real(eig_d(1)), imag(eig_d(1)), abs(eig_d(1)));
    fprintf('  z2 = %.6f + %.6fi  (|z2| = %.6f)\n', ...
        real(eig_d(2)), imag(eig_d(2)), abs(eig_d(2)));

    fprintf('\nTransmitancje dyskretne G(z):\n');
    G_d = tf(sys_d);
    G_d
end

%% ============================================================
%  POROWNANIE JAKOSCI DYSKRETYZACJI
%  Odpowiedzi skokowe: ciagly vs dyskretny (rozne Tp)
% =============================================================

t_end = 20;
t_fine = linspace(0, t_end, 5000);

in_names  = {'C_{Ain}', 'F_C', 'T_{in}', 'T_{Cin}'};
out_names = {'C_A', 'T'};
colors    = {'r--', 'g-.', 'm-', 'k--'};

% --- Sterowania (CAin, FC) ---
figure('Name','Jakosc dyskretyzacji - sterowania', ...
    'NumberTitle','off', 'Position',[50 50 1200 700]);

idx = 0;
for in_i = 1:2  % CAin, FC
    % Odpowiedz ciagla
    du_c = zeros(length(t_fine), 4);
    du_c(:, in_i) = 1;
    [y_c, ~] = lsim(sys_c, du_c, t_fine);

    for out_i = 1:2
        idx = idx + 1;
        subplot(2, 2, idx);
        plot(t_fine, y_c(:, out_i), 'b-', 'LineWidth', 2.5); hold on;

        leg = {'ciagly'};
        for i = 1:length(Tp_test)
            Tp = Tp_test(i);
            sys_d = c2d(sys_c, Tp, 'zoh');
            t_d = (0:Tp:t_end)';
            du_d = zeros(length(t_d), 4);
            du_d(:, in_i) = 1;
            [y_d, ~] = lsim(sys_d, du_d, t_d);
            stairs(t_d, y_d(:, out_i), colors{i}, 'LineWidth', 1.2);
            leg{end+1} = sprintf('Tp=%.2f min', Tp);
        end

        grid on;
        title(sprintf('%s -> %s', in_names{in_i}, out_names{out_i}));
        xlabel('t [min]');
        ylabel(sprintf('\\Delta %s', out_names{out_i}));
        legend(leg, 'Location','best', 'FontSize', 8);
    end
end
sgtitle('Jakosc dyskretyzacji (ZOH) - sterowania');

% --- Zaklocenia (Tin, TCin) ---
figure('Name','Jakosc dyskretyzacji - zaklocenia', ...
    'NumberTitle','off', 'Position',[100 100 1200 700]);

idx = 0;
for in_i = 3:4  % Tin, TCin
    du_c = zeros(length(t_fine), 4);
    du_c(:, in_i) = 1;
    [y_c, ~] = lsim(sys_c, du_c, t_fine);

    for out_i = 1:2
        idx = idx + 1;
        subplot(2, 2, idx);
        plot(t_fine, y_c(:, out_i), 'b-', 'LineWidth', 2.5); hold on;

        leg = {'ciagly'};
        for i = 1:length(Tp_test)
            Tp = Tp_test(i);
            sys_d = c2d(sys_c, Tp, 'zoh');
            t_d = (0:Tp:t_end)';
            du_d = zeros(length(t_d), 4);
            du_d(:, in_i) = 1;
            [y_d, ~] = lsim(sys_d, du_d, t_d);
            stairs(t_d, y_d(:, out_i), colors{i}, 'LineWidth', 1.2);
            leg{end+1} = sprintf('Tp=%.2f min', Tp);
        end

        grid on;
        title(sprintf('%s -> %s', in_names{in_i}, out_names{out_i}));
        xlabel('t [min]');
        ylabel(sprintf('\\Delta %s', out_names{out_i}));
        legend(leg, 'Location','best', 'FontSize', 8);
    end
end
sgtitle('Jakosc dyskretyzacji (ZOH) - zaklocenia');

%% ============================================================
%  WYBRANY MODEL DYSKRETNY (Tp = 0.1 min)
% =============================================================

Tp = 0.1;
sys_d_final = c2d(sys_c, Tp, 'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sys_d_final);

fprintf('\n========================================\n');
fprintf('  WYBRANY MODEL DYSKRETNY (Tp = %.2f min)\n', Tp);
fprintf('========================================\n');
fprintf('Ad =\n'); disp(Ad);
fprintf('Bd =\n'); disp(Bd);
fprintf('Cd = I(2x2),  Dd = 0(2x4)\n\n');

fprintf('Transmitancje dyskretne (Tp = %.2f min):\n', Tp);
G_d_final = tf(sys_d_final);
G_d_final

%% ============================================================
%  WERYFIKACJA: dyskretny (lsim) vs nieliniowy (ode15s)
% =============================================================

t_d_lsim = (0:Tp:t_end)';
du_lsim = zeros(length(t_d_lsim), 4);
du_lsim(50:end, 1) = 0.1;   % skok dCAin=+0.1 w k=50
du_lsim(100:end, 2) = 1.0;  % skok dFC=+1 w k=100
[y_lsim, ~] = lsim(sys_d_final, du_lsim, t_d_lsim);

% Nieliniowy z tym samym scenariuszem
t_step1 = 50 * Tp;
t_step2 = 100 * Tp;
u_funs = {@(t) u0(1) + 0.1*(t>=t_step1), ...
          @(t) u0(2) + 1.0*(t>=t_step2), ...
          @(t) u0(3), @(t) u0(4)};

opts_nl = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.01);
[t_nl, Y_nl] = ode15s(@(t,y) model_nl(t, y, ...
    u_funs{1}, u_funs{2}, u_funs{3}, u_funs{4}, p), ...
    [0, t_end], x0, opts_nl);

figure('Name','Weryfikacja: dyskretny vs nieliniowy', ...
    'NumberTitle','off', 'Position',[150 50 1000 500]);

subplot(1,2,1);
plot(t_nl, Y_nl(:,1), 'b-', 'LineWidth', 2); hold on;
stairs(t_d_lsim, y_lsim(:,1) + x0(1), 'r--', 'LineWidth', 1.5);
xlabel('t [min]'); ylabel('C_A [kmol/m^3]'); title('C_A');
legend('nieliniowy', 'dyskretny (lsim)', 'Location','best');
grid on;

subplot(1,2,2);
plot(t_nl, Y_nl(:,2), 'b-', 'LineWidth', 2); hold on;
stairs(t_d_lsim, y_lsim(:,2) + x0(2), 'r--', 'LineWidth', 1.5);
xlabel('t [min]'); ylabel('T [K]'); title('T');
legend('nieliniowy', 'dyskretny (lsim)', 'Location','best');
grid on;

sgtitle(sprintf('Weryfikacja: nieliniowy vs dyskretny (Tp = %.2f min)', Tp));
exportgraphics(gcf, 'wykresy/1d_weryfikacja.pdf', 'ContentType', 'vector');

%% ============================================================
%  FUNKCJE LOKALNE
% =============================================================

function dy = model_nl(t, y, u1, u2, z1, z2, p)
    C_Ain = u1(t); F_C = u2(t); T_in = z1(t); T_Cin = z2(t);
    C_A = y(1); T = y(2);
    exp_term = exp(-p.E_R / T);
    K = (p.a * F_C^(p.b+1)) / (F_C + (p.a * F_C^p.b)/(2*p.q_c*p.c_pc));
    dy = [
        p.F_in*C_Ain/p.V - p.F*C_A/p.V - p.k0*exp_term*C_A;
        p.F_in*T_in/p.V - p.F*T/p.V + (p.h*p.k0/(p.q*p.c_p))*exp_term*C_A ...
        - K*(T - T_Cin)/(p.V*p.q*p.c_p);
    ];
end

function A = jacobian_A(x, u, p)
    C_A = x(1); T = x(2); F_C = u(2);
    exp_term = exp(-p.E_R / T);
    K = (p.a*F_C^(p.b+1)) / (F_C + (p.a*F_C^p.b)/(2*p.q_c*p.c_pc));
    A = zeros(2,2);
    A(1,1) = -p.F/p.V - p.k0*exp_term;
    A(1,2) = -p.k0*C_A*exp_term*(p.E_R/T^2);
    A(2,1) = (p.h*p.k0/(p.q*p.c_p))*exp_term;
    A(2,2) = -p.F/p.V + (p.h*p.k0/(p.q*p.c_p))*C_A*exp_term*(p.E_R/T^2) ...
             - K/(p.V*p.q*p.c_p);
end

function B = jacobian_B(x, u, p)
    T = x(2); F_C = u(2); T_Cin = u(4);
    K = (p.a*F_C^(p.b+1)) / (F_C + (p.a*F_C^p.b)/(2*p.q_c*p.c_pc));
    N = p.a * F_C^(p.b+1);
    D = F_C + (p.a*F_C^p.b)/(2*p.q_c*p.c_pc);
    dN = p.a * (p.b+1) * F_C^p.b;
    dD = 1 + (p.a*p.b*F_C^(p.b-1))/(2*p.q_c*p.c_pc);
    dK = (dN*D - N*dD) / D^2;
    B = zeros(2,4);
    B(1,1) = p.F_in / p.V;
    B(2,2) = -(T - T_Cin) * dK / (p.V*p.q*p.c_p);
    B(2,3) = p.F_in / p.V;
    B(2,4) = K / (p.V*p.q*p.c_p);
end
