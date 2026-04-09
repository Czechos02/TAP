clear; clc; close all;

%% Punkt pracy
[x0, u0, p] = punkt_pracy();

%% Model ciagly
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

%% Dyskretyzacja - rozne Tp
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

%% Porownanie jakosci dyskretyzacji
t_end = 20;
t_fine = linspace(0, t_end, 5000);

in_names  = {'C_{Ain}', 'F_C', 'T_{in}', 'T_{Cin}'};
out_names = {'C_A', 'T'};
colors    = {'r--', 'g-.', 'm-', 'k--'};

figure('Name','Jakosc dyskretyzacji - sterowania', ...
    'NumberTitle','off', 'Position',[50 50 1200 700]);

idx = 0;
for in_i = 1:2
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
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end
exportgraphics(gcf, 'wykresy/1d_jakosc_sterowania.pdf', 'ContentType', 'vector');

figure('Name','Jakosc dyskretyzacji - zaklocenia', ...
    'NumberTitle','off', 'Position',[100 100 1200 700]);

idx = 0;
for in_i = 3:4
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
exportgraphics(gcf, 'wykresy/1d_jakosc_zaklocenia.pdf', 'ContentType', 'vector');

%% Wybrany model dyskretny Tp = 0.1
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

%% Weryfikacja: NL vs LIN vs dyskretny
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end

opts_nl = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.01);
step_time = 5;
t_lin = (0:0.01:t_end)';
t_d = (0:Tp:t_end)';

test_skoki = {1, +0.1; 2, +1};
input_names_w = {'C_{Ain}', 'F_C'};
input_short_w = {'CAin', 'FC'};

for ti = 1:size(test_skoki, 1)
    inp = test_skoki{ti, 1};
    amp = test_skoki{ti, 2};

    u_funs = {@(t) u0(1), @(t) u0(2), @(t) u0(3), @(t) u0(4)};
    u_funs{inp} = @(t) u0(inp) + amp*(t >= step_time);
    [t_nl, Y_nl] = ode15s(@(t,y) model_nl(t,y, ...
        u_funs{1},u_funs{2},u_funs{3},u_funs{4},p), ...
        [0 t_end], x0, opts_nl);

    du_c = zeros(length(t_lin), 4);
    du_c(t_lin >= step_time, inp) = amp;
    [y_lin, ~] = lsim(sys_c, du_c, t_lin);
    y_lin = y_lin + x0;

    du_d = zeros(length(t_d), 4);
    du_d(t_d >= step_time, inp) = amp;
    [y_disc, ~] = lsim(sys_d_final, du_d, t_d);
    y_disc = y_disc + x0;

    figure('Name', sprintf('1d weryfikacja %s=%+g', input_short_w{ti}, amp), ...
        'NumberTitle','off', 'Position',[50 50 900 400]);

    subplot(1,2,1);
    plot(t_nl, Y_nl(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t_lin, y_lin(:,1), 'r--', 'LineWidth', 1.5);
    stairs(t_d, y_disc(:,1), 'g-.', 'LineWidth', 1.2);
    xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
    xlabel('t [min]'); ylabel('C_A [kmol/m^3]'); title('C_A');
    legend('nieliniowy','liniowy','dyskretny','Location','best');
    grid on;

    subplot(1,2,2);
    plot(t_nl, Y_nl(:,2), 'b-', 'LineWidth', 2); hold on;
    plot(t_lin, y_lin(:,2), 'r--', 'LineWidth', 1.5);
    stairs(t_d, y_disc(:,2), 'g-.', 'LineWidth', 1.2);
    xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
    xlabel('t [min]'); ylabel('T [K]'); title('T');
    legend('nieliniowy','liniowy','dyskretny','Location','best');
    grid on;

    sgtitle(sprintf('NL vs LIN vs dyskretny (Tp=%.1f min) -- skok %s = %+g', ...
        Tp, input_names_w{ti}, amp));
    exportgraphics(gcf, sprintf('wykresy/1d_weryfikacja_%s.pdf', input_short_w{ti}), ...
        'ContentType', 'vector');
end

%% Porownanie: transmitancje dyskretne vs przestrzen stanow
G_d_tf = tf(sys_d_final);

du_test = zeros(length(t_d), 4);
du_test(t_d >= step_time, 1) = 0.1;

[y_ss, ~] = lsim(sys_d_final, du_test, t_d);
[y_tf, ~] = lsim(G_d_tf, du_test, t_d);

figure('Name','SS vs TF dyskretne','NumberTitle','off', ...
    'Position',[50 50 1000 500]);

subplot(1,2,1);
stairs(t_d, y_ss(:,1), 'b-', 'LineWidth', 2); hold on;
stairs(t_d, y_tf(:,1), 'r--', 'LineWidth', 1.5);
xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
xlabel('t [min]'); ylabel('\Delta C_A'); title('C_A');
legend('przestrzen stanow','transmitancje','Location','best');
grid on;

subplot(1,2,2);
stairs(t_d, y_ss(:,2), 'b-', 'LineWidth', 2); hold on;
stairs(t_d, y_tf(:,2), 'r--', 'LineWidth', 1.5);
xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
xlabel('t [min]'); ylabel('\Delta T'); title('T');
legend('przestrzen stanow','transmitancje','Location','best');
grid on;

sgtitle('Porownanie: przestrzen stanow vs transmitancje dyskretne (skok C_{Ain}=+0.1)');
exportgraphics(gcf, 'wykresy/1d_ss_vs_tf.pdf', 'ContentType', 'vector');

fprintf('Max roznica SS vs TF: CA=%.2e, T=%.2e\n', ...
    max(abs(y_ss(:,1)-y_tf(:,1))), max(abs(y_ss(:,2)-y_tf(:,2))));

%%

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
