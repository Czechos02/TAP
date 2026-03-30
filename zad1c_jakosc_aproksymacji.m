clear; clc; close all;
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end

% ==========================================================
% PUNKT 1c: Dyskusja jakosci przyblizenia liniowego
%   w zaleznosci od wielkosci zmian sygnalow wejsciowych
%   - Bledy bezwzgledne i wzgledne
%   - Wykresy bledu vs amplituda skoku
%   - Wnioski
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

%% LINEARYZACJA
A = jacobian_A(x0, u0, p);
B = jacobian_B(x0, u0, p);
sys_lin = ss(A, B, eye(2), zeros(2,4));

%% KONFIGURACJA
sim_end = 20;
step_time = 5;
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.01);
t_lin = (0:0.01:sim_end)';

% Amplitudy do przetestowania (gestsza siatka)
amps_CAin = [-1.0, -0.5, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.5, 1.0];
amps_FC   = [-10, -5, -2, -1, -0.5, 0.5, 1, 2, 5, 10];
amps_Tin  = [-20, -10, -5, -2, -1, 1, 2, 5, 10, 20];
amps_TCin = [-20, -10, -5, -2, -1, 1, 2, 5, 10, 20];

all_amps = {amps_CAin, amps_FC, amps_Tin, amps_TCin};
input_names = {'C_{Ain}', 'F_C', 'T_{in}', 'T_{Cin}'};
input_units = {'kmol/m^3', 'm^3/min', 'K', 'K'};

%% ============================================================
%  OBLICZENIE BLEDOW DLA KAZDEGO WEJSCIA I AMPLITUDY
% =============================================================

% Metryki bledu: max|e|, RMSE, max|e%| - dla CA i T osobno
results = struct();

for inp = 1:4
    amps = all_amps{inp};
    n = length(amps);

    err_max_CA  = zeros(1, n);
    err_max_T   = zeros(1, n);
    err_rmse_CA = zeros(1, n);
    err_rmse_T  = zeros(1, n);
    err_ss_CA   = zeros(1, n);  % blad stanu ustalonego
    err_ss_T    = zeros(1, n);

    for s = 1:n
        amp = amps(s);

        % --- Nieliniowy ---
        u_funs = cell(1,4);
        for j = 1:4
            if j == inp
                v0 = u0(j); v1 = u0(j) + amp;
                u_funs{j} = @(t) v0 + (v1 - v0) * (t >= step_time);
            else
                v = u0(j);
                u_funs{j} = @(t) v;
            end
        end

        [t_nl, Y_nl] = ode15s(@(t,y) model_nl(t, y, ...
            u_funs{1}, u_funs{2}, u_funs{3}, u_funs{4}, p), ...
            [0 sim_end], x0, opts);

        % --- Liniowy ---
        du = zeros(length(t_lin), 4);
        du(t_lin >= step_time, inp) = amp;
        [y_lin_out, ~] = lsim(sys_lin, du, t_lin);
        y_lin_out(:,1) = y_lin_out(:,1) + x0(1);
        y_lin_out(:,2) = y_lin_out(:,2) + x0(2);

        % Interpolacja nieliniowego na siatke liniowego
        Y_nl_interp = interp1(t_nl, Y_nl, t_lin, 'linear', 'extrap');

        % Bledy (tylko po skoku: t >= step_time)
        idx = t_lin >= step_time;
        e_CA = Y_nl_interp(idx, 1) - y_lin_out(idx, 1);
        e_T  = Y_nl_interp(idx, 2) - y_lin_out(idx, 2);

        err_max_CA(s)  = max(abs(e_CA));
        err_max_T(s)   = max(abs(e_T));
        err_rmse_CA(s) = sqrt(mean(e_CA.^2));
        err_rmse_T(s)  = sqrt(mean(e_T.^2));

        % Blad stanu ustalonego (ostatnie 2 min)
        idx_ss = t_lin >= (sim_end - 2);
        err_ss_CA(s) = abs(mean(Y_nl_interp(idx_ss,1) - y_lin_out(idx_ss,1)));
        err_ss_T(s)  = abs(mean(Y_nl_interp(idx_ss,2) - y_lin_out(idx_ss,2)));
    end

    results(inp).amps = amps;
    results(inp).err_max_CA  = err_max_CA;
    results(inp).err_max_T   = err_max_T;
    results(inp).err_rmse_CA = err_rmse_CA;
    results(inp).err_rmse_T  = err_rmse_T;
    results(inp).err_ss_CA   = err_ss_CA;
    results(inp).err_ss_T    = err_ss_T;
end

%% ============================================================
%  WYKRESY BLEDOW vs AMPLITUDA
% =============================================================

% --- Figura 1: Blad max vs amplituda ---
figure('Name','Blad max vs amplituda', ...
    'NumberTitle','off', 'Position',[50 50 1200 700]);

for inp = 1:4
    subplot(2, 4, inp);
    bar(results(inp).amps, results(inp).err_max_CA, 0.6, 'FaceColor', [0.2 0.4 0.8]);
    xlabel(sprintf('\\Delta %s [%s]', input_names{inp}, input_units{inp}));
    ylabel('max|e_{CA}|');
    title(sprintf('skok %s', input_names{inp}));
    grid on;

    subplot(2, 4, 4 + inp);
    bar(results(inp).amps, results(inp).err_max_T, 0.6, 'FaceColor', [0.8 0.2 0.2]);
    xlabel(sprintf('\\Delta %s [%s]', input_names{inp}, input_units{inp}));
    ylabel('max|e_T| [K]');
    title(sprintf('skok %s', input_names{inp}));
    grid on;
end
sgtitle('Blad maksymalny |e| = |y_{NL} - y_{LIN}| vs amplituda skoku');
exportgraphics(gcf, 'wykresy/1c_blad_max.pdf', 'ContentType', 'vector');

% --- Figura 2: RMSE vs amplituda ---
figure('Name','RMSE vs amplituda', ...
    'NumberTitle','off', 'Position',[100 100 1200 700]);

for inp = 1:4
    subplot(2, 4, inp);
    bar(results(inp).amps, results(inp).err_rmse_CA, 0.6, 'FaceColor', [0.2 0.6 0.4]);
    xlabel(sprintf('\\Delta %s [%s]', input_names{inp}, input_units{inp}));
    ylabel('RMSE_{CA}');
    title(sprintf('skok %s', input_names{inp}));
    grid on;

    subplot(2, 4, 4 + inp);
    bar(results(inp).amps, results(inp).err_rmse_T, 0.6, 'FaceColor', [0.6 0.4 0.2]);
    xlabel(sprintf('\\Delta %s [%s]', input_names{inp}, input_units{inp}));
    ylabel('RMSE_T [K]');
    title(sprintf('skok %s', input_names{inp}));
    grid on;
end
sgtitle('RMSE = sqrt(mean(e^2)) vs amplituda skoku');
exportgraphics(gcf, 'wykresy/1c_rmse.pdf', 'ContentType', 'vector');

% --- Figura 3: Blad stanu ustalonego ---
figure('Name','Blad stanu ustalonego vs amplituda', ...
    'NumberTitle','off', 'Position',[150 50 1200 700]);

for inp = 1:4
    subplot(2, 4, inp);
    bar(results(inp).amps, results(inp).err_ss_CA, 0.6, 'FaceColor', [0.5 0.2 0.7]);
    xlabel(sprintf('\\Delta %s [%s]', input_names{inp}, input_units{inp}));
    ylabel('|e_{ss,CA}|');
    title(sprintf('skok %s', input_names{inp}));
    grid on;

    subplot(2, 4, 4 + inp);
    bar(results(inp).amps, results(inp).err_ss_T, 0.6, 'FaceColor', [0.7 0.5 0.2]);
    xlabel(sprintf('\\Delta %s [%s]', input_names{inp}, input_units{inp}));
    ylabel('|e_{ss,T}| [K]');
    title(sprintf('skok %s', input_names{inp}));
    grid on;
end
sgtitle('Blad stanu ustalonego vs amplituda skoku');
exportgraphics(gcf, 'wykresy/1c_blad_ustalony.pdf', 'ContentType', 'vector');

% --- Figura 4: RMSE vs |amplituda| (wykres liniowy) ---
figure('Name','RMSE vs |amplituda|', ...
    'NumberTitle','off', 'Position',[200 50 1000 500]);

colors_inp = {'b', [0 0.6 0], 'r', [0.6 0 0.6]};

subplot(1,2,1); hold on;
for inp = 1:4
    pos_idx = results(inp).amps > 0;
    plot(abs(results(inp).amps(pos_idx)), results(inp).err_rmse_CA(pos_idx), ...
        '-o', 'Color', colors_inp{inp}, 'LineWidth', 1.5, 'MarkerSize', 6);
end
xlabel('|amplituda skoku|'); ylabel('RMSE_{CA}');
title('RMSE C_A vs |amplituda|');
legend(input_names, 'Location','northwest', 'FontSize', 8); grid on;

subplot(1,2,2); hold on;
for inp = 1:4
    pos_idx = results(inp).amps > 0;
    plot(abs(results(inp).amps(pos_idx)), results(inp).err_rmse_T(pos_idx), ...
        '-o', 'Color', colors_inp{inp}, 'LineWidth', 1.5, 'MarkerSize', 6);
end
xlabel('|amplituda skoku|'); ylabel('RMSE_T [K]');
title('RMSE T vs |amplituda|');
legend(input_names, 'Location','northwest', 'FontSize', 8); grid on;

sgtitle('RMSE vs |amplituda| - wzrost bledu z wielkoscia skoku');
exportgraphics(gcf, 'wykresy/1c_rmse_vs_amplituda.pdf', 'ContentType', 'vector');

% --- Figura 5: Blad procentowy vs zmiana procentowa wejscia ---
figure('Name','Blad procentowy vs zmiana % wejscia', ...
    'NumberTitle','off', 'Position',[250 50 1000 500]);

subplot(1,2,1); hold on;
for inp = 1:4
    pos_idx = results(inp).amps > 0;
    amp_pct = abs(results(inp).amps(pos_idx)) / u0(inp) * 100;
    err_pct = results(inp).err_rmse_CA(pos_idx) / x0(1) * 100;
    plot(amp_pct, err_pct, '-o', 'Color', colors_inp{inp}, 'LineWidth', 1.5, 'MarkerSize', 6);
end
xlabel('Zmiana wejscia [% wartosci nominalnej]');
ylabel('RMSE_{CA} [% wartosci nominalnej]');
title('Blad % C_A');
legend(input_names, 'Location','northwest', 'FontSize', 8); grid on;

subplot(1,2,2); hold on;
for inp = 1:4
    pos_idx = results(inp).amps > 0;
    amp_pct = abs(results(inp).amps(pos_idx)) / u0(inp) * 100;
    err_pct = results(inp).err_rmse_T(pos_idx) / x0(2) * 100;
    plot(amp_pct, err_pct, '-o', 'Color', colors_inp{inp}, 'LineWidth', 1.5, 'MarkerSize', 6);
end
xlabel('Zmiana wejscia [% wartosci nominalnej]');
ylabel('RMSE_T [% wartosci nominalnej]');
title('Blad % T');
legend(input_names, 'Location','northwest', 'FontSize', 8); grid on;

sgtitle('Blad procentowy aproksymacji liniowej vs procentowa zmiana wejscia');
exportgraphics(gcf, 'wykresy/1c_blad_procentowy.pdf', 'ContentType', 'vector');

%% ============================================================
%  PRZYKLADOWE PRZEBIEGI: male vs duze skoki (CAin i FC)
% =============================================================

% CAin: maly (0.1) vs duzy (1.0)
figure('Name','Porownanie: maly vs duzy skok CAin', ...
    'NumberTitle','off', 'Position',[50 400 1100 500]);

amps_show = [0.1, 0.5, 1.0];
colors_show = {[0 0.5 0], [0.8 0.5 0], [0.8 0 0]};

for s = 1:length(amps_show)
    amp = amps_show(s);

    u_funs = {@(t) u0(1) + amp*(t>=step_time), @(t) u0(2), @(t) u0(3), @(t) u0(4)};
    [t_nl, Y_nl] = ode15s(@(t,y) model_nl(t,y,u_funs{1},u_funs{2},u_funs{3},u_funs{4},p), ...
        [0 sim_end], x0, opts);

    du = zeros(length(t_lin), 4);
    du(t_lin >= step_time, 1) = amp;
    [y_l, ~] = lsim(sys_lin, du, t_lin);
    y_l(:,1) = y_l(:,1) + x0(1);
    y_l(:,2) = y_l(:,2) + x0(2);

    subplot(1,2,1);
    plot(t_nl, Y_nl(:,1), '-', 'Color', colors_show{s}, 'LineWidth', 2); hold on;
    plot(t_lin, y_l(:,1), '--', 'Color', colors_show{s}, 'LineWidth', 1.5);

    subplot(1,2,2);
    plot(t_nl, Y_nl(:,2), '-', 'Color', colors_show{s}, 'LineWidth', 2); hold on;
    plot(t_lin, y_l(:,2), '--', 'Color', colors_show{s}, 'LineWidth', 1.5);
end

subplot(1,2,1);
xlabel('t [min]'); ylabel('C_A [kmol/m^3]'); title('C_A');
legend('NL +0.1','LIN +0.1','NL +0.5','LIN +0.5','NL +1.0','LIN +1.0','Location','best');
grid on;

subplot(1,2,2);
xlabel('t [min]'); ylabel('T [K]'); title('T');
legend('NL +0.1','LIN +0.1','NL +0.5','LIN +0.5','NL +1.0','LIN +1.0','Location','best');
grid on;

sgtitle('Skok C_{Ain}: ciagla = NL, przerywana = LIN');
exportgraphics(gcf, 'wykresy/1c_przebiegi_CAin.pdf', 'ContentType', 'vector');

% FC: maly (1) vs duzy (10)
figure('Name','Porownanie: maly vs duzy skok FC', ...
    'NumberTitle','off', 'Position',[100 400 1100 500]);

amps_show_fc = [1, 5, 10];

for s = 1:length(amps_show_fc)
    amp = amps_show_fc(s);

    u_funs = {@(t) u0(1), @(t) u0(2)+amp*(t>=step_time), @(t) u0(3), @(t) u0(4)};
    [t_nl, Y_nl] = ode15s(@(t,y) model_nl(t,y,u_funs{1},u_funs{2},u_funs{3},u_funs{4},p), ...
        [0 sim_end], x0, opts);

    du = zeros(length(t_lin), 4);
    du(t_lin >= step_time, 2) = amp;
    [y_l, ~] = lsim(sys_lin, du, t_lin);
    y_l(:,1) = y_l(:,1) + x0(1);
    y_l(:,2) = y_l(:,2) + x0(2);

    subplot(1,2,1);
    plot(t_nl, Y_nl(:,1), '-', 'Color', colors_show{s}, 'LineWidth', 2); hold on;
    plot(t_lin, y_l(:,1), '--', 'Color', colors_show{s}, 'LineWidth', 1.5);

    subplot(1,2,2);
    plot(t_nl, Y_nl(:,2), '-', 'Color', colors_show{s}, 'LineWidth', 2); hold on;
    plot(t_lin, y_l(:,2), '--', 'Color', colors_show{s}, 'LineWidth', 1.5);
end

subplot(1,2,1);
xlabel('t [min]'); ylabel('C_A [kmol/m^3]'); title('C_A');
legend('NL +1','LIN +1','NL +5','LIN +5','NL +10','LIN +10','Location','best');
grid on;

subplot(1,2,2);
xlabel('t [min]'); ylabel('T [K]'); title('T');
legend('NL +1','LIN +1','NL +5','LIN +5','NL +10','LIN +10','Location','best');
grid on;

sgtitle('Skok F_C: ciagla = NL, przerywana = LIN');
exportgraphics(gcf, 'wykresy/1c_przebiegi_FC.pdf', 'ContentType', 'vector');

%% ============================================================
%  WYDRUK TABELKI BLEDOW
% =============================================================

fprintf('\n========================================\n');
fprintf('  PODSUMOWANIE JAKOSCI APROKSYMACJI\n');
fprintf('========================================\n\n');

for inp = 1:4
    fprintf('--- Skok %s ---\n', input_names{inp});
    fprintf('%12s | %12s | %12s | %12s | %12s\n', ...
        'Amplituda', 'maxE_CA', 'RMSE_CA', 'maxE_T', 'RMSE_T');
    fprintf('%s\n', repmat('-', 1, 70));
    for s = 1:length(results(inp).amps)
        fprintf('%+12.2f | %12.6f | %12.6f | %12.4f | %12.4f\n', ...
            results(inp).amps(s), ...
            results(inp).err_max_CA(s), ...
            results(inp).err_rmse_CA(s), ...
            results(inp).err_max_T(s), ...
            results(inp).err_rmse_T(s));
    end
    fprintf('\n');
end

fprintf('========================================\n');
fprintf('  WNIOSKI\n');
fprintf('========================================\n');
fprintf('1. Blad aproksymacji liniowej rosnie wraz ze wzrostem\n');
fprintf('   amplitudy skoku - zaleznosc jest w przyblizeniu kwadratowa\n');
fprintf('   (wynika z zaniedbanego czlonu drugiego rzedu w szeregu Taylora).\n\n');
fprintf('2. Dla malych skokow (CAin +/-0.1, FC +/-1) bledy sa pomijalne\n');
fprintf('   - model liniowy dobrze przybliza nieliniowy.\n\n');
fprintf('3. Dla duzych skokow (CAin +/-1.0, FC +/-10) bledy rosna\n');
fprintf('   znaczaco, szczegolnie w temperaturze T, ze wzgledu na\n');
fprintf('   silna nieliniowosc czlonu Arrheniusa exp(-E/(RT)).\n\n');
fprintf('4. Asymetria bledow: skoki w gore i w dol daja rozne bledy,\n');
fprintf('   co wynika z nieliniowego charakteru obiektu.\n\n');
fprintf('5. Model liniowy jest odpowiedni do regulacji wokol punktu\n');
fprintf('   pracy przy malych odchylkach (typowe dla MPC/PID).\n');

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
