
clear; clc; close all;
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end

% ==========================================================
% PUNKT 1b: Porownanie modelu liniowego z nieliniowym
%   - Osobny wykres dla kazdego wejscia i amplitudy
%   - Na kazdym: NL vs LIN, 2 subploty (CA, T)
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
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.01);
sim_end = 20;
step_time = 5;
t_lin = (0:0.01:sim_end)';

% Amplitudy dla kazdego wejscia
amps_all = {
    [-0.1, -0.5, +0.1, +0.5],      ... % CAin
    [-1,   -5,   +1,   +5],        ... % FC
    [-2,   -10,  +2,   +10],       ... % Tin
    [-2,   -10,  +2,   +10]        ... % TCin
};

input_names = {'C_{Ain}', 'F_C', 'T_{in}', 'T_{Cin}'};
input_units = {'kmol/m^3', 'm^3/min', 'K', 'K'};
input_short = {'CAin', 'FC', 'Tin', 'TCin'};

%% SYMULACJA
for inp = 1:4
    amps = amps_all{inp};

    for ai = 1:length(amps)
        amp = amps(ai);

        % --- Nieliniowy (ode15s) ---
        u_funs = cell(1,4);
        for j = 1:4
            if j == inp
                v0 = u0(j);
                u_funs{j} = @(t) v0 + amp * (t >= step_time);
            else
                v = u0(j);
                u_funs{j} = @(t) v;
            end
        end

        [t_nl, Y_nl] = ode15s(@(t,y) model_nl(t, y, ...
            u_funs{1}, u_funs{2}, u_funs{3}, u_funs{4}, p), ...
            [0 sim_end], x0, opts);

        % --- Liniowy (lsim) ---
        du = zeros(length(t_lin), 4);
        du(t_lin >= step_time, inp) = amp;
        [y_lin, ~] = lsim(sys_lin, du, t_lin);
        y_lin(:,1) = y_lin(:,1) + x0(1);
        y_lin(:,2) = y_lin(:,2) + x0(2);

        % --- Wykres ---
        figure('Name', sprintf('%s = %+g', input_short{inp}, amp), ...
            'NumberTitle','off', 'Position',[50 50 900 400]);

        subplot(1,2,1);
        plot(t_nl, Y_nl(:,1), 'b-', 'LineWidth', 2); hold on;
        plot(t_lin, y_lin(:,1), 'r--', 'LineWidth', 1.5);
        xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
        xlabel('t [min]'); ylabel('C_A [kmol/m^3]');
        title('C_A'); legend('nieliniowy','liniowy','Location','best');
        grid on;

        subplot(1,2,2);
        plot(t_nl, Y_nl(:,2), 'b-', 'LineWidth', 2); hold on;
        plot(t_lin, y_lin(:,2), 'r--', 'LineWidth', 1.5);
        xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
        xlabel('t [min]'); ylabel('T [K]');
        title('T'); legend('nieliniowy','liniowy','Location','best');
        grid on;

        sgtitle(sprintf('Skok %s = %+g %s', ...
            input_names{inp}, amp, input_units{inp}));

        % Nazwa pliku: 1b_CAin_p0.5.pdf, 1b_FC_m2.pdf, etc.
        if amp >= 0
            amp_str = sprintf('p%.4g', amp);
        else
            amp_str = sprintf('m%.4g', abs(amp));
        end
        fname = sprintf('wykresy/1b_%s_%s.pdf', input_short{inp}, amp_str);
        exportgraphics(gcf, fname, 'ContentType', 'vector');
        close(gcf);
    end
end

fprintf('Wygenerowano %d wykresow.\n', ...
    sum(cellfun(@length, amps_all)));

%% ============================================================
%  FUNKCJE LOKALNE
% =============================================================

function dy = model_nl(t, y, u1, u2, z1, z2, p)
    C_Ain = u1(t); F_C = u2(t); T_in = z1(t); T_Cin = z2(t);
    C_A = y(1); T = y(2);
    exp_term = exp(-p.E_R / T);
    K = (p.a * F_C^(p.b+1)) / (F_C + (p.a * F_C^p.b) / (2*p.q_c*p.c_pc));
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
    N  = p.a * F_C^(p.b+1);
    D  = F_C + (p.a*F_C^p.b)/(2*p.q_c*p.c_pc);
    dN = p.a * (p.b+1) * F_C^p.b;
    dD = 1 + (p.a*p.b*F_C^(p.b-1))/(2*p.q_c*p.c_pc);
    dK = (dN*D - N*dD) / D^2;
    B = zeros(2,4);
    B(1,1) = p.F_in / p.V;
    B(2,2) = -(T - T_Cin) * dK / (p.V*p.q*p.c_p);
    B(2,3) = p.F_in / p.V;
    B(2,4) = K / (p.V*p.q*p.c_p);
end
