clear; clc; close all;
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end

% ==========================================================
% PUNKT 1b: Porownanie modelu liniowego z nieliniowym
%   - Jedna dluga symulacja z sekwencja skokow
%   - Osobna figura dla kazdego wejscia
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

% Sekwencje skokow: [amplituda1, amplituda2, ...] w odchylkach od pp
% Kazdy skok trwa 10 min, start po 5 min stabilizacji
steps_seq = {
    [+0.1, +0.5, -0.1, -0.5],   ... % CAin
    [+1,   +5,   -1,   -5],     ... % FC
    [+5,   +10,  -5,   -10],    ... % Tin
    [+5,   +10,  -5,   -10]     ... % TCin
};

input_names = {'C_{Ain}', 'F_C', 'T_{in}', 'T_{Cin}'};
input_units = {'[kmol/m^3]', '[m^3/min]', '[K]', '[K]'};

dt = 10;         % czas trwania jednego skoku
t_start = 5;     % poczatek pierwszego skoku

%% SYMULACJA - osobna figura dla kazdego wejscia
for inp = 1:4
    amps = steps_seq{inp};
    n_steps = length(amps);
    sim_end = t_start + n_steps * dt + 5;  % +5 min na koncu
    t_lin = (0:0.01:sim_end)';

    % --- Budowanie sygnalu wejsciowego (schodkowy) ---
    du_signal = zeros(length(t_lin), 1);
    for s = 1:n_steps
        t_on  = t_start + (s-1)*dt;
        t_off = t_start + s*dt;
        du_signal(t_lin >= t_on & t_lin < t_off) = amps(s);
    end

    % Funkcja wejscia dla ode15s
    amp_vec = [0, amps, 0];
    t_breaks = [0];
    for s = 1:n_steps+1
        t_breaks(end+1) = t_start + (s-1)*dt;
    end
    t_breaks(end+1) = sim_end + 1;

    inp_fun = @(t) amp_vec(sum(t >= t_breaks));

    % --- Model nieliniowy (ode15s) ---
    u_funs = cell(1,4);
    for j = 1:4
        if j == inp
            base = u0(j);
            u_funs{j} = @(t) base + inp_fun(t);
        else
            val = u0(j);
            u_funs{j} = @(t) val;
        end
    end

    [t_nl, Y_nl] = ode15s(@(t,y) model_nl(t, y, ...
        u_funs{1}, u_funs{2}, u_funs{3}, u_funs{4}, p), ...
        [0 sim_end], x0, opts);

    % --- Model liniowy (lsim) ---
    du_lin = zeros(length(t_lin), 4);
    du_lin(:, inp) = du_signal;
    [y_lin, ~] = lsim(sys_lin, du_lin, t_lin);
    y_lin(:,1) = y_lin(:,1) + x0(1);
    y_lin(:,2) = y_lin(:,2) + x0(2);

    % Sygnal wejsciowy (wartosci bezwzgledne)
    u_abs = du_signal + u0(inp);

    % --- Wykres: 3 subploty (wejscie, CA, T) ---
    figure('Name', sprintf('1b: %s', input_names{inp}), ...
        'NumberTitle','off', 'Position',[50+60*(inp-1), 50, 900, 700]);

    % Wejscie
    subplot(3,1,1);
    stairs(t_lin, u_abs, 'k', 'LineWidth', 1.5);
    ylabel(sprintf('%s %s', input_names{inp}, input_units{inp}));
    title(sprintf('Sygnal wejsciowy %s', input_names{inp}));
    grid on;

    % CA
    subplot(3,1,2);
    plot(t_nl, Y_nl(:,1), 'b', 'LineWidth', 2); hold on;
    plot(t_lin, y_lin(:,1), 'r--', 'LineWidth', 1.5);
    ylabel('C_A [kmol/m^3]'); title('Wyjscie C_A');
    legend('nieliniowy', 'liniowy', 'Location','best'); grid on;

    % T
    subplot(3,1,3);
    plot(t_nl, Y_nl(:,2), 'b', 'LineWidth', 2); hold on;
    plot(t_lin, y_lin(:,2), 'r--', 'LineWidth', 1.5);
    xlabel('t [min]'); ylabel('T [K]'); title('Wyjscie T');
    legend('nieliniowy', 'liniowy', 'Location','best'); grid on;

    sgtitle(sprintf('Porownanie NL vs LIN - sekwencja skokow %s', input_names{inp}));
    exportgraphics(gcf, sprintf('wykresy/1b_%d.pdf', inp), 'ContentType', 'vector');
end

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
