clear; clc; close all;
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end

%% Punkt pracy
[x0, u0, p] = punkt_pracy();

%% Linearyzacja
A = jacobian_A(x0, u0, p);
B = jacobian_B(x0, u0, p);
sys_lin = ss(A, B, eye(2), zeros(2,4));

%% Konfiguracja
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.01);
sim_end = 20;
step_time = 5;
t_lin = (0:0.01:sim_end)';

amps_all = {
    [-1, -0.5, -0.1, +0.1, +0.5, +1],      ... % CAin
    [-10, -5, -1, +1, +5, +10],             ... % FC
    [-20, -10, -2, +2, +10, +20],           ... % Tin
    [-20, -10, -2, +2, +10, +20]            ... % TCin
};

input_names = {'C_{Ain}', 'F_C', 'T_{in}', 'T_{Cin}'};
input_units = {'kmol/m^3', 'm^3/min', 'K', 'K'};
input_short = {'CAin', 'FC', 'Tin', 'TCin'};

colors = [
    0.0  0.0  0.8;
    0.2  0.5  1.0;
    0.0  0.7  0.3;
    0.9  0.7  0.0;
    1.0  0.4  0.0;
    0.8  0.0  0.0;
];

%% Fan-ploty
for inp = 1:4
    amps = amps_all{inp};

    figure('Name', sprintf('1b - %s', input_short{inp}), ...
        'NumberTitle','off', 'Position',[50 50 1000 500]);

    leg_entries = {};

    for ai = 1:length(amps)
        amp = amps(ai);
        col = colors(ai, :);

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

        du = zeros(length(t_lin), 4);
        du(t_lin >= step_time, inp) = amp;
        [y_lin, ~] = lsim(sys_lin, du, t_lin);
        y_lin(:,1) = y_lin(:,1) + x0(1);
        y_lin(:,2) = y_lin(:,2) + x0(2);

        subplot(1,2,1); hold on;
        plot(t_nl, Y_nl(:,1), '-', 'Color', col, 'LineWidth', 1.5);
        plot(t_lin, y_lin(:,1), '--', 'Color', col, 'LineWidth', 1.2);

        subplot(1,2,2); hold on;
        plot(t_nl, Y_nl(:,2), '-', 'Color', col, 'LineWidth', 1.5);
        plot(t_lin, y_lin(:,2), '--', 'Color', col, 'LineWidth', 1.2);

        leg_entries{end+1} = sprintf('NL, \\Delta=%+g', amp);
        leg_entries{end+1} = sprintf('LIN, \\Delta=%+g', amp);
    end

    subplot(1,2,1);
    xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
    xlabel('t [min]'); ylabel('C_A [kmol/m^3]'); title('C_A');
    legend(leg_entries, 'Location','best', 'FontSize', 6);
    grid on;

    subplot(1,2,2);
    xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
    xlabel('t [min]'); ylabel('T [K]'); title('T');
    legend(leg_entries, 'Location','best', 'FontSize', 6);
    grid on;

    sgtitle(sprintf('NL vs LIN -- skoki %s [%s]', input_names{inp}, input_units{inp}));
    exportgraphics(gcf, sprintf('wykresy/1b_%s.pdf', input_short{inp}), 'ContentType', 'vector');
end

fprintf('Wygenerowano 4 wykresy.\n');

%%

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
