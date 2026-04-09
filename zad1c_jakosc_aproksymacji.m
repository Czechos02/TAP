clear; clc; close all;
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end

%% Punkt pracy
[x0, u0, p] = punkt_pracy();

%% Linearyzacja
A = jacobian_A(x0, u0, p);
B = jacobian_B(x0, u0, p);
sys_lin = ss(A, B, eye(2), zeros(2,4));

%% Konfiguracja
sim_end = 20;
step_time = 5;
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.01);
t_lin = (0:0.01:sim_end)';

amps_all = {
    [-0.5, -0.2, -0.1, -0.05, +0.05, +0.1, +0.2, +0.5],  ... % CAin
    [-5, -2, -1, -0.5, +0.5, +1, +2, +5],                  ... % FC
    [-10, -5, -2, -1, +1, +2, +5, +10],                     ... % Tin
    [-10, -5, -2, -1, +1, +2, +5, +10]                      ... % TCin
};

input_names     = {'C_{Ain}', 'F_C', 'T_{in}', 'T_{Cin}'};
input_names_tex = {'$C_{Ain}$', '$F_C$', '$T_{in}$', '$T_{Cin}$'};
input_units     = {'kmol/m^3', 'm^3/min', 'K', 'K'};

%% Obliczenie RMSE
results = struct();

for inp = 1:4
    amps = amps_all{inp};
    n = length(amps);

    rmse_CA = zeros(1, n);
    rmse_T  = zeros(1, n);

    for s = 1:n
        amp = amps(s);

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

        du = zeros(length(t_lin), 4);
        du(t_lin >= step_time, inp) = amp;
        [y_lin_out, ~] = lsim(sys_lin, du, t_lin);
        y_lin_out(:,1) = y_lin_out(:,1) + x0(1);
        y_lin_out(:,2) = y_lin_out(:,2) + x0(2);

        Y_nl_interp = interp1(t_nl, Y_nl, t_lin, 'linear', 'extrap');

        idx = t_lin >= step_time;
        e_CA = Y_nl_interp(idx, 1) - y_lin_out(idx, 1);
        e_T  = Y_nl_interp(idx, 2) - y_lin_out(idx, 2);

        rmse_CA(s) = sqrt(mean(e_CA.^2));
        rmse_T(s)  = sqrt(mean(e_T.^2));
    end

    results(inp).amps    = amps;
    results(inp).rmse_CA = rmse_CA;
    results(inp).rmse_T  = rmse_T;
end

%% Zapis do CSV
fid_csv = fopen('wyniki_1c.csv', 'w');
fprintf(fid_csv, 'input,amplitude,rmse_CA,rmse_T\n');
for inp = 1:4
    for s = 1:length(results(inp).amps)
        fprintf(fid_csv, '%d,%.4f,%.8f,%.8f\n', ...
            inp, results(inp).amps(s), ...
            results(inp).rmse_CA(s), results(inp).rmse_T(s));
    end
end
fclose(fid_csv);
fprintf('Wyniki zapisane do: wyniki_1c.csv\n');

%% Wydruk tabelek w konsoli
fprintf('\n');
for inp = 1:4
    fprintf('=== Skok %s [%s] ===\n', input_names{inp}, input_units{inp});
    fprintf('%12s | %14s | %14s\n', 'Amplituda', 'RMSE CA', 'RMSE T [K]');
    fprintf('%s\n', repmat('-', 1, 46));
    for s = 1:length(results(inp).amps)
        fprintf('%+12.2f | %14.6f | %14.6f\n', ...
            results(inp).amps(s), ...
            results(inp).rmse_CA(s), ...
            results(inp).rmse_T(s));
    end
    fprintf('\n');
end

%% Generacja tabelek LaTeX
fid = fopen('tabele_1c.tex', 'w');

for inp = 1:4
    fprintf(fid, '\\begin{table}[H]\n\\centering\n');
    fprintf(fid, '\\caption{RMSE aproksymacji liniowej -- skok %s}\n', input_names_tex{inp});
    fprintf(fid, '\\label{tab:1c_%d}\n', inp);
    fprintf(fid, '\\begin{tabular}{r|c|c}\n');
    fprintf(fid, '\\hline\n');
    fprintf(fid, '$\\Delta %s$ [%s] & RMSE $C_A$ & RMSE $T$ [K] \\\\\n', ...
        input_names{inp}, input_units{inp});
    fprintf(fid, '\\hline\n');

    for s = 1:length(results(inp).amps)
        amp = results(inp).amps(s);
        r_ca = results(inp).rmse_CA(s);
        r_t  = results(inp).rmse_T(s);

        if r_ca < 0.001
            str_ca = sprintf('%.2e', r_ca);
        else
            str_ca = sprintf('%.4f', r_ca);
        end
        if r_t < 0.01
            str_t = sprintf('%.2e', r_t);
        else
            str_t = sprintf('%.4f', r_t);
        end

        fprintf(fid, '%+.2f & %s & %s \\\\\n', amp, str_ca, str_t);
    end

    fprintf(fid, '\\hline\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n\n');
end

fclose(fid);
fprintf('Tabelki LaTeX zapisane do: tabele_1c.tex\n');

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
