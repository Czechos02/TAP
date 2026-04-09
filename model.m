clear; clc; close all;
if ~exist('wykresy', 'dir'), mkdir('wykresy'); end

%% Punkt pracy
[y0, u0, p] = punkt_pracy();
sim_time = [0, 20];

%% Symulacja
step_time = 5;
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.01);

scenarios = {
    [0.1,  0,   0,  0],  'skok C_{Ain} = +0.1';
    [0,    1,   0,  0],  'skok F_C = +1';
    [0,    0,   5,  0],  'skok T_{in} = +5';
    [0,    0,   0,  5],  'skok T_{Cin} = +5';
};

colors = {'b', [0 0.6 0], 'r', [0.6 0 0.6]};

% Punkt pracy
u1 = @(t) u0(1); u2 = @(t) u0(2);
z1 = @(t) u0(3); z2 = @(t) u0(4);
[t_pp, Y_pp] = ode15s(@(t,y) model_nl(t,y,u1,u2,z1,z2,p), sim_time, y0, opts);

figure('Name','Punkt pracy','NumberTitle','off','Position',[50 400 900 400]);
subplot(1,2,1);
plot(t_pp, Y_pp(:,1), 'b', 'LineWidth', 1.5);
xlabel('t [min]'); ylabel('C_A [kmol/m^3]'); title('C_A w punkcie pracy');
grid on;
subplot(1,2,2);
plot(t_pp, Y_pp(:,2), 'r', 'LineWidth', 1.5);
xlabel('t [min]'); ylabel('T [K]'); title('T w punkcie pracy');
grid on;
sgtitle('Symulacja modelu nieliniowego - punkt pracy');
exportgraphics(gcf, 'wykresy/model_punkt_pracy.pdf', 'ContentType', 'vector');

% Odpowiedzi skokowe
figure('Name','Odpowiedzi skokowe','NumberTitle','off','Position',[50 50 1000 450]);

for sc = 1:size(scenarios, 1)
    du = scenarios{sc, 1};

    u1 = @(t) u0(1) + du(1)*(t >= step_time);
    u2 = @(t) u0(2) + du(2)*(t >= step_time);
    z1 = @(t) u0(3) + du(3)*(t >= step_time);
    z2 = @(t) u0(4) + du(4)*(t >= step_time);

    [t_sc, Y_sc] = ode15s(@(t,y) model_nl(t,y,u1,u2,z1,z2,p), sim_time, y0, opts);

    subplot(1,2,1); hold on;
    plot(t_sc, Y_sc(:,1), 'Color', colors{sc}, 'LineWidth', 1.5);

    subplot(1,2,2); hold on;
    plot(t_sc, Y_sc(:,2), 'Color', colors{sc}, 'LineWidth', 1.5);
end

leg_names = scenarios(:,2)';

subplot(1,2,1);
xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
xlabel('t [min]'); ylabel('C_A [kmol/m^3]'); title('Wyjscie C_A');
legend(leg_names, 'Location','best', 'FontSize', 8); grid on;

subplot(1,2,2);
xline(step_time, ':', 'Color', [0.5 0.5 0.5]);
xlabel('t [min]'); ylabel('T [K]'); title('Wyjscie T');
legend(leg_names, 'Location','best', 'FontSize', 8); grid on;

sgtitle('Model nieliniowy - odpowiedzi na skoki poszczegolnych wejsc (t_{skok} = 5 min)');
exportgraphics(gcf, 'wykresy/model_odpowiedzi_skokowe.pdf', 'ContentType', 'vector');

%%

function dy_dt = model_nl(t, y, u1, u2, z1, z2, p)

uu1 = u1(t);
uu2 = u2(t);
zz1 = z1(t);
zz2 = z2(t);

dy_dt = [
    p.F_in * uu1 / p.V - p.F * y(1) / p.V - p.k0*exp(-p.E_R/y(2)) * y(1);
    zz1*p.F_in/p.V - y(2)*p.F/p.V + (p.h*p.k0/(p.q*p.c_p))*exp(-p.E_R/y(2)) * y(1) - ( y(2) - zz2)*((p.a*uu2^(p.b+1))/(uu2+(p.a*uu2^p.b)/(2*p.q*p.c_pc)))/(p.V*p.q*p.c_p);
];

end


function out = step_function(t, before, after, step_time)
    if t < step_time
        out = before;
    else
        out = after;
    end
end
