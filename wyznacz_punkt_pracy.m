% Wyznaczanie punktu pracy reaktora CSTR

% Parametry modelu
q = 1e6; q_c = 1e6; c_p = 1; c_pc = 1;
k0 = 1e10; E_R = 8330.1; h = 130e6; a = 1.678e6; b = 0.5;
V = 1; F = 1; F_in = 1;

% Wejscia w punkcie pracy
C_Ain = 2; F_C = 15; T_in = 323; T_Cin = 365;

% Wspolczynnik wymiany ciepla
K = (a * F_C^(b+1)) / (F_C + (a * F_C^b) / (2 * q_c * c_pc));

% Rownania stanu (dy/dt = 0)
model = @(y) [
    F_in*C_Ain/V - F*y(1)/V - k0*exp(-E_R/y(2))*y(1);
    F_in*T_in/V - F*y(2)/V + (h*k0/(q*c_p))*exp(-E_R/y(2))*y(1) - K*(y(2)-T_Cin)/(V*q*c_p)
];

% Rozwiazanie ukladu rownan
options = optimoptions('fsolve', 'Display', 'iter', 'TolFun', 1e-15, 'TolX', 1e-15);
y0 = [0.26, 393.9];
[y_eq, fval, exitflag] = fsolve(model, y0, options);

% Wyniki
fprintf('\n=== Punkt pracy reaktora CSTR ===\n');
fprintf('C_A = %.10f [mol/l]\n', y_eq(1));
fprintf('T   = %.10f [K]\n', y_eq(2));

% Weryfikacja
fprintf('\n=== Weryfikacja (dy/dt w punkcie rownowagi) ===\n');
fprintf('dC_A/dt = %.2e\n', fval(1));
fprintf('dT/dt   = %.2e\n', fval(2));
fprintf('exitflag = %d\n', exitflag);
