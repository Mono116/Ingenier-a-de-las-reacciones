%%% PROYECTO DE INGENIERIA DE LAS REACCIONES
%%% Análisis de Sensibilidad y Evaluación Técnico-Económica de un Proceso Haber-Bosch Multilecho

% INTEGRANTES: Luis Mario Paredes Narvaez, Valentina Villalobos Agresott,
% Sergio Andrés Barrios Altamiranda, Juan Esteban Suarez Buendía.
% Docente: Deyler Rafael Castilla Caballero

% =========================================================================
% SCRIPT PRINCIPAL: SIMULACIÓN DE REACTORES EN SERIE PROCESO HABER-BOSCH
% =========================================================================
clear; clc; close all;

% 1. DATOS INICIALES Y PARÁMETROS   
params.R  = 8.314;          % Constante de los gases ideales [J/mol K]
params.T0 = 400 + 273.15;   % Temperatura inicial [K]
params.P0 = 150 * 101325;   % Presión inicial [Pa]
params.H  = -92400;         % Calor de reacción [J/mol]

% Parámetros cinéticos
params.alpha_TP = 0.5;      % Parámetro de Temkin-Pyzhev
params.E = 96232;           % Energía de activación [J/mol]
params.A = (9.02e8) / 3600; % Factor preexponencial [mol/kg s]

% Parámetros del Lecho
params.phi = 0.45;          % Porosidad del lecho
params.Dp  = 1500e-6;        % Diámetro de partícula [m] (1500 micras)
params.rho_b = 2800;        % Densidad bulk [kg/m3]
params.rho_iron = params.rho_b / (1 - params.phi); % Densidad del hierro [kg/m3]

% Propiedades del gas
params.mu = 2.0e-5;         % Viscosidad del gas [Pa*s]
params.G  = 26.7;           % Flujo másico superficial [kg/m2*s]

% Constantes para las capacidades caloríficas
params.a_N2  = 28.90;   params.b_N2  = -1.571e-3; params.c_N2  = 8.081e-6; params.d_N2  = -2.873e-9;
params.a_H2  = 28.84;   params.b_H2  = 7.650e-5;  params.c_H2  = 3.288e-6; params.d_H2  = -8.698e-10;
params.a_NH3 = 35.15;   params.b_NH3 = 2.954e-2;  params.c_NH3 = 4.421e-6; params.d_NH3 = -6.686e-9;

% Geometría del reactor
D_t = 3.0;                  % Diámetro del reactor industrial [m]
params.Ac = (pi * D_t^2)/4; % Área transversal [m2]

% 2. CONDICIONES INICIALES DE FLUJO
params.FA0 = 100000;            % Flujo inicial de N2 [mol/s]
params.FB0 = 300000;            % Flujo inicial de H2 [mol/s]
params.FC0 = 0.05 * params.FA0; % Flujo inicial de NH3 [mol/s] (5% de FA0)

% Vector de estado inicial: y0 = [F_N2, F_H2, F_NH3, p, T]
y0 = [params.FA0, params.FB0, params.FC0, 1, params.T0];

% 3. TERMODINÁMICA DE ENTRADA (PENG-ROBINSON) Y DENSIDAD REAL
FT0 = params.FA0 + params.FB0 + params.FC0;
y_entrada = [params.FA0 / FT0, params.FB0 / FT0, params.FC0 / FT0];

% Llamada única a la función externa calcularPR
[Z0, f_hat_entrada] = calcularPR(params.T0, params.P0, y_entrada, params);
params.Z0 = Z0;

% Peso molecular promedio en la entrada [kg/mol]
MW = [0.0280134; 0.00201588; 0.017031];
MW_prom0 = y_entrada * MW;

% Densidad real inicial [kg/m3]
params.rho_0 = (params.P0 * MW_prom0) / (params.Z0 * params.R * params.T0);

% Mostrar resultados iniciales en la ventana de comandos
disp(['Factor de compresibilidad inicial (Z0): ', num2str(params.Z0)]);
disp(['Densidad real inicial [kg/m3]: ', num2str(params.rho_0)]);
disp('Fugacidades iniciales [atm]:');
disp(['  f_N2  = ', num2str(f_hat_entrada(1)), ' atm']);
disp(['  f_H2  = ', num2str(f_hat_entrada(2)), ' atm']);
disp(['  f_NH3 = ', num2str(f_hat_entrada(3)), ' atm']);

% 4. CÁLCULO DE ALPHA PARA ERGUN (Usando params.rho_0 ya definido)
gc = 1; % Factor de conversión (1 en SI)
beta0 = (params.G * (1 - params.phi) / (params.rho_0 * gc * params.Dp * params.phi^3)) * ...
    (150 * ((1 - params.phi) * params.mu / params.Dp) + 1.75 * params.G);
params.alpha_Ergun = 2 * beta0 / (params.Ac * params.rho_iron * (1 - params.phi) * params.P0);

%% SIMULACIONES DE LECHOS %%
opciones = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% --- SIMULACIÓN DEL LECHO 1 ---
W_lecho1 = 20000; % kg de catalizador
W_span1  = [0, W_lecho1];
[W_out1, y_out1] = ode15s(@(W, y) funcion_reactor(W, y, params), W_span1, y0, opciones);

% Intercooling 1
y0_lecho2 = y_out1(end, :);
y0_lecho2(5) = 673; % Enfriamiento a 673 K

% --- SIMULACIÓN DEL LECHO 2 ---
W_lecho2 = 30000; % kg de catalizador
W_span2  = [W_out1(end), W_out1(end) + W_lecho2];
[W_out2, y_out2] = ode15s(@(W, y) funcion_reactor(W, y, params), W_span2, y0_lecho2, opciones);

% Intercooling 2
y0_lecho3 = y_out2(end, :);
y0_lecho3(5) = 673; % Enfriamiento a 673 K

% --- SIMULACIÓN DEL LECHO 3 ---
W_lecho3 = 30000; % kg de catalizador
W_span3  = [W_out2(end), W_out2(end) + W_lecho3];
[W_out3, y_out3] = ode15s(@(W, y) funcion_reactor(W, y, params), W_span3, y0_lecho3, opciones); 

% Ensamblaje de los perfiles multilecho
W_out = [W_out1; W_out2; W_out3];
y_out = [y_out1; y_out2; y_out3];

%% CÁLCULO DE CALOR EN INTERCOOLERS (Q1 y Q2)
% Función anónima para la integral analítica de Cp(T) desde T1 hasta T2
int_Cp = @(a, b, c, d, T1, T2) a*(T2 - T1) + (b/2)*(T2^2 - T1^2) + ...
    (c/3)*(T2^3 - T1^3) + (d/4)*(T2^4 - T1^4);

T_objetivo = 673; % Temperatura de entrada al siguiente lecho [K]

% Intercooler 1 (Salida del Lecho 1)
F_sal1 = y_out1(end, 1:3); % [F_N2, F_H2, F_NH3]
T_sal1 = y_out1(end, 5);

Q1 = (F_sal1(1) * int_Cp(params.a_N2,  params.b_N2,  params.c_N2,  params.d_N2,  T_objetivo, T_sal1) + ...
    F_sal1(2) * int_Cp(params.a_H2,  params.b_H2,  params.c_H2,  params.d_H2,  T_objetivo, T_sal1) + ...
    F_sal1(3) * int_Cp(params.a_NH3, params.b_NH3, params.c_NH3, params.d_NH3, T_objetivo, T_sal1)) / 1e6; % [MW]

% Intercooler 2 (Salida del Lecho 2)
F_sal2 = y_out2(end, 1:3);
T_sal2 = y_out2(end, 5);

Q2 = (F_sal2(1) * int_Cp(params.a_N2,  params.b_N2,  params.c_N2,  params.d_N2,  T_objetivo, T_sal2) + ...
    F_sal2(2) * int_Cp(params.a_H2,  params.b_H2,  params.c_H2,  params.d_H2,  T_objetivo, T_sal2) + ...
    F_sal2(3) * int_Cp(params.a_NH3, params.b_NH3, params.c_NH3, params.d_NH3, T_objetivo, T_sal2)) / 1e6; % [MW]

Q_total_MW = Q1 + Q2;
disp(['Calor Intercooler 1: ', num2str(Q1, '%.2f'), ' MW']);
disp(['Calor Intercooler 2: ', num2str(Q2, '%.2f'), ' MW']);
disp(['Calor Total Retirado: ', num2str(Q_total_MW, '%.2f'), ' MW']);

% EXTRACCIÓN Y POST-PROCESAMIENTO
F_N2  = y_out(:, 1); 
F_H2  = y_out(:, 2); 
F_NH3 = y_out(:, 3); 
p     = y_out(:, 4); 
T     = y_out(:, 5); 

P_real = p .* params.P0; 
X_N2 = (params.FA0 - F_N2) ./ params.FA0;

% 6. GRÁFICAS
figure(1);
subplot(2,2,1);
plot(W_out, F_N2, 'b-', W_out, F_H2, 'r-', W_out, F_NH3, 'g-', 'LineWidth', 1.5);
title('Perfiles Molares'); xlabel('W [kg]'); ylabel('Flujos [mol/s]');
legend('N_2', 'H_2', 'NH_3'); grid on;

subplot(2,2,2);
plot(W_out, T, 'y-', 'LineWidth', 1.5);
title('Temperatura'); xlabel('W [kg]'); ylabel('T [K]'); grid on;

subplot(2,2,3);
plot(W_out, P_real/101325, 'm-', 'LineWidth', 1.5);
title('Caída de Presión'); xlabel('W [kg]'); ylabel('P [atm]'); grid on;

subplot(2,2,4);
plot(W_out, X_N2 * 100, 'c-', 'LineWidth', 1.5);
title('Conversión N_2'); xlabel('W [kg]'); ylabel('X [%]'); grid on;

%% ANÁLISIS DE SENSIBILIDAD: EFECTO DE LA PRESIÓN DE OPERACIÓN (P0)
P_vector_atm = [120, 150, 180, 210, 240, 270]; % Escenarios de presión [atm]
n_esc = length(P_vector_atm);

% Vectores para almacenar resultados de cada escenario
X_final_vec   = zeros(1, n_esc);
T_pico_vec    = zeros(1, n_esc);
dP_total_vec  = zeros(1, n_esc);
Q_recup_vec   = zeros(1, n_esc); % [MW]
W_comp_vec    = zeros(1, n_esc); % [MW]
Prod_NH3_vec  = zeros(1, n_esc); % [ton/día]

% Constantes para compresión multietapa (4 etapas desde 1 atm) y calor
N_etapas = 4;
eta_is   = 0.78;
gamma    = 1.4;
MW       = [0.0280134; 0.00201588; 0.017031];
int_Cp   = @(a, b, c, d, T1, T2) a*(T2 - T1) + (b/2)*(T2^2 - T1^2) + ...
    (c/3)*(T2^3 - T1^3) + (d/4)*(T2^4 - T1^4);

for k = 1:n_esc
    % 1. Actualizar P0 y recalcular termodinámica/hidrodinámica de entrada
    params.P0 = P_vector_atm(k) * 101325;
    FT0 = params.FA0 + params.FB0 + params.FC0;
    y_ent = [params.FA0/FT0, params.FB0/FT0, params.FC0/FT0];

    [Z0_k, ~] = calcularPR(params.T0, params.P0, y_ent, params);
    params.Z0 = Z0_k;
    params.rho_0 = (params.P0 * (y_ent * MW)) / (params.Z0 * params.R * params.T0);

    beta0 = (params.G * (1 - params.phi) / (params.rho_0 * 1 * params.Dp * params.phi^3)) * ...
        (150 * ((1 - params.phi) * params.mu / params.Dp) + 1.75 * params.G);
    params.alpha_Ergun = 2 * beta0 / (params.Ac * params.rho_iron * (1 - params.phi) * params.P0);

    % 2. Resolver los 3 lechos en serie
    y0_k = [params.FA0, params.FB0, params.FC0, 1, params.T0];
    [~, y1] = ode15s(@(W, y) funcion_reactor(W, y, params), [0, 20000], y0_k, opciones);

    y0_2 = y1(end, :); y0_2(5) = params.T0;
    [~, y2] = ode15s(@(W, y) funcion_reactor(W, y, params), [20000, 50000], y0_2, opciones);

    y0_3 = y2(end, :); y0_3(5) = params.T0;
    [~, y3] = ode15s(@(W, y) funcion_reactor(W, y, params), [50000, 80000], y0_3, opciones);

    % 3. Extraer KPIs técnicos
    X_final_vec(k)  = ((params.FA0 - y3(end, 1)) / params.FA0) * 100;
    T_pico_vec(k)   = max([y1(:,5); y2(:,5); y3(:,5)]);
    dP_total_vec(k) = P_vector_atm(k) * (1 - y3(end, 4));

    % Producción neta de NH3 [ton/día]
    F_NH3_neto = y3(end, 3) - params.FC0; % [mol/s]
    Prod_NH3_vec(k) = F_NH3_neto * 0.017031 * (86400 / 1000);

    % 4. Calor recuperable (Intercoolers 1 y 2 hasta T0, y caldera salida hasta 523 K)
    Q1 = (y1(end,1)*int_Cp(params.a_N2,params.b_N2,params.c_N2,params.d_N2,params.T0,y1(end,5)) + ...
        y1(end,2)*int_Cp(params.a_H2,params.b_H2,params.c_H2,params.d_H2,params.T0,y1(end,5)) + ...
        y1(end,3)*int_Cp(params.a_NH3,params.b_NH3,params.c_NH3,params.d_NH3,params.T0,y1(end,5))) / 1e6;

    Q2 = (y2(end,1)*int_Cp(params.a_N2,params.b_N2,params.c_N2,params.d_N2,params.T0,y2(end,5)) + ...
        y2(end,2)*int_Cp(params.a_H2,params.b_H2,params.c_H2,params.d_H2,params.T0,y2(end,5)) + ...
        y2(end,3)*int_Cp(params.a_NH3,params.b_NH3,params.c_NH3,params.d_NH3,params.T0,y2(end,5))) / 1e6;

    Q3 = (y3(end,1)*int_Cp(params.a_N2,params.b_N2,params.c_N2,params.d_N2,523.15,y3(end,5)) + ...
        y3(end,2)*int_Cp(params.a_H2,params.b_H2,params.c_H2,params.d_H2,523.15,y3(end,5)) + ...
        y3(end,3)*int_Cp(params.a_NH3,params.b_NH3,params.c_NH3,params.d_NH3,523.15,y3(end,5))) / 1e6;

    Q_recup_vec(k) = Q1 + Q2 + Q3;

    % 5. Potencia de compresión desde 1 atm hasta P_vector_atm(k) [MW]
    rp = (P_vector_atm(k) / 1)^(1 / N_etapas);
    Z_prom = (1 + Z0_k) / 2;
    W_comp_vec(k) = (N_etapas / eta_is) * FT0 * (gamma / (gamma - 1)) * ...
        Z_prom * params.R * 313.15 * (rp^((gamma - 1)/gamma) - 1) / 1e6;
end

%% GRÁFICAS DE SENSIBILIDAD TÉCNICO-ENERGÉTICA
figure(2);
subplot(2,2,1);
plot(P_vector_atm, X_final_vec, 'o-b', 'LineWidth', 1.8);
title('Conversión Global vs. Presión'); xlabel('P_0 [atm]'); ylabel('X_{N2} [%]'); grid on;

subplot(2,2,2);
plot(P_vector_atm, T_pico_vec, 's-r', 'LineWidth', 1.8);
yline(780, '--k', 'Límite Sinterización (780 K)');
title('Temperatura Pico vs. Presión'); xlabel('P_0 [atm]'); ylabel('T_{max} [K]'); grid on;

subplot(2,2,3);
plot(P_vector_atm, W_comp_vec, '^-m', P_vector_atm, Q_recup_vec, 'v-g', 'LineWidth', 1.8);
title('Balance Energético (Gasto vs. Recuperación)'); xlabel('P_0 [atm]'); ylabel('Potencia [MW]');
legend('Compresión (Gasto)', 'Calor Recuperado (Crédito)', 'Location', 'northwest'); grid on;

subplot(2,2,4);
plot(P_vector_atm, Prod_NH3_vec, 'd-c', 'LineWidth', 1.8);
title('Producción Neta de Amoníaco'); xlabel('P_0 [atm]'); ylabel('NH_3 [ton/día]'); grid on;

%% ANÁLISIS DE SENSIBILIDAD: EFECTO DE LA TEMPERATURA DE ENTRADA (T0)
params.P0 = 150 * 101325; % Restaurar presión base [Pa]
T_vector_K = [623.15, 643.15, 663.15, 673.15, 693.15, 713.15, 733.15]; % [350 a 460 °C]
n_esc_T = length(T_vector_K);

% Vectores para almacenar resultados
X_final_T  = zeros(1, n_esc_T);
T_pico_T   = zeros(1, n_esc_T);
Q_recup_T  = zeros(1, n_esc_T); % [MW]
W_comp_T   = zeros(1, n_esc_T); % [MW]
Prod_NH3_T = zeros(1, n_esc_T); % [ton/día]

for k = 1:n_esc_T
    % 1. Actualizar T0 y recalcular propiedades termodinámicas/hidrodinámicas
    params.T0 = T_vector_K(k);
    FT0 = params.FA0 + params.FB0 + params.FC0;
    y_ent = [params.FA0/FT0, params.FB0/FT0, params.FC0/FT0];

    [Z0_k, ~] = calcularPR(params.T0, params.P0, y_ent, params);
    params.Z0 = Z0_k;
    params.rho_0 = (params.P0 * (y_ent * MW)) / (params.Z0 * params.R * params.T0);

    beta0 = (params.G * (1 - params.phi) / (params.rho_0 * 1 * params.Dp * params.phi^3)) * ...
        (150 * ((1 - params.phi) * params.mu / params.Dp) + 1.75 * params.G);
    params.alpha_Ergun = 2 * beta0 / (params.Ac * params.rho_iron * (1 - params.phi) * params.P0);

    % 2. Resolver los 3 lechos en serie (con intercooling hasta params.T0)
    y0_k = [params.FA0, params.FB0, params.FC0, 1, params.T0];
    [~, y1] = ode15s(@(W, y) funcion_reactor(W, y, params), [0, 20000], y0_k, opciones);

    y0_2 = y1(end, :); y0_2(5) = params.T0;
    [~, y2] = ode15s(@(W, y) funcion_reactor(W, y, params), [20000, 50000], y0_2, opciones);

    y0_3 = y2(end, :); y0_3(5) = params.T0;
    [~, y3] = ode15s(@(W, y) funcion_reactor(W, y, params), [50000, 80000], y0_3, opciones);

    % 3. Extraer KPIs técnicos
    X_final_T(k)  = ((params.FA0 - y3(end, 1)) / params.FA0) * 100;
    T_pico_T(k)   = max([y1(:,5); y2(:,5); y3(:,5)]);

    F_NH3_neto    = y3(end, 3) - params.FC0;
    Prod_NH3_T(k) = F_NH3_neto * 0.017031 * (86400 / 1000);

    % 4. Calor recuperable (Intercoolers 1 y 2 hasta params.T0, caldera hasta 523.15 K)
    Q1 = (y1(end,1)*int_Cp(params.a_N2,params.b_N2,params.c_N2,params.d_N2,params.T0,y1(end,5)) + ...
        y1(end,2)*int_Cp(params.a_H2,params.b_H2,params.c_H2,params.d_H2,params.T0,y1(end,5)) + ...
        y1(end,3)*int_Cp(params.a_NH3,params.b_NH3,params.c_NH3,params.d_NH3,params.T0,y1(end,5))) / 1e6;

    Q2 = (y2(end,1)*int_Cp(params.a_N2,params.b_N2,params.c_N2,params.d_N2,params.T0,y2(end,5)) + ...
        y2(end,2)*int_Cp(params.a_H2,params.b_H2,params.c_H2,params.d_H2,params.T0,y2(end,5)) + ...
        y2(end,3)*int_Cp(params.a_NH3,params.b_NH3,params.c_NH3,params.d_NH3,params.T0,y2(end,5))) / 1e6;

    Q3 = (y3(end,1)*int_Cp(params.a_N2,params.b_N2,params.c_N2,params.d_N2,523.15,y3(end,5)) + ...
        y3(end,2)*int_Cp(params.a_H2,params.b_H2,params.c_H2,params.d_H2,523.15,y3(end,5)) + ...
        y3(end,3)*int_Cp(params.a_NH3,params.b_NH3,params.c_NH3,params.d_NH3,523.15,y3(end,5))) / 1e6;

    Q_recup_T(k) = Q1 + Q2 + Q3;

    % 5. Potencia de compresión (constante a P0 = 150 atm)
    rp = (150 / 1)^(1 / N_etapas);
    Z_prom = (1 + Z0_k) / 2;
    W_comp_T(k) = (N_etapas / eta_is) * FT0 * (gamma / (gamma - 1)) * ...
        Z_prom * params.R * 313.15 * (rp^((gamma - 1)/gamma) - 1) / 1e6;
end

%% GRÁFICAS DE SENSIBILIDAD POR TEMPERATURA (FIGURE 3)
figure(3);
subplot(2,2,1);
plot(T_vector_K, X_final_T, 'o-b', 'LineWidth', 1.8);
title('Conversión Global vs. T_0'); xlabel('T_0 [K]'); ylabel('X_{N2} [%]'); grid on;

subplot(2,2,2);
plot(T_vector_K, T_pico_T, 's-r', 'LineWidth', 1.8);
yline(780, '--k', 'Límite Sinterización (780 K)');
title('Temperatura Pico vs. T_0'); xlabel('T_0 [K]'); ylabel('T_{max} [K]'); grid on;

subplot(2,2,3);
plot(T_vector_K, W_comp_T, '^-m', T_vector_K, Q_recup_T, 'v-g', 'LineWidth', 1.8);
title('Balance Energético vs. T_0'); xlabel('T_0 [K]'); ylabel('Potencia [MW]');
legend('Compresión (Gasto)', 'Calor Recuperado (Crédito)', 'Location', 'west'); grid on;

subplot(2,2,4);
plot(T_vector_K, Prod_NH3_T, 'd-c', 'LineWidth', 1.8);
title('Producción Neta de Amoníaco vs. T_0'); xlabel('T_0 [K]'); ylabel('NH_3 [ton/día]'); grid on;
%% MODELADO %%

function dydW = funcion_reactor(~, y, params)
dydW = zeros(5, 1);

% 1. Desempaquetar el vector de estado 'y'
F_N2  = y(1);
F_H2  = y(2);
F_NH3 = y(3);
p     = y(4); % Presión reducida (P/P0)
T     = y(5); % Temperatura actual [K]

% 2. Cálculos hidrodinámicos y termodinámicos locales (Peng-Robinson)
F_T  = F_N2 + F_H2 + F_NH3;
F_T0 = params.FA0 + params.FB0 + params.FC0;

P_local = p * params.P0; % [Pa]
y_local = [F_N2 / F_T, F_H2 / F_T, F_NH3 / F_T];

% Evaluar Z y fugacidades parciales [atm] en cada paso dW
[Z, f_hat] = calcularPR(T, P_local, y_local, params);

f_N2  = f_hat(1);
f_H2  = f_hat(2);
f_NH3 = f_hat(3);

if f_NH3 < 1e-10
    f_NH3 = 1e-10; 
end

% 3. Constante cinética inversa (Arrhenius) y de equilibrio (Ka)
k_reversa = params.A * exp(-params.E / (params.R * T));

log10_Ka = -2.691122 * log10(T) ...
    - 5.51965e-5 * T ...
    + 1.848863e-6 * T^2 ...
    + (2001.3 / T) ...
    + 2.689;
Ka = 10^(log10_Ka);

% 4. Ley de velocidad de Temkin-Pyzhev usando FUGACIDADES
r_NH3 = k_reversa * ( Ka^2 * f_N2 * (f_H2^3 / f_NH3^2)^params.alpha_TP - ...
    (f_NH3^2 / f_H2^3)^(1 - params.alpha_TP) );

r_N2 = -0.5 * r_NH3;
r_H2 = -1.5 * r_NH3;

% 5. Balances molares
dydW(1) = r_N2;
dydW(2) = r_H2;
dydW(3) = r_NH3;

% 6. Balance de momentum (Ecuación de Ergun corregida con Z / Z0)
dydW(4) = (-params.alpha_Ergun / (2 * p)) * (F_T / F_T0) * (T / params.T0) * (Z / params.Z0);

% 7. Balance de energía
Cp_N2  = params.a_N2  + params.b_N2*T  + params.c_N2*T^2  + params.d_N2*T^3;
Cp_H2  = params.a_H2  + params.b_H2*T  + params.c_H2*T^2  + params.d_H2*T^3;
Cp_NH3 = params.a_NH3 + params.b_NH3*T + params.c_NH3*T^2 + params.d_NH3*T^3;

Delta_Cp  = 2*Cp_NH3 - Cp_N2 - 3*Cp_H2;
DeltaH_Rx = params.H + Delta_Cp * (T - 298.15);

sum_Fi_Cpi = F_N2*Cp_N2 + F_H2*Cp_H2 + F_NH3*Cp_NH3;

dydW(5) = (r_N2 * DeltaH_Rx) / sum_Fi_Cpi;
end