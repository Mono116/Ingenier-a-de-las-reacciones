%% PROYECTO DE INGENIERIA DE LAS REACCIONES
%% Análisis de Sensibilidad y Evaluación Técnico-Económica de un Proceso Haber-Bosch Multilecho

% INTEGRANTES: Luis Mario Paredes Narvaez, Valentina Villalobos Agresott,
% Sergio Andrés Barrios Altamiranda, Juan Esteban Suarez Buendía.
% Docente: Deyler Rafael Castilla Caballero


% Algoritmo inicial:
% Previamente se tiene en cuenta:
% Definir los parametros iniciales a partír de la literatura:
% alpha_TP, alpha_Ergun, beta, P0, T0, a, b, c, d, calor de reacción
% Definir flujos molares: FA0, FB0, FC0, FD0
% Definir las constantes de Arrhenius: (A, E), Definir el Wtotal
%  Balances molares
% A. Ecuaciones de diseño para cada especie: dFj/dW = r'j

% B. Usar la ley de velocidad
% de Temkin-Pyzhep

% 1. Para cada paso de integración dW MATLAB proporciona los valores
% numéricos actuales de FN2, FH2, FNH3, T y P
% 2. Cálculo del flujo Total
% FT = FN2 + FH2 + FNH3 + FInertes

% 3. Cálculo de las presiones parciales
% Recordar la definición de fracción molar yi = Fj/FT y
% de presión parcial Pj = yi*P

% 4. Evaluación de
% Temkin-Pyzhep empleando las presiones parciales calculadas
% k1 y k2 dependen de la temperatura por Arrhenius

% 5. Aplicación de las leyes de velocidad relativas
% Es decir, se obtendrá la ley de velocidad para cada especie

% 6. Retorno de las derivadas en MATLAB
% En este paso se agrupan los balances molares y las leyes de velocidad
% para cada especie, respectivamente. También se hace uso de la ecuación
% de caida de presión

% Aplicación del modelo no ideal:
% a) Se usa COCO para generar una matriz de datos calculando Z y \Phi_i
% para la mezcla a diferentes temperaturas (400 y 450 °C  y 150, 200 y 250
% bares según la segunda entrega del proyecto) (se evaluará si los rangos 
% deben ser o no cambiados)

% b) Exportar la tabla y cargar en MATLAB

% c) Crear correlaciones polinomiales (se usa fit o polyfitn) (Consultar a
% fondo con IA)

% d) Llamar a las funciones polinomiales dentro de las EDOS y resolver con
% el solver.

% NOTA: La primera parte del script asume el comportamiento ideal Z = 1. Despues
% se hacen las correciones con el modelo de Peng-Robinson.


% =========================================================================
% SCRIPT PRINCIPAL: SIMULACIÓN REACTOR HABER-BOSCH
% =========================================================================

% 1. DATOS INICIALES Y PARÁMETROS 
params.R = 8.314;           % Constante de los gases ideales [J/mol K]
params.T0 = 400 + 273.15;   % Temperatura inicial [K]
params.P0 = 150 * 101325;   % Presión inicial [Pa]
params.H = -92400;          % Calor de reacción [J/mol]

% Parámetros cinéticos
params.alpha_TP = 0.5;      % Parámetro de Temkin-Pyzhev
params.E = 96232;           % Energía de activación [J/mol]
params.A = (9.02e8) / 3600; % Factor preexponencial [mol/kg s]

% Parámetros del Lecho y Ecuación de Ergun
params.phi = 0.45;          % Porosidad del lecho
params.Dp = 300e-6;         % Diámetro de partícula [m] (300 micras)
params.rho_b = 2500;        % Densidad bulk [kg/m3]
params.rho_iron = params.rho_b / (1 - params.phi); % Densidad del hierro [kg/m3]

% Propiedades del gas a la entrada
params.mu = 2.0e-5;         % Viscosidad del gas [Pa*s]
params.rho_0 = 26.7;        % Densidad del gas en la entrada [kg/m3]
params.G = 26.7;            % Flujo másico superficial [kg/m2*s]
% Constantes para las capacidades calorificas
params.a_N2 = 28.90;
params.b_N2 = -1.571e-3;
params.c_N2 = 8.081e-6;
params.d_N2 = -2.873e-9;

params.a_H2 = 28.84;
params.b_H2 = 7.650e-5;
params.c_H2 = 3.288e-6;
params.d_H2 = -8.698e-10;

params.a_NH3 = 35.15;
params.b_NH3 = 2.954e-2;
params.c_NH3 = 4.421e-6;
params.d_NH3 = -6.686e-9;
% Geometría del reactor
D_t = 2.0;                  % Diámetro del reactor industrial [m]
params.Ac = (pi * D_t^2)/4; % Área transversal [m2]
L_lecho = 2.5;              % Longitud de este primer lecho [m]

% Cálculo de alpha para Ergun
gc = 1;                     % Factor de conversión (1 en SI)
beta0 = (params.G * (1 - params.phi) / (params.rho_0 * gc * params.Dp * params.phi^3)) * ...
    (150 * ((1 - params.phi) * params.mu / params.Dp) + 1.75 * params.G);
params.alpha_Ergun = 2 * beta0 / (params.Ac * params.rho_iron * (1 - params.phi) * params.P0);

% 2. CONDICIONES INICIALES DE FLUJO
params.FA0 = 100000; %1          % Flujo inicial de N2 [mol/s]
params.FB0 = 300000; %3         % Flujo inicial de H2 [mol/s]
params.FC0 = 0.05*params.FA0 ; % 5% de FA0           % Flujo inicial de NH3 [mol/s]

% Vector de estado inicial: y0 = [F_N2, F_H2, F_NH3, p, T]
y0 = [params.FA0, params.FB0, params.FC0, 1, params.T0];

% 3. LÍMITES DE INTEGRACIÓN (Masa total de catalizador)
W_total = params.Ac * L_lecho * params.rho_iron * (1 - params.phi);
W_span = [0 W_total];

% 4. SOLUCIÓN NUMÉRICA
opciones = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[W_out, y_out] = ode15s(@(W, y) funcion_reactor(W, y, params), W_span, y0, opciones);

% 5. EXTRACCIÓN Y POST-PROCESAMIENTO
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
plot(W_out, P_real/101325, 'm-', 'LineWidth', 1.5); % Graficado en atm para legibilidad
title('Caída de Presión'); xlabel('W [kg]'); ylabel('P [atm]'); grid on;

subplot(2,2,4);
plot(W_out, X_N2 * 100, 'c-', 'LineWidth', 1.5);
title('Conversión N_2'); xlabel('W [kg]'); ylabel('X [%]'); grid on;

%% MODELADO %%

function dydW = funcion_reactor(~, y, params)
% Inicializar el vector columna de derivadas
dydW = zeros(5, 1);

% 1. Desempaquetar el vector de estado 'y'
F_N2  = y(1);
F_H2  = y(2);
F_NH3 = y(3);
p     = y(4); % Presión reducida (P/P0)
T     = y(5); % Temperatura actual [K]

% 2. Cálculos hidrodinámicos y de flujo
F_T = F_N2 + F_H2 + F_NH3;
F_T0 = params.FA0 + params.FB0 + params.FC0; % Flujo total inicial

% Presión real local en Pascales (o la unidad de tu P0)
P_local = p * params.P0; 

% Presiones parciales/locales
P_local_atm = (p * params.P0) / 101325; 
% Presiones parciales en ATMÓSFERAS
P_N2  = (F_N2 / F_T) * P_local_atm;
P_H2  = (F_H2 / F_T) * P_local_atm;
P_NH3 = (F_NH3 / F_T) * P_local_atm;

% Trampa numérica: Evitar división por cero en Temkin-Pyzhev al inicio del lecho
if P_NH3 < 1e-10
    P_NH3 = 1e-10; 
end

% 3. Propiedades termodinámicas y cinéticas dependientes de la temperatura (T)
% Constante cinética inversa (Arrhenius)
k_reversa = params.A * exp(-params.E / (params.R * T));

% Constante de equilibrio termodinámico (Ka)
log10_Ka = -2.691122 * log10(T) ...
    - 5.51965e-5 * T ...
    + 1.848863e-6 * T^2 ...
    + (2001.3 / T) ...
    + 2.689;
Ka = 10^(log10_Ka);

% 4. Ley de velocidad clásica (Temkin-Pyzhev)
% Tasa de reacción respecto a la formación de NH3
r_NH3 = k_reversa * ( Ka^2 * P_N2 * (P_H2^3 / P_NH3^2)^params.alpha_TP - ...
    (P_NH3^2 / P_H2^3)^(1 - params.alpha_TP) );

% Leyes de velocidad relativa por estequiometría
r_N2 = -0.5 * r_NH3;
r_H2 = -1.5 * r_NH3;

% 5. Ensamblaje del sistema de ecuaciones diferenciales (dydW)
% Balances molares
dydW(1) = r_N2;
dydW(2) = r_H2;
dydW(3) = r_NH3;

% Balance de momentum (Ecuación de Ergun acoplada con T y Flujo)
dydW(4) = (-params.alpha_Ergun / (2 * p)) * (F_T / F_T0) * (T / params.T0);

% Balance de energía (¡Por completar!)
% Necesitas definir los Cp de cada gas y acoplar el calor de reacción
% dydW(5) = ... 

% Balance de energía
% a. Capacidades caloríficas individuales [J/(mol*K)]
Cp_N2  = params.a_N2  + params.b_N2*T  + params.c_N2*T^2  + params.d_N2*T^3;
Cp_H2  = params.a_H2  + params.b_H2*T  + params.c_H2*T^2  + params.d_H2*T^3;
Cp_NH3 = params.a_NH3 + params.b_NH3*T + params.c_NH3*T^2 + params.d_NH3*T^3;

% b. Variación de la capacidad calorífica (Delta Cp) para la reacción base
% Basado en la estequiometría de 1 mol de N2: N2 + 3H2 -> 2NH3
Delta_Cp = 2*Cp_NH3 - Cp_N2 - 3*Cp_H2;

% c. Calor de reacción a la temperatura T (Ley de Kirchhoff)
% Se asume T_ref = 298.15 K para el params.H (-92400 J/mol N2)
T_ref = 298.15;
DeltaH_Rx = params.H + Delta_Cp * (T - T_ref);

% d. Sumatoria de flujos molares por sus capacidades caloríficas
sum_Fi_Cpi = F_N2*Cp_N2 + F_H2*Cp_H2 + F_NH3*Cp_NH3;

% e. Ecuación diferencial del balance de energía (dT/dW)
% r_N2 es negativo (consumo) y DeltaH_Rx es negativo (exotérmico).
% El producto resultante es positivo, lo que incrementa la temperatura.
dydW(5) = (r_N2 * DeltaH_Rx) / sum_Fi_Cpi;
end

