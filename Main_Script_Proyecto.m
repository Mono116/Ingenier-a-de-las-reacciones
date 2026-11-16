%%% PROYECTO DE INGENIERIA DE LAS REACCIONES
%%% Análisis de Sensibilidad y Evaluación Técnico-Económica de un Proceso Haber-Bosch Multilecho

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
% SCRIPT PRINCIPAL: SIMULACIÓN DE REACTORES EN SERIE PROCESO HABER-BOSCH
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
%%%%                        % Densidad del gas en la entrada [kg/m3]
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
D_t = 3.0;                  % Diámetro del reactor industrial [m]
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

% Mapeo estricto de índices cruzados: 
% y(1)=N2, y(2)=H2, y(3)=NH3
term_cruzados = 2 * ( ...
    y(1)*y(2)*sqrt(a(1)*a(2))*(1 - k_N2_H2)  + ... % N2 - H2
    y(1)*y(3)*sqrt(a(1)*a(3))*(1 - k_N2_NH3) + ... % N2 - NH3
    y(2)*y(3)*sqrt(a(2)*a(3))*(1 - k_H2_NH3)   ... % H2 - NH3
);

a_m = term_puros + term_cruzados;
b_m = y(1)*v(1) + y(2)*v(2) + y(3)*v(3);

% 4. Ecuación cúbica de Estado en Z (Peng-Robinson)
A = a_m * params.P0 / ((params.R * params.T0)^2);
B = b_m * params.P0 / (params.R * params.T0);

raices_reales = roots([1, -(1 - B), (A - 3*B^2 - 2*B), -(A*B - B^2 - B^3)]);

% Seleccionar la raíz física real mayor (fase gas)
Z_entrada = max(real(raices_reales(imag(raices_reales) == 0)));
params.Z0 = Z_entrada; % Guardamos Z0 en tu estructura de parámetros
disp(['Factor de compresibilidad inicial (Z0): ', num2str(params.Z0)]);

%%%%% Densidad a partir de Z0
% Calcular primero FT0, y_entrada y llamar a calcularPR para obtener Z0
FT0 = params.FA0 + params.FB0 + params.FC0;
y_entrada = [params.FA0 / FT0, params.FB0 / FT0, params.FC0 / FT0];
[Z0, f_hat_entrada] = calcularPR(params.T0, params.P0, y_entrada, params);
params.Z0 = Z0;

% Peso molecular promedio en la entrada [kg/mol]
MW = [0.0280134; 0.00201588; 0.017031];
MW_prom0 = y_entrada * MW;

% Densidad real inicial [kg/m3]
params.rho_0 = (params.P0 * MW_prom0) / (params.Z0 * params.R * params.T0);
disp(["Densidad real inicial:", params.rho_0]);
% 2. Construir el vector de fracciones molares iniciales [y_N2, y_H2, y_NH3]
y_entrada = [params.FA0 / FT0, params.FB0 / FT0, params.FC0 / FT0];

% 3. Llamar a la función local para obtener Z0 y el vector de fugacidades (f_hat)
[Z0, f_hat_entrada] = calcularPR(params.T0, params.P0, y_entrada, params);

% 4. Guardar Z0 en tu estructura de parámetros para la ecuación de Ergun
params.Z0 = Z0;

% Mostrar resultados iniciales en la ventana de comandos
disp(['Factor de compresibilidad inicial (Z0): ', num2str(Z0)]);
disp('Fugacidades iniciales [atm]:');
disp(['  f_N2  = ', num2str(f_hat_entrada(1)), ' atm']);
disp(['  f_H2  = ', num2str(f_hat_entrada(2)), ' atm']);
disp(['  f_NH3 = ', num2str(f_hat_entrada(3)), ' atm']);

%% SIMULACIONES DE LECHOS %%

% --- SIMULACIÓN DEL LECHO 1 ---
W_lecho1 = 20000; % kg de catalizador
W_span1  = [0, W_lecho1];
[W_out1, y_out1] = ode15s(@(W, y) funcion_reactor(W, y, params), W_span1, y0, opciones);

% Intercooling 1
y0_lecho2 = y_out1(end, :);         % Estado final del gas 1 = Estado inicial gas 2
y0_lecho2(5) = 673;                 % Enfriamiento del gas de salida de nuevo a 673 K

% --- SIMULACIÓN DEL LECHO 2 ---
W_lecho2 = 30000; % kg de catalizador
% La integración inicia donde terminó el lecho anterior
W_span2  = [W_out1(end), W_out1(end) + W_lecho2];
% CORRECCIÓN: Usar y0_lecho2 en lugar de y0
[W_out2, y_out2] = ode15s(@(W, y) funcion_reactor(W, y, params), W_span2, y0_lecho2, opciones);

% Intercooling 2
y0_lecho3 = y_out2(end, :);         % Estado final del gas 2 = estado inicial del gas 3
y0_lecho3(5) = 673; 

% --- SIMULACIÓN DEL LECHO 3 ---
W_lecho3 = 30000; % kg de catalizador
% Limites de integración
W_span3 = [W_out2(end), W_out2(end) + W_lecho3];
% CORRECCIÓN: Usar y0_lecho3 en lugar de y0
[W_out3, y_out3] = ode15s(@(W, y) funcion_reactor(W, y, params), W_span3, y0_lecho3, opciones); 

% Ensamblaje de los perfiles multilecho
W_out = [W_out1; W_out2; W_out3];
y_out = [y_out1; y_out2; y_out3];

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

% 2. Cálculos hidrodinámicos y termodinámicos locales (Peng-Robinson)
F_T = F_N2 + F_H2 + F_NH3;
F_T0 = params.FA0 + params.FB0 + params.FC0;

P_local = p * params.P0; % [Pa]
y_local = [F_N2 / F_T, F_H2 / F_T, F_NH3 / F_T];

% Evaluar Z y fugacidades parciales [atm] en cada paso dW
[Z, f_hat] = calcularPR(T, P_local, y_local, params);

f_N2  = f_hat(1);
f_H2  = f_hat(2);
f_NH3 = f_hat(3);

% Trampa numérica para evitar división por cero
if f_NH3 < 1e-10
    f_NH3 = 1e-10; 
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

% 4. Ley de velocidad de Temkin-Pyzhev usando FUGACIDADES
r_NH3 = k_reversa * ( Ka^2 * f_N2 * (f_H2^3 / f_NH3^2)^params.alpha_TP - ...
    (f_NH3^2 / f_H2^3)^(1 - params.alpha_TP) );

r_N2 = -0.5 * r_NH3;
r_H2 = -1.5 * r_NH3;

% 5. Balances molares (ESTAS TRES LÍNEAS FALTABAN)
dydW(1) = r_N2;
dydW(2) = r_H2;
dydW(3) = r_NH3;

% Balance de momentum (Ecuación de Ergun corregida con Z / Z0)
dydW(4) = (-params.alpha_Ergun / (2 * p)) * (F_T / F_T0) * (T / params.T0) * (Z / params.Z0);
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

