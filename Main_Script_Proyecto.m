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


% DATOS INICIALES Y PARÁMETEROS

% Parámetros termodinámicos
R = 8.314;  % Constante de los gases ideales [J/mol K]
T0 = 400 + 273.15; % Temperatura inicial [K]
P0 = 150 * 101325; % Presión inicial [Pa]
H = -92400; % J/mol
% Parametros de Temkin-Pyzhev
alpha_TP = 0.5;
% Parámetros de Arrhenius
E = 96232; % Energía de activación [J/mol]
A = 9.02*10^8; % Factor A

% Parámetros de la ecuación de Ergún
phi = 0.45; % Porosidad del hierro promovido
gc = 1; % Factor de conversión
Dp = 300 / 10^6; % Diametro de partícula [m]
rho_b = 2500; % Densidad bulk [kg/m3]
mu = 2.0*10^-5; % viscocidad del gas [Pa*s]
rho_0 = 26.7; % Densidad del gas en la entrada [kg/m3]
v = 1; % Velocidad superficial [m/s]
G = 26.7; % Flujo másico superficial [m2/s]
rho_iron = rho_b/(1-phi);   % Dnsidad del hierro
beta0 = (G*(1-phi)/(rho_0*gc*Dp*phi^3))*(150*((1-phi)*mu/Dp)+1.75*G);

D_t = 2 / 39.37; %Diametro de tuberia [m]
Ac = (pi*D_t^2)/4; % Área de la sección transversal [m2]

alpha_Ergun = 2*beta0/(Ac*rho_iron*(1-phi)*P0);

% Datos de la reaccion
a = 1;
b = 3;
c = 2;

% Datos de flujo

FA0 = 100;      %Flujo de A [mol/s]
FB0 = 300;      %Flujo de B [mol/s]
FC0 = 0;
