%%% REACTOR PBR -- INFORME 6 INGENIERÍA DE LAS REACCIONES -- %%%
% INTEGRANTES: Luis Paredes, Juan Buendia, Valentina Villalobos, Sergio Barrios %

clc; clear; close all;

% DATOS %
P0 = 10;        % atm
FA0 = 2;        % mol/s (o mol/min)
CA0 = 0.2;      % mol/dm3
kA = 10;        % dm9/(kg*mol2*s)
alpha = 0.05;   % 1/kg

% Alimentación estequiométrica %
CB0 = CA0 * 2;
Theta_B = CB0 / CA0;

% Sistema de ecuaciones diferenciales acopladas con efecto de presión en la cinética
PBR = @(W, Y) [
    (kA * CA0^3 / FA0) * (1 - Y(1)) * (Theta_B - 2*Y(1)).^2 * (Y(2)^3) ; % dX/dW
    -alpha / (2 * Y(2))                                                 % dy/dW
    ];

% Configurar evento para
% detener la integración cuando y = 0 (Presión = 0)%
opciones = odeset('Events', @eventoPresionCero);

% Condiciones iniciales e intervalo %
Wspan = [0 30];  % Masa de catalizador de 0 a W_final
Y0    = [0; 1];   % En W=0: X=0 (conversión inicial) e y=1 (P=P0)

% Resolver el sistema %
[W, Y] = ode45(PBR, Wspan, Y0, opciones);

X = Y(:, 1); % Extraer conversión
y = Y(:, 2); % Extraer presión relativa (P/P0)

% Graficar ambos perfiles %
figure;
subplot(2,1,1);
plot(W, X, 'b-', 'LineWidth', 1.5);
grid on;
ylabel('Conversión (X)');
title('Perfil de Conversión y Caída de Presión vs Masa de Catalizador (W)');

subplot(2,1,2);
plot(W, y, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Masa de Catalizador, W (kg)');
ylabel('Presión Relativa (P/P0)');

% Función local para detectar cuAndo y llega a 0
function [value, isterminal, direction] = eventoPresionCero(~, Y)
value      = Y(2) - 0.01; % Se detiene cuando y <= 0.01 (1% de la presión inicial)%
isterminal = 1;           % Detener la integración de ode45 %
direction  = -1;          % Detectar solo cuando va decreciendo %
end
%% RELACIÓN DE FLUJO VS W %%

% Extraer la presión relativa de la matriz Y (segunda columna)
y = Y(:, 2); 

% Calcular v/v0 
% No hay variacion en los moles F_T / F_T0 = 1 %
v_ratio = 1 ./ y;

% Graficar la relación de flujo volumetrico %
figure;
plot(W, v_ratio, 'r-', 'LineWidth', 2);
xlabel('Masa de catalizador, W (kg)');
ylabel('Relación de flujo volumétrico (v / v_0)');
title('Expansión del gas debido a la caída de presión');
grid on;