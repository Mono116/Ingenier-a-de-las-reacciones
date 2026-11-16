function [Z, f_hat] = calcularPR(T, P, y, params)
    y = y(:); % (fuerza a 'y' a ser vector columna 3x1)
    % Propiedades críticas y factores acéntricos (1: N2, 2: H2, 3: NH3)
    Tc = [126.2; 33.16; 405.4];                    % [K]
    Pc = [33.5356526; 13.13890945; 112.0453984703] * 101325; % [Pa]
    omega = [0.037; -0.215993; 0.256];
    
    % Matriz de Parámetros de Interacción Binaria (kij)
    k_ij = [
        0.0000,  0.0711,  0.2193;  
        0.0711,  0.0000,  0.0000;  
        0.2193,  0.0000,  0.0000   
        ];
    
    R = params.R;
    
    % Parámetros puros a_i y b_i evaluados a la temperatura local T
    a_i = zeros(3, 1);
    b_i = zeros(3, 1);
    
    for i = 1:3
        kappa = 0.37464 + 1.54226 * omega(i) - 0.26992 * omega(i)^2;
        alpha = (1 + kappa * (1 - sqrt(T / Tc(i))))^2;
        a_i(i) = 0.45724 * ((R^2 * Tc(i)^2) / Pc(i)) * alpha;
        b_i(i) = 0.07780 * R * Tc(i) / Pc(i);
    end
    
    % Reglas de mezcla de van der Waals
    a_m = 0;
    for i = 1:3
        for j = 1:3
            a_ij = sqrt(a_i(i) * a_i(j)) * (1 - k_ij(i, j));
            a_m = a_m + y(i) * y(j) * a_ij;
        end
    end
    b_m = sum(y .* b_i);
    
    % Coeficientes adimensionales A y B
    A = a_m * P / (R * T)^2;
    B = b_m * P / (R * T);
    
    % Resolución de la Ecuación Cúbica de Z
    coeffs = [1, -(1 - B), (A - 3*B^2 - 2*B), -(A*B - B^2 - B^3)];
    raices = roots(coeffs);
    Z = max(real(raices(imag(raices) == 0)));
    
    % Coeficientes de fugacidad (phi_hat)
    ln_phi = zeros(3, 1);
    term_common = A / (2 * sqrt(2) * B);
    
    for i = 1:3
        sum_ay = 0;
        for j = 1:3
            a_ij = sqrt(a_i(i) * a_i(j)) * (1 - k_ij(i, j));
            sum_ay = sum_ay + y(j) * a_ij;
        end
    
        term1 = (b_i(i) / b_m) * (Z - 1);
        term2 = -log(Z - B);
        term3 = -term_common * ((2 * sum_ay / a_m) - (b_i(i) / b_m)) * ...
            log((Z + (1 + sqrt(2)) * B) / (Z + (1 - sqrt(2)) * B));
    
        ln_phi(i) = term1 + term2 + term3;
    end
    
    phi_hat = exp(ln_phi);
    
    % Fugacidades parciales en atmósferas
    P_atm = P / 101325;
    f_hat = phi_hat .* y * P_atm; 
end