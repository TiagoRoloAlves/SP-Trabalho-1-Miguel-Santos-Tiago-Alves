% Constants
R = 8.314; % Universal gas constant [J/mol-K]

% Temperature range and equivalence ratio range
phi_range = 0.5:0.01:2.0;

% Initial temperature
T0 = 298.15; % [K]

% NASA coefficients for reactants and products
coeffs_CH4 = [0.849032208, 0.141650903E-01, -0.680005317E-04, 0.196106159E-06, -0.286146342E-10, -0.102810669E+05, 0.155875849E+01];
coeffs_O2  = [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00];
coeffs_CO2 = [2.35677352E+00, 8.98459677E-03, -7.12356269E-06, 2.45919022E-09, -1.43699548E-13, -4.83719697E+04, 9.90105222E+00];
coeffs_H2O = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00];
coeffs_H2  = [2.34433112E+00, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01];
coeffs_CO  = [3.57953347E+00, -6.10353680E-04, 1.01681433E-06, -9.07005884E-10, 3.22944658E-13, -1.43440860E+04, 3.50840928E+00];

% Preallocate storage for results
T_ad = zeros(size(phi_range));

% Enthalpy_balance_dissociation function
function res = enthalpy_balance_dissociation(T, H_reactants, phi, coeffs_CH4, coeffs_O2, coeffs_CO2, coeffs_H2O, coeffs_H2, coeffs_CO)
    % Estimate equilibrium fractions using simplified equilibrium constants
    alpha = equilibrium_fraction_CO2(T); % Fraction of CO2 dissociated to CO
    beta = equilibrium_fraction_H2O(T);  % Fraction of H2O dissociated to H2

    % Product composition
    if phi <= 1
        % Lean or stoichiometric
        n_CO2 = (1 - alpha); % CO2 not dissociated
        n_H2O = 2 * phi * (1 - beta); % H2O not dissociated
        n_CO = alpha; % CO from dissociation
        n_H2 = beta;  % H2 from dissociation
        n_O2 = (2 * (1 - phi)); % Excess oxygen
    else
        % Rich mixture
        n_CO2 = (1 - alpha); % CO2 not dissociated
        n_H2O = 2 * (1 - beta); % H2O not dissociated
        n_CO = alpha; % CO from dissociation
        n_H2 = beta + (phi - 1) * 2; % H2 from dissociation and excess fuel
        n_O2 = 0; % No excess oxygen
    end

    % Product enthalpy
    h_CO2 = enthalpy(T, coeffs_CO2);
    h_H2O = enthalpy(T, coeffs_H2O);
    h_CO = enthalpy(T, coeffs_CO);
    h_H2 = enthalpy(T, coeffs_H2);
    h_O2 = enthalpy(T, coeffs_O2);
    H_products = n_CO2 * h_CO2 + n_H2O * h_H2O + n_CO * h_CO + n_H2 * h_H2 + n_O2 * h_O2;

    % Residual (should be zero at correct T)
    res = H_reactants - H_products;
end

% Function to estimate equilibrium fraction of CO2 dissociation (alpha)
function alpha = equilibrium_fraction_CO2(T)
    R = 8.314;  % Universal gas constant [J/mol-K]
    % Simplified temperature-dependent equilibrium constant for CO2 dissociation
    % This is just a generic placeholder expression for demonstration
    K_CO2 = exp(-15000 / (R * T));  % Placeholder: K_CO2 depends on T
    alpha = K_CO2 / (1 + K_CO2);  % Fraction dissociated
end

% Function to estimate equilibrium fraction of H2O dissociation (beta)
function beta = equilibrium_fraction_H2O(T)
    R = 8.314;  % Universal gas constant [J/mol-K]
    % Simplified temperature-dependent equilibrium constant for H2O dissociation
    % This is just a generic placeholder expression for demonstration
    K_H2O = exp(-17000 / (R * T));  % Placeholder: K_H2O depends on T
    beta = K_H2O / (1 + K_H2O);  % Fraction dissociated
end

% Function to compute enthalpy at a given temperature
function h = enthalpy(T, coeffs)
    % Compute enthalpy using NASA polynomial
    R = 8.314; % Universal gas constant [J/mol-K]
    h = R * T * (coeffs(1) + coeffs(2) * T / 2 + coeffs(3) * T^2 / 3 + ...
                 coeffs(4) * T^3 / 4 + coeffs(5) * T^4 / 5 + coeffs(6) / T);
end

% Loop over equivalence ratios
for i = 1:length(phi_range)
    phi = phi_range(i);

    % Reactant moles (scaled by equivalence ratio)
    n_CH4 = phi; % Methane
    n_O2 = 2; % Oxygen (stoichiometric with CH4)

    % Air composition (21% O2, 79% N2), assume N2 constant
    n_N2 = n_O2 * (79 / 21);

    % Initial enthalpy of reactants
    h_CH4 = enthalpy(T0, coeffs_CH4);
    h_O2 = enthalpy(T0, coeffs_O2);
    h_N2 = enthalpy(T0, coeffs_O2); % Approximate using O2 for simplicity
    H_reactants = n_CH4 * h_CH4 + n_O2 * h_O2 + n_N2 * h_N2;

    % Solve for adiabatic flame temperature using enthalpy balance
    T_guess = 2000; % Initial guess
    T_ad(i) = fsolve(@(T) enthalpy_balance_dissociation(T, H_reactants, phi, coeffs_CH4, coeffs_O2, coeffs_CO2, coeffs_H2O, coeffs_H2, coeffs_CO), T_guess);
end

% Plot results
figure;
plot(phi_range, T_ad, 'LineWidth', 1.5);
xlabel('Equivalence Ratio (\phi)');
ylabel('Adiabatic Flame Temperature (K)');
title('Adiabatic Flame Temperature with Product Dissociation vs Equivalence Ratio');
grid on;



