% Constants
R = 8.314; % Universal gas constant [J/mol-K]

% Temperature range and equivalence ratio range
phi_range = 0.5:0.01:2.0;

% Initial temperature
T0 = 298.15; % [K]

% Coefficients for NASA polynomials (valid for T >= 1000 K)
coeffs_CH4 = [7.48514950E-02, 1.33909467E-02, -5.73285809E-06, 1.22292535E-09, -1.01815230E-13, -9.46834459E+03, 1.84373180E+01];
coeffs_O2  = [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00];
coeffs_N2  = [2.95257626E+00, 1.39690000E-03, -4.92631603E-07, 7.86010367E-11, -4.60755321E-15, -9.23948688E+02, 5.87188762E+00];
coeffs_CO2 = [3.85746029E+00, 4.41437026E-03, -2.21481404E-06, 5.23490188E-10, -4.72084164E-14, -4.87591660E+04, 2.27163806E+00];
coeffs_H2O = [3.03399249E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00];

% Preallocate result storage
T_ad = zeros(size(phi_range));

% Enthalpy_balance function
function res = enthalpy_balance(T, H_reactants, phi, coeffs_CH4, coeffs_CO2, coeffs_H2O, coeffs_O2, coeffs_N2)
    % Determine products based on equivalence ratio
    if phi <= 1
        % Lean or stoichiometric: complete combustion
        n_CO2 = 1;
        n_H2O = 2;
        n_O2 = 2 / phi - 2;
        n_N2 = (79 / 21) * (2 / phi);
        n_CH4 = 0; % No excess fuel
    else
        % Rich: excess methane
        n_CO2 = 1;
        n_H2O = 2;
        n_O2 = 0; % No excess oxygen
        n_N2 = (79 / 21) * 2; % Air contribution
        n_CH4 = phi - 1; % Excess methane
    end

    % Compute product enthalpy
    h_CO2 = enthalpy(T, coeffs_CO2);
    h_H2O = enthalpy(T, coeffs_H2O);
    h_O2 = enthalpy(T, coeffs_O2);
    h_N2 = enthalpy(T, coeffs_N2);
    h_CH4 = enthalpy(T, coeffs_CH4);
    H_products = n_CO2 * h_CO2 + n_H2O * h_H2O + n_O2 * h_O2 + n_N2 * h_N2 + n_CH4 * h_CH4;

    % Residual (should be zero at correct T)
    res = H_reactants - H_products;
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

    % Reactant moles 
    n_CH4 = 1; % 1 mol of CH4
    n_O2 = (2 / phi) * n_CH4; % Adjusted for equivalence ratio
    n_N2 = (79 / 21) * n_O2; % Air composition (21% O2, 79% N2)

    % Initial enthalpy of reactants
    h_CH4 = enthalpy(T0, coeffs_CH4);
    h_O2 = enthalpy(T0, coeffs_O2);
    h_N2 = enthalpy(T0, coeffs_N2);
    H_reactants = n_CH4 * h_CH4 + n_O2 * h_O2 + n_N2 * h_N2;

    % Solve for adiabatic flame temperature using enthalpy balance
    T_guess = 2000; % Initial guess
    T_ad(i) = fsolve(@(T) enthalpy_balance(T, H_reactants, phi, coeffs_CH4, coeffs_CO2, coeffs_H2O, coeffs_O2, coeffs_N2), T_guess);
end

% Plot results
figure;
plot(phi_range, T_ad, 'LineWidth', 1.5);
xlabel('Equivalence Ratio (\phi)');
ylabel('Adiabatic Flame Temperature (K)');
title('Adiabatic Flame Temperature vs Equivalence Ratio');
grid on;



