% Constants
h_CO2 = -393.52; % Enthalpy of formation of CO2 (kJ/mol)
h_H2O = -241.83; % Enthalpy of formation of H2O (kJ/mol)
h_N2 = 0;        % Enthalpy of formation of N2 (kJ/mol)
h_O2 = 0;        % Enthalpy of formation of O2 (kJ/mol)
h_CH4 = -74.85;  % Enthalpy of formation of CH4 (kJ/mol)

MW_CH4 = 16.04;  % Molecular weight of CH4 (kg/kmol)
N2_O2_ratio = 79 / 21; % Mole ratio of N2 to O2 in air

% Equivalence ratio range
phi_values = 0.5:0.01:2.0;

% Initialize results
combustion_enthalpy_per_kmol = [];
combustion_enthalpy_per_kg = [];

% Dissociation constants (simple model)
alpha = 0.1;  % Fraction of CO2 dissociated to CO
beta = 0.1;   % Fraction of H2O dissociated to H2

% Loop through equivalence ratios
for phi = phi_values
    % Reactant moles
    moles_CH4 = 1;
    moles_O2 = 2 / phi;
    moles_N2 = moles_O2 * N2_O2_ratio;

    % Product composition based on equivalence ratio
    if phi < 1
        % Lean combustion
        moles_CO2 = (1 - alpha);    % Partially dissociated CO2
        moles_H2O = 2 * (1 - beta); % Partially dissociated H2O
        moles_O2_product = (2 - 2 * phi);
        moles_N2_product = moles_N2;
        moles_CH4_unburnt = 0;
    else
        % Rich combustion
        moles_CO2 = (1 - alpha);    % Partially dissociated CO2
        moles_H2O = 2 * (1 - beta); % Partially dissociated H2O
        moles_O2_product = 0;
        moles_N2_product = moles_N2;
        moles_CH4_unburnt = (phi - 1);  % Unburnt CH4 for rich mixture
    end

    % Combustion enthalpy calculation
    H_products = moles_CO2 * h_CO2 + moles_H2O * h_H2O + moles_N2_product * h_N2 + ...
                 moles_O2_product * h_O2 + moles_CH4_unburnt * h_CH4;
    H_reactants = moles_CH4 * h_CH4 + moles_O2 * h_O2 + moles_N2 * h_N2;
    delta_H = H_products - H_reactants; % Combustion enthalpy per kmol

    % Store results
    combustion_enthalpy_per_kmol = [combustion_enthalpy_per_kmol, delta_H];
    combustion_enthalpy_per_kg = [combustion_enthalpy_per_kg, delta_H / MW_CH4];
end

% Plot results
figure;
plot(phi_values, combustion_enthalpy_per_kmol, 'b-', 'LineWidth', 2);
hold on;
plot(phi_values, combustion_enthalpy_per_kg, 'r--', 'LineWidth', 2);
xlabel('Equivalence Ratio (\phi)');
ylabel('Combustion Enthalpy');
legend('\DeltaH (kJ/kmol)', '\DeltaH (kJ/kg)');
title('Combustion Enthalpy vs Equivalence Ratio (with Dissociation)');
grid on;

