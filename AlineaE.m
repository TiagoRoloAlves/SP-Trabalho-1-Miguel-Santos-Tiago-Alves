clear; clc;

import cantera.*;

% Constants
P = 101325; % Atmospheric pressure in Pa
T_initial = 298.15; % Initial temperature in Kelvin
T_standard = 298.15; % Standard conditions temperature in Kelvin
phi_range = 0.5:0.01:2.0; % Equivalence ratio range
species_list = {'CH4', 'O2', 'N2', 'CO2', 'H2O', 'H2', 'CO', 'OH'}; % Relevant species
fuel = 'CH4'; % Fuel
oxidizer = {'O2', 'N2'}; % Oxidizer
molar_mass_fuel = 16.04; % CH4 molar mass in g/mol

% Results storage
adiabatic_flame_temp = zeros(size(phi_range));
product_mole_fractions = zeros(length(phi_range), length(species_list));
combustion_enthalpy_kmol = zeros(size(phi_range)); % Combustion enthalpy per kmol of CH4
combustion_enthalpy_kg = zeros(size(phi_range));   % Combustion enthalpy per kg of CH4

% Loop over equivalence ratios
for i = 1:length(phi_range)
    phi = phi_range(i);
    
    % Fuel and oxidizer stoichiometry
    % Air composition: 21% O2 and 79% N2 (by volume)
    oxidizer_mixture = {oxidizer{1}, 0.21, oxidizer{2}, 0.79};
    
    % Reactant gas mixture at initial conditions
    gas = Solution('gri30.yaml'); % Load GRI-Mech 3.0 mechanism
    reactants = sprintf('%s:%f, %s:%f, %s:%f', ...
        fuel, phi, oxidizer_mixture{1}, 2/phi, oxidizer_mixture{3}, 2/phi*3.76); % Adjust oxidizer for stoichiometry
    set(gas, 'T', T_initial, 'P', P, 'X', reactants);
    
    % Enthalpy of reactants (initial conditions)
    h_reactants = enthalpy_mass(gas);
    
    % Equilibrate gas mixture (for flame temperature and products)
    equilibrate(gas, 'HP'); % Adiabatic flame temp using constant enthalpy and pressure
    
    % Store flame temperature
    adiabatic_flame_temp(i) = temperature(gas);
    
    % Mole fractions of products
    mole_fractions = moleFractions(gas);
    for j = 1:length(species_list)
        product_mole_fractions(i, j) = mole_fractions(speciesIndex(gas, species_list{j}));
    end
    
    % Reset gas to standard conditions
    set(gas, 'T', T_standard, 'P', P, 'X', reactants);
    
    % Enthalpy of reactants (standard conditions)
    h_reactants_std = enthalpy_mass(gas);
    
    % Equilibrate gas mixture at standard conditions
    equilibrate(gas, 'TP'); % Equilibrium at constant temperature and pressure
    h_products_std = enthalpy_mass(gas);
    
    % Combustion enthalpy (standard conditions, per kmol and per kg of CH4)
    combustion_enthalpy_kmol(i) = (h_reactants_std - h_products_std) * molar_mass_fuel; % J/kmol
    combustion_enthalpy_kg(i) = combustion_enthalpy_kmol(i) / molar_mass_fuel * 1e3; % J/kg
end

% Plot 1: Adiabatic Flame Temperature vs Equivalence Ratio
figure;
plot(phi_range, adiabatic_flame_temp, 'LineWidth', 1.5);
xlabel('Equivalence Ratio (\phi)');
ylabel('Adiabatic Flame Temperature (K)');
title('Adiabatic Flame Temperature vs. Equivalence Ratio');
grid on;

% Plot 2: Product Mole Fractions vs Equivalence Ratio
figure;
hold on;
for j = 1:length(species_list)
    plot(phi_range, product_mole_fractions(:, j), 'LineWidth', 1.5, 'DisplayName', species_list{j});
end
xlabel('Equivalence Ratio (\phi)');
ylabel('Mole Fractions');
title('Product Mole Fractions vs. Equivalence Ratio');
legend('Location', 'best');
grid on;
hold off;

% Plot 3: Combustion Enthalpy vs Equivalence Ratio (both kmol and kg in one plot)
figure;
yyaxis left;
plot(phi_range, combustion_enthalpy_kmol, 'b-', 'LineWidth', 1.5);
ylabel('Combustion Enthalpy (J/kmol of CH4)');
yyaxis right;
plot(phi_range, combustion_enthalpy_kg, 'r--', 'LineWidth', 1.5);
ylabel('Combustion Enthalpy (J/kg of CH4)');
xlabel('Equivalence Ratio (\phi)');
title('Combustion Enthalpy vs. Equivalence Ratio');
legend({'J/kmol of CH4', 'J/kg of CH4'}, 'Location', 'best');
grid on;

disp('Simulation complete. Figures displayed.');
