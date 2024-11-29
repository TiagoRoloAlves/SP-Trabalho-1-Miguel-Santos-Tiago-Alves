% Equivalence ratio range
phi = 0.5:0.01:2.0;

% Mole fractions
chi_CO2 = zeros(size(phi));
chi_H2O = zeros(size(phi));
chi_O2 = zeros(size(phi));
chi_N2 = zeros(size(phi));

% Mole fractions based on phi
for i = 1:length(phi)
    if phi(i) < 1  % Lean mixture
        chi_CO2(i) = 1 / (10.52 - 2 * phi(i));
        chi_H2O(i) = chi_CO2(i);
        chi_O2(i) = (2 - 2 * phi(i)) / (10.52 - 2 * phi(i));
        chi_N2(i) = 7.52 / (10.52 - 2 * phi(i));
   
    else  % Rich mixture
        chi_CO2(i) = 1 / (7.52 + phi(i));
        chi_H2O(i) = chi_CO2(i);
        chi_O2(i) = 0;
        chi_N2(i) = 7.52 / (7.52 + phi(i));
    end
end

% Plot mole fractions
figure;
plot(phi, chi_CO2, 'r-', 'LineWidth', 1.5); hold on;
plot(phi, chi_H2O, 'b--', 'LineWidth', 1.5);
plot(phi, chi_N2, 'g:', 'LineWidth', 1.5);
if any(phi < 1)  % Include O2 only if phi < 1 exists
    plot(phi, chi_O2, 'k-.', 'LineWidth', 1.5);
end

% labels and legend
xlabel('Equivalence Ratio (\phi)', 'FontSize', 12);
ylabel('Mole Fraction', 'FontSize', 12);
title('Mole Fractions of Combustion Products', 'FontSize', 14);
legend({'CO2', 'H2O', 'N2', 'O2'}, 'Location', 'best', 'FontSize', 10);
grid on;

% axis limits for better readability
xlim([min(phi), max(phi)]);
ylim([0, 1]);

