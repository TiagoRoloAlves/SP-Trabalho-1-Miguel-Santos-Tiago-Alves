% Constants and Input
function AlineaDi_equilibrium()
  R = 8.314;  % Universal gas constant [J/mol/K]
  P = 1;      % Pressure in atm
  P0 = 1;     % Reference pressure in atm

  % Equivalence ratio range
  phi = 0.5:0.01:2.0;
  T = 2000;   % Temperature in Kelvin (assumed constant for now)

  % Initial species (moles)
  CH4_input = phi;
  O2_input = 2.0;
  N2_input = 7.52;

  % Preallocate arrays for results
  mole_fractions = zeros(length(phi), 6);  % For CO2, H2O, CO, H2, O2, N2

  for i = 1:length(phi)
    phi_i = phi(i);

    % Initial guess for moles of products
    n = [1, 1, 0, 0, 0, 7.52];  % [n_CO2, n_H2O, n_CO, n_H2, n_O2, n_N2]

    % Solve equilibrium equations using fsolve
    n_eq = fsolve(@(n) equilibrium_equations_kn(n, CH4_input(i), O2_input, N2_input, P, T, P0), n);

    % Calculate mole fractions
    n_tot = sum(n_eq);
    mole_fractions(i, :) = n_eq / n_tot;
  end

  % Plot mole fractions
  plot(phi, mole_fractions(:, 1), 'r', 'DisplayName', 'CO2');
  hold on;
  plot(phi, mole_fractions(:, 2), 'b', 'DisplayName', 'H2O');
  plot(phi, mole_fractions(:, 3), 'g', 'DisplayName', 'CO');
  plot(phi, mole_fractions(:, 4), 'm', 'DisplayName', 'H2');
  plot(phi, mole_fractions(:, 5), 'k', 'DisplayName', 'O2');
  plot(phi, mole_fractions(:, 6), 'c', 'DisplayName', 'N2');
  hold off;
  xlabel('\phi (Equivalence Ratio)');
  ylabel('Mole Fraction');
  title('Equilibrium Mole Fractions vs Equivalence Ratio');
  legend('show');
  grid on;

end

% Equilibrium equations
function F = equilibrium_equations_kn(n, CH4_input, O2_input, N2_input, P, T, P0)
  R = 8.314;  % Universal gas constant [J/mol/K]

  % Extract variables
  n_CO2 = n(1);
  n_H2O = n(2);
  n_CO = n(3);
  n_H2 = n(4);
  n_O2 = n(5);
  n_N2 = n(6);

  % Total moles
  n_tot = sum(n);

  % Equilibrium constants (Kn)
  Kn_H2O = (P / (n_tot * P0))^0.5 * 10^(T/1000 - 1);  % Example temperature-dep. relation
  Kn_CO2 = (P / (n_tot * P0))^0.5 * 10^(T/2000 - 1);  % Example temperature-dep. relation

  % Mass balance equations
  F(1) = CH4_input - (n_CO2 + n_CO);
  F(2) = 2 * CH4_input - (2 * n_H2O + 2 * n_H2);
  F(3) = 2 * O2_input - (2 * n_O2 + n_CO2 + n_H2O + n_CO);
  F(4) = N2_input - n_N2;

  % Equilibrium relations
  F(5) = Kn_H2O - (n_H2 * n_O2^(0.5)) / n_H2O;
  F(6) = Kn_CO2 - (n_CO * n_O2^(0.5)) / n_CO2;
end


AlineaDi_equilibrium()

