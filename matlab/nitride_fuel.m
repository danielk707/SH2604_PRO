lambda_UN = @(T) 1.864 .* T .^(0.361);
lambda_PuN= @(T) 7.446 + 0.663E-2 .* T - 1.55E-6 .* T.^2;
lambda_Av = @(T) 0.8*lambda_UN(T) + 0.2*lambda_PuN(T);

LTE = @(b,T) b(1).*(T-293) + b(2).*(T-293).^2 + b(3).*(T-293).^3;

b_UN = [7.505E-4 1.407E-7 0.0];
b_PuN= [9.973E-4 1.288E-7 -9.445E-12];

rho_Pb = @(T) 11441 - 1.2795 .* T; % kg/m^3

eta_Pb = @(T) 4.55E-4 .* exp (1069./T);
lambda_Pb = @(T) 9.2 + 0.011 .* T;

fprintf('q´_max = %d kW/m\n', 4*pi*quad(lambda_Av, 900, 3100)/1E3);