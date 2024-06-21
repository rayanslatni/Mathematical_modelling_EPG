%Find lienar R and non linear R
% The intersections of the curves below define the fixed points of the system
%In the code we solve the system of liner R and non linear R in this way we
%obtain the 

%FIND THE INTERSECTION BETWEEN LINEAR AND NN LINEAR R
alpha = 7 / (8 * 0.9);
k_DR = 0.0005;
k_BR = 1;
D = 0;

% Define the nonlinear function (Eq. S19)
nonlinear_R = @(B) -B.^2 + 1/2 * log((alpha + B) / (alpha - B));

% Define the linear function (Eq. S12)
linear_R = @(B, D) k_BR * B + k_DR * D;

%Define the polynomial obtian from their intersection (Eq. S25)
poli = @(B) 2*B.^3 +B.^2 -2*B*alpha.^2 +alpha-alpha.^2;

% Define the range for B (0 to 0.95 with a step of 0.01)
B = 0:0.01:0.99;

% Evaluate the functions
r_nonlinear = nonlinear_R(B);
r_linear = linear_R(B, D);
polinomial = poli(B);


% Solve the polynomial equation to find the intersections
coefficients = [2, 1, -2 * alpha^2, alpha - alpha^2]; %coeff of the polinomial
solutions = roots(coefficients);

% Filter positive solutions
positive_id = solutions >= 0;
solutions_positive = solutions(positive_id);

% Calculate corresponding values for B, R, and D
B1 = solutions_positive(1);
R1 = nonlinear_R(B1);
D1 = (R1 - B1) / k_DR;

B2 = solutions_positive(2);
R2 = nonlinear_R(B2);
D2 = (R2 - B2) / k_DR;


disp(['B1 = ', num2str(B1), ' gives R1 = ', num2str(R1), ' gives D1 = ', num2str(D1)]);
disp(['B2 = ', num2str(B2), ' gives R2 = ', num2str(R2), ' gives D2 = ', num2str(D2)]);

% Plot the functions
figure;
plot(B, r_nonlinear, 'LineWidth', 2, 'DisplayName', 'Nonlinear R');
hold on;
plot(B, r_linear, 'LineWidth', 2, 'DisplayName', 'Linear R');
plot(B, polinomial, 'LineWidth',2, 'DisplayName','Polinomial' );
plot([0.41 0.41],ylim, 'r--', 'LineWidth', 1.5, 'DisplayName','Neurotoxicity threshold');  
xlabel('B');
ylabel('R');
title('Intersection between linear R and non linear R')
legend('show');
grid on;

 