% Compute the eigenvalues

K_SB = 0.875;
k_DR = 0.0005;
k_BR = 1;
k_BI = 1;
k_IB = 0.1;
alpha = k_BI - k_IB;

% D can take the values 0; 0.2; 0.3; 0.34; 0.38 ~0.4103036; 1.
D = 0;

% Define Functions
curve = @(B) -alpha * B + K_SB * (exp(2 * (B.^2 + k_BR * B + k_DR * D)) - 1) /...
    (exp(2 * (B.^2 + k_BR * B + k_DR * D)) + 1);

linear_R = @(B) k_BR * B + k_DR * D;

b = -0.0001:0.00001:1.0001;
out = curve(b);

% Plot
figure;
plot(b, out);
xlabel('B');
ylabel('Curve(B)');
title('Curve Plot');

ind = find(out < 0);
jump_ind = find(diff(ind) > 1);

b1 = 0.5 * (b(ind(1)) + b(ind(1) - 1));
b2 = 0.5 * (b(ind(jump_ind) + 1) + b(ind(jump_ind)));
b3 = 0.5 * (b(ind(jump_ind + 1) - 1) + b(ind(jump_ind + 1)));

r1 = linear_R(b1);
r2 = linear_R(b2);
r3 = linear_R(b3);

% Jacobian Calculation
Jac = @(b, r) [
     -alpha + 7 * b .* exp(2 * (b.^2 + r)) ./ ((exp(2 * (b.^2 + r)) + 1).^2),7 * exp(2 * (b.^2 + r)) ./ (2 * (exp(2 * (b.^2 + r)) + 1).^2);
    1, -1];

eig1 = eig(Jac(b1, r1));

eig2 = eig(Jac(b2, r2));

eig3 = eig(Jac(b3, r3));

fprintf('Fixed point 1 is (B, R) = (%f, %f) and the eigenvalues at this point are %f, %f\n', b1, r1, real(eig1(1)), real(eig1(2)));
fprintf('Fixed point 2 is (B, R) = (%f, %f) and the eigenvalues at this point are %f, %f\n', b2, r2, real(eig2(1)), real(eig2(2)));
fprintf('Fixed point 3 is (B, R) = (%f, %f) and the eigenvalues at this point are %f, %f\n', b3, r3, real(eig3(1)), real(eig3(2)));
