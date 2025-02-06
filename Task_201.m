% Lennard-Jones Potential
% Parameter nach White, J. Chem. Phys. 1999, 111, 9352–9356
m = 12;
n = 6; 
sigma = 0.3345; % /nm, V(r=sigma) = 0
E = -125.7; % /k_B; /K, Potentialmulde
a = 4;

% r_m = nthroot(2, 6)*sigma

r = 0.32 : 0.01 : 0.88;
for i = 1:length(r)
    V_LJ(i) = a*E*((sigma/r(i))^n-(sigma/r(i))^m);
end
figure;
plot(r, V_LJ, '-b','LineWidth', 2);
title(['Lennard-Jones-Potential']);
xlabel('r / nm'); ylabel('E/k_B /K'); grid on;

