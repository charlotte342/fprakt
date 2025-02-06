% Lennard-Jones Potential
% Parameter nach White, J. Chem. Phys. 1999, 111, 9352–9356
m = 12;
n = 6; 
sigma = 0.3345; % /nm, V(r=sigma) = 0
E = -125.7; % /k_B; /K, Potentialmulde
a = 4;
% 
% % r_m = nthroot(2, 6)*sigma
% 
% r = 0.32 : 0.01 : 0.88;
% for i = 1:length(r)
%     V_LJ(i) = a*E*((sigma/r(i))^n-(sigma/r(i))^m);
% end
% figure;
% plot(r, V_LJ, '-b','LineWidth', 2);
% title(['Lennard-Jones-Potential']);
% xlabel('r / nm'); ylabel('E/k_B /K'); grid on;

% Gradient



% x_1 = [1, 1];
% x_2 = [0, 0];
% r_dir = x_1 - x_2;
% r = norm(r_dir);
% F_LJ = 24*E/r^2*(sigma/r)^6*(1-2*(sigma/r)^6)*r_dir;

% 
spacing = 0.1;
x = -1:0.1:1;
y = -1:0.1:1;
[X,Y] = meshgrid(x,y);

for i=1:length(X)
    F_x = 24.*E/((X(i,:)).^2).*(sigma/X(i,:)).^6.*(1-2*(sigma/X(i,:)).^6);
%     F_y = 24.*E/((y.').^2).*(sigma/y.').^6.*(1-2*(sigma/y.').^6);
end

figure;
quiver(X,Y,F_x,F_y);
hold on
% contour(X,Y,Z)
axis equal; grid on
hold off

% n = 0.25:0.01:0.88;
% x_1 = zeros(length(n), 2);
% for i = 1:length(n)
%     x_1(i,:) = x_1(i,:) + n(i);
%     r_dir(i,:) = x_1(i,:) - x_2(1,:);
%     r(i) = norm(r_dir(i,:));
%     F_LJ(i,:) = 24*E/(r(i)^2)*(sigma/r(i))^6*(1-2*(sigma/r(i))^6)*r_dir(i,:);
% end
% 
% figure;
% quiver(r_dir, F_LJ, '-b','LineWidth', 2);
% title(['Gradient']);
% % xlabel('r / nm'); ylabel('E/k_B /K'); 
% grid on;
