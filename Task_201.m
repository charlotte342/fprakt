% LJ_Potential
m = 12;
n = 6; 
a = 4;

%% Parameter nach White, J. Chem. Phys. 1999, 111, 9352–9356
sigma = 0.3345; % /nm, V(r=sigma) = 0
E = -125.7; % /k_B; /K, Potentialmulde

r = 0.32 : 0.01 : 0.88;
for i = 1:length(r)
    V_LJ(i) = a*E*((sigma/r(i))^n-(sigma/r(i))^m);
    F_LJ(i) = LJ_Grad(E, sigma, 0, r(i));
end
figure;
plot(r, V_LJ, '-b','LineWidth', 2);
hold on
plot(r, F_LJ, '-g', 'LineWidth', 2);
title(['Lennard-Jones-Potential'],['und Lennard-Jones-Gradient']);
xlabel('r / nm'); ylabel('E/k_B /K'); grid on;

%% Lennard-Jones Potential
% Variation von sigma und E
sigma = 0.3 : 0.01 : 0.36;
E = -140 : 5 : -110;
r = 0.32 : 0.01 : 0.88;
V_LJ = zeros(length(r), length(E), length(sigma));

figure;
C = {'-k','-b','-r','-g', '-c', '-m', '-y'};
for i = 1:length(sigma)
    for j = 1:length(E)
        for k = 1:length(r)
            V_LJ(k,j,i) = a*E(j)*((sigma(i)/r(k))^n-(sigma(i)/r(k))^m);
        end
        plot(r, V_LJ(:,j,i), C{i},'LineWidth', 1.5);
        hold on
    end
end

title(['Lennard-Jones-Potential']);
xlabel('r / nm'); ylabel('E/k_B /K'); grid on;

%% LJ-Gradient für ein Atom

% r_m = nthroot(2, 6)*sigma, bei ungefähr 1.12 sigma
% Gradient ist null bei GGW-Abstand r_m
sigma = 3.3345; % jetzt in Angström
E = -125.7;
x_N = [0,0];
x = -10:0.2:10;
y = -10:0.2:10;
[X,Y] = meshgrid(x,y);
xy = [X(:), Y(:)];
n = length(xy);

r_betrag = zeros(n,1);

r = bsxfun(@minus, xy, x_N);
for j=1:n
    r_betrag(j,1) = norm(r(j,:));
end
F = bsxfun(@times, (24*E./(r_betrag.^2).*(sigma./r_betrag).^6.*(1-2*(sigma./r_betrag).^6)), r);
for j = 1:n
    if r_betrag(j,1) < sigma
        F(j,:) = zeros(1,2);
    end
end

U = reshape(F(:,1), length(X), length(X));
V = reshape(F(:,2), length(Y), length(Y));

figure;
quiver(X,Y,U,V);
title(['Lennard-Jones Gradient für ein Atom']);
xlabel('x / A'); ylabel('y / A'); axis equal; grid on;

%% LJ-Gradient für N Atome

sigma = 3.3345; % jetzt in Angström
E = -125.7;
x_N = [-3, 0; 2, 0; 0,7];

x = -10:0.5:10;
y = -10:0.5:10;
[X,Y] = meshgrid(x,y);
xy = [X(:), Y(:)];
n = length(xy);

figure
for i=1:length(x_N)
    plot(x_N(i,1), x_N(i,2), 'b+', 'LineWidth', 2);
    hold on
end

F = zeros(n,2,length(x_N));
r = zeros(n,2,length(x_N));
F_total = zeros(n,2);
r_betrag = zeros(n,1,length(x_N));
for i=1:length(x_N)
    r(:,:,i) = bsxfun(@minus, x_N(i,:), xy);
    for j=1:n
        r_betrag(j,1,i) = norm(r(j,:,i));
    end
    F(:,:,i) = bsxfun(@times, (24*E./(r_betrag(:,:,i).^2).*(sigma./r_betrag(:,:,i)).^6.*(1-2*(sigma./r_betrag(:,:,i)).^6)), r(:,:,i));
    for j = 1:n
        if r_betrag(j,1,i) < sigma
            F(j,:,i) = zeros(1,2);
        end
    end
end
F_total(:,:) = sum(F, 3);

U = reshape(F_total(:,1), length(X), length(X));
V = reshape(F_total(:,2), length(Y), length(Y));

quiver(X,Y,U,V);
title(['Lennard-Jones Gradient für N Atome']);
xlabel('x / A'); ylabel('y / A'); axis equal; grid on;



