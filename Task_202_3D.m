%% Initialisierung
% Parameter
sigma = 1;
E = 1;
m = 1; % Masse für H-Atom, 1.00784 u

% Ausganbsbedingungen und Speicherzuordnung
t = 0;
delta_t = 0.01;
t_end = delta_t*100;
iteration = 0;
iteration_total = t_end/delta_t;

%% Kraftberechnung und Moleküldynamik

x = 0 : 1.5*sigma : 1.5*sigma*5;
y = 0 : 1.5*sigma : 1.5*sigma*5;
z = 0 : 1.5*sigma : 1.5*sigma*5;

[X,Y,Z] = meshgrid(x,y,z);
Positions = [X(:), Y(:), Z(:)];
n = length(Positions);

v = zeros(n,3);
xyz_positions = zeros(n,3,int32(iteration_total));

F = LJ_Kraft(Positions, sigma, E);

while t < t_end
    t = t + delta_t;
    iteration = iteration + 1;
    Positions = Positions + delta_t*(v + F.*delta_t*0.5/m);
    F_neu = LJ_Kraft(Positions, sigma, E);
    v = v + bsxfun(@rdivide, (F+F_neu), 2*m)*delta_t;
    xyz_positions(:, :, iteration) = Positions(:, :);
    F = F_neu;
end

% Generate_xyz(xyz_positions, '202_3D.xyz')

%% Skalierungsverhalten für unterschiedlich viele Atome


a = 1:1:9;
Rechenzeit = zeros(length(a),1);

for j = 1:length(a) % a = Anzahl Atome in einer Dimension
    x = 0 : 1.5*sigma : 1.5*sigma*j;
    y = 0 : 1.5*sigma : 1.5*sigma*j;
    z = 0 : 1.5*sigma : 1.5*sigma*j;

    [X,Y,Z] = meshgrid(x,y,z);
    Positions = [X(:), Y(:), Z(:)];
    n = length(Positions);

    v = zeros(n,3);
    xyz_positions = zeros(n,3,int32(iteration_total));

    F = LJ_Kraft(Positions, sigma, E);

    tstart = tic;
    while t < t_end
        t = t + delta_t;
        iteration = iteration + 1;
        Positions = Positions + delta_t*(v + F.*delta_t*0.5/m);
        F_neu = LJ_Kraft(Positions, sigma, E);
        v = v + bsxfun(@rdivide, (F+F_neu), 2*m)*delta_t;
        xyz_positions(:, :, iteration) = Positions(:, :);
        F = F_neu;
    end
    Rechenzeit(j) = toc;
end

figure;
plot(a,Rechenzeit,'-o', 'Linewidth', 2);
title(['Skalierungsverhalten in 3D']);
xlabel('Anzahl Atome pro Dimension'); ylabel('Zeit / ms'); grid on;

%% Funktionen

function F_total = LJ_Kraft(x, sigma, E)
% Berechnung der Kräfte
n = length(x);
F_neu = zeros(n,3,n);
r = zeros(n,3,n);
F_total = zeros(n,3);
r_betrag = zeros(n,1,n);
for i=1:n
    r(:,:,i) = bsxfun(@minus, x(i,:), x);
    for j=1:n
        r_betrag(j,1,i) = norm(r(j,:,i));
    end
    F_neu(:,:,i) = bsxfun(@times, (24*E./(r_betrag(:,:,i).^2).*(sigma./r_betrag(:,:,i)).^6.*(1-2*(sigma./r_betrag(:,:,i)).^6)), r(:,:,i));
    F = F_neu;
    F(isnan(F_neu)) = 0;
end
F_total(:,:) = sum(F, 3);
end

function Generate_xyz(M, Dateiname)
% Output als xyz Datei
n = length(M);
Atomanzahl = repmat({'H'},n,1);
fileID = fopen(Dateiname, 'w'); 

for i=1:size(M,3)
    fprintf(fileID, '%d\n', n);
    fprintf(fileID, 'Moleküldynamik_3D\n');
    for j = 1:n
        fprintf(fileID, '%s %6.4f %6.4f %6.4f\n', Atomanzahl{j}, M(j,1,i), M(j,2,i), M(j,3,i));
    end
end

fclose(fileID);
end