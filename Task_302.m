%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementierung von periodischen Randbedingunten %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% restricting particle coordinates to the simulation box

% x - Teilchenposition in 3D
% x_size - Größe der Simulationsbos (hier: zentrierte Box als
% Vorraussetzung)

% if (x < -x_size * 0.5)
%     x = x + x_size;
% elseif (x >= x_size * 0.5)
%     x = x - x_size
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code-Strukturierung: Implementierung der Routinen eines MD-Programms %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hauptprogramm
% 
% tic
% 
% % Initialisierung der Daten mit np Punkt-Teilchen in einer nd-D Box mit den
% % Boxmaßen d_x,y,z = 30
% 
% [a, sigma, E, m, t, delta_t, n_steps, t_end] = Initialisierung_Parameter(4, 6.4345, 3.7973e4, 72820.75, 0, 5e-21/2.42e-17, 10000);
% % m = 72820.75 a.u (Argon)
% 
% % Teilchen in Box
% [xyz, v, d] = Initialisierung_PBC(30, 30, 30, 10);
% 
% % Kraftberechnung - z. B. LJ_Kraft oder Coulomb ...
% F_0 = LJ_Kraft(xyz, sigma, E);
% 
% % % Zeitintegration - z. B. Velocity-Verlet
% xyz_all = Zeitintegration(xyz, t, delta_t, F_0, n_steps, v, m, sigma, E, d);
% % 
% % Visualisierung als Film in matlab oder output als xyz Datei für vmd
% % Visualisierung(xyz_all, 'Film.avi')
% % Generate_xyz(xyz_all, 'Box.xyz')
% % 
% toc

figure
d_x =15; d_y =15; d_z =15;
x = [d_x, d_x, d_x, d_x, -d_x, -d_x, -d_x, -d_x];
y = [d_y, d_y, -d_y, -d_y, -d_y, -d_y, d_y, d_y];
z = [-d_z, d_z, d_z, -d_z, -d_z, d_z, d_z, -d_z];
% plot3(x,y,z, '-b')
% patch (x,y,z, 'red')
fill3(x,y,z, 'b')

% xmin=0;
% xmax=1;
% ymin=0;
% ymax=1;
% plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin])
% plot3([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin], [z])

%% Funktionen

function [a, sigma, E, m, t, delta_t, n_steps, t_end] = Initialisierung_Parameter (a, sigma, E, m, t, delta_t, n_steps)
t_end = delta_t*n_steps;
end

function [xyz, v, d] = Initialisierung_PBC(d_x, d_y, d_z, np_Teilchen)
% Initialisierung
% zentrierte Box
d = [d_x, d_y, d_z];
% np zufällig verteilte Punkt-Teilchen in den Grenzen der Box
% Initialisierte Zufallszahlen mit 1239465719 für Vergleichbarkeit mit rng

rng(1239465719);
xyz_0 = rand (np_Teilchen,3); % rand: Zufallszahlen zwischen 0 und 1
xyz = bsxfun(@minus, bsxfun(@times, xyz_0, d), 0.5*d);
v = zeros(size(xyz, 1), 3);
end


function F_total = LJ_Kraft(xyz, sigma, E)
% Berechnung der Kräfte
n = length(xyz);
F_neu = zeros(n,3,n);
r = zeros(n,3,n);
F_total = zeros(n,3);
r_betrag = zeros(n,1,n);
for i=1:n
    r(:,:,i) = bsxfun(@minus, xyz(i,:), xyz);
    for j=1:n
        r_betrag(j,1,i) = norm(r(j,:,i));
    end
    F_neu(:,:,i) = bsxfun(@times, (24*E./(r_betrag(:,:,i).^2).*(sigma./r_betrag(:,:,i)).^6.*(1-2*(sigma./r_betrag(:,:,i)).^6)), r(:,:,i));
    F = F_neu;
    F(isnan(F_neu)) = 0;
end
F_total(:,:) = sum(F, 3);
end

function xyz_all = Zeitintegration(xyz, t, delta_t, F_0, n_steps, v, m, sigma, E, d)
% Velocity-Verlet
xyz_all = zeros(length(xyz), 3, n_steps);
iteration = 0;
while t < delta_t*n_steps
    t = t + delta_t;
    iteration = iteration + 1;
    xyz = xyz + delta_t*(v + F_0.*delta_t*0.5/m);
    % neu: Periodische Randbedingungen
    for i = 1:size(xyz, 1)
        for j = 1:size(xyz, 2)
            if xyz(i,j) >= 0.5*d(j)
                xyz(i,j) = xyz(i,j) - d(j);
            elseif xyz(i,j) < -0.5*d(j)
                xyz(i,j) = xyz(i,j) + d(j);
            end
        end
    end
    F_neu = LJ_Kraft(xyz, sigma, E);
    v = v + bsxfun(@rdivide, (F_0+F_neu), 2*m)*delta_t;
    xyz_all(:, :, iteration) = xyz(:, :);
    F_0 = F_neu;
end
end

function Visualisierung(xyz_all, Dateiname)
% Visualisierung in matlab mithilfe von VideoWriter
v = VideoWriter(Dateiname);
open(v);
figure;
for i = 1:size(xyz_all,3)
%     xyz_i = xyz_all(:, :, i);
    clf;
    plot3(xyz_all(:, 1, i), xyz_all(:, 2, i), xyz_all(:, 3, i), 'o', 'MarkerSize', 6);
    title('Moleküldynamik');
    xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;
    M = getframe(gcf);
    writeVideo(v, M);
end
close(v);
end

function Generate_xyz(M, Dateiname)
% Output als xyz Datei
n = size(M, 1);
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