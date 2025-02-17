%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code-Strukturierung: Implementierung der Routinen eines MD-Programms %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hauptprogramm

tic

% Einlesen der Daten oder Eingabe in Matlab
[a, sigma, E, m, t, delta_t, n_steps, t_end, v] = Initialisierung_Parameter(4, 1, 1, 1, 0, 0.01, 1000);

% [xyz, a, sigma, E, t, delta_t, n_steps, t_end] = Initialisierung('file.txt', a, E, delta_t, n_steps);

xyz = Initialisierung_Matlab;

% Kraftberechnung - z. B. LJ_Kraft oder Coulomb ...
F_0 = LJ_Kraft(xyz, sigma, E);

% Zeitintegration - z. B. Velocity-Verlet
xyz_all = Zeitintegration(xyz, t, delta_t, F_0, n_steps, v, m, sigma, E);

% Visualisierung als Film in matlab oder output als xyz Datei für vmd
Visualisierung(xyz_all, 'Film.avi')
Generate_xyz(xyz_all, 'Dateiname.xyz')

toc

%% Funktionen
function [a, sigma, E, m, t, delta_t, n_steps, t_end] = Initialisierung_Parameter (a, sigma, E, m, t, delta_t, n_steps)
t_end = delta_t*n_steps;
end

function [xyz, v] = Datenimport (Textdatei)
% Initialisierung
xyz = importdata(Textdatei); % Input muss als txt Datei vorliegen mit nur 3 Spalten
v = zeros (length(xyz), 3);
end

function [xyz, v] = Initialisierung_Matlab
% Initialisierung
prompt = "Eingabe der xyz-Positionen als n x 3 Matrix: ";
xyz = input(prompt); % Inputeingabe in Matlab
v = zeros (length(xyz), 3);
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

function xyz_all = Zeitintegration(xyz, t, delta_t, F_0, n_steps, v, m, sigma, E)
% Velocity-Verlet
xyz_all = zeros(length(xyz), 3, n_steps);
iteration = 0;
while t < delta_t*n_steps
    t = t + delta_t;
    iteration = iteration + 1;
    xyz = xyz + delta_t*(v + F_0.*delta_t*0.5/m);
    F_neu = LJ_Kraft(xyz, sigma, E);
    v = v + bsxfun(@rdivide, (F_0+F_neu), 2*m)*delta_t;
    xyz_all(:, :, iteration) = xyz(:, :);
    F_0 = F_neu;
end
end

function Visualisierung(xyz_all, Dateiname)
% Visualisierung in matlab mithilfe von VideoWriter
vw = VideoWriter(Dateiname);
open(vw);
figure;
for i = 1:size(xyz_all,3)
%     xyz_i = xyz_all(:, :, i);
    clf;
    plot3(xyz_all(:, 1, i), xyz_all(:, 2, i), xyz_all(:, 3, i), 'o', 'MarkerSize', 6);
    title('Moleküldynamik');
    xlabel('X'); ylabel('Y'); zlabel('Z'); grid on;
    M = getframe(gcf);
    writeVideo(vw, M);
end
close(vw);
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