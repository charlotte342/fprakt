%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermostat zur Temperaturkontrolle %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Skalierung der Geschwindigkeit 

%% Hauptprogramm

tic

% Initialisierung der Daten mit np Punkt-Teilchen in einer nd-D Box mit den
% Boxmaßen d_x,y,z = 30

% Parameter importieren
[a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter('Parameter_302.txt');

% Teilchen in Box
[xyz, v, d, np_Teilchen] = Initialisierung_PBC(30, 30, 30, 10, sigma);

% Kraftberechnung - z. B. LJ_Kraft oder Coulomb ...
F = LJ_Kraft(xyz, sigma, E, d);

% % Zeitintegration - z. B. Velocity-Verlet
[xyz_all, v, E_kin_all, E_pot_all, E_tot, T_all] = Zeitintegration(xyz, t, delta_t, F, n_steps, v, m, sigma, E, d, np_Teilchen, a, delta_t*1e3); % letzte ist tau fürs Thermostat


% Visualisierung als Film in matlab oder output als xyz Datei für vmd
% Visualisierung(xyz_all, 'Film.avi')
Generate_xyz(xyz_all, 'Thermostat.xyz')

figure;
plot(T_all, '-b', 'LineWidth', 2);
axis equal; grid on;

figure;
plot(E_kin_all, '-b','LineWidth',2);
hold on
plot(E_pot_all, '-r', 'LineWidth', 2);
grid on;

toc
%% Funktionen

function [a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter(Textdatei1)
% Daten importieren
fileID = fopen(Textdatei1);
while ~feof(fileID)
    line = fgets(fileID); % Zeile für Zeile Datei durchgehen
    % 1. Integrator bestimmen
    if contains (line, 'Integrator')
        integrator = strtrim(strrep(line, 'Integrator: ', '')); % string trim, strrep: find and replace substring
    elseif contains(line, 'dt')
        delta_t = str2double(strtrim(strrep(line, 'dt: ', '')));
    elseif contains(line, 'n_steps')
        n_steps = str2double(strtrim(strrep(line, 'n_steps: ', '')));
        t_end = delta_t*n_steps;
    end
    if strcmp(integrator,'Runge-Kutta') == 1 % strcmp: compare strings
        if contains(line, 'Order')
            Order = str2double(strtrim(strrep(line, 'Order: ', '')));
        end
    end
    % 2. Force-Field
    if contains(line, 'Force-Field')
        FF = strtrim(strrep(line, 'Force-Field: ', ''));
    elseif contains(line, 'alpha')
        a = str2double(strtrim(strrep(line, 'alpha: ', '')));
    elseif contains(line, 'sigma')
        sigma = str2double(strtrim(strrep(line, 'sigma: ', '')));
    elseif contains(line, 'E')
        E = str2double(strtrim(strrep(line, 'E: ', '')));
    elseif contains(line, 'Masse')
        m = str2double(strtrim(strrep(line, 'Masse: ', '')));
    end
    
end
t = 0;
fileID = fclose(fileID);
end

function [xyz, v, d, np_Teilchen] = Initialisierung_PBC(d_x, d_y, d_z, np_Teilchen, sigma)
% Initialisierung
% zentrierte Box
d = [d_x, d_y, d_z];
% np zufällig verteilte Punkt-Teilchen in den Grenzen der Box

r = zeros(np_Teilchen,3,np_Teilchen);
r_betrag = zeros(np_Teilchen,1,np_Teilchen);
% Initialisierte Zufallszahlen mit 1239465719 für Vergleichbarkeit mit rng
rng(1239465719);
xyz_0 = rand (np_Teilchen, 3); % rand: Zufallszahlen zwischen 0 und 1
n = length(xyz_0);
xyz = bsxfun(@minus, bsxfun(@times, xyz_0, d), 0.5*d); % Zufallszahlen im Intervall von d

% r_min implementieren
for i = 1:n
    r(:,:,i) = bsxfun(@minus, xyz, xyz(i,:));
    r(:,:,i) = r(:,:,i) - d.* round(r(:,:,i)./d); % PBC: wenn Abstand zu Teilchen in nächster Box kürzer
end
for i = 1:n
    for j = i:n
        r_betrag(j,1,i) = norm(r(j,:,i));
        r_betrag(i,1,j) = r_betrag(j,1,i);
        r_betrag(j,1,j) = sigma;
        while r_betrag(j,1,i) < sigma
            xyz_0(j,:) = rand(1,3);
            xyz(j,:) = bsxfun(@minus, bsxfun(@times, xyz_0(j,:), d), 0.5*d);
            r(:,:,i) = bsxfun(@minus, xyz(i,:), xyz);
            r_betrag(j,1,i) = norm(r(j,:,i));
        end
    end
end     
v = zeros(size(xyz, 1), 3);
end


function F_total = LJ_Kraft(xyz, sigma, E, d)
% Berechnung der Kräfte
n = length(xyz);
F_neu = zeros(n,3,n);
r = zeros(n,3,n);
F_total = zeros(n,3);
r_betrag = zeros(n,1,n);

for i=1:n
    r(:,:,i) = bsxfun(@minus, xyz(i,:), xyz);
    r(:,:,i) = r(:,:,i) - d.*round(r(:,:,i)./d);
    for j=i:n
        r_betrag(j,1,i) = norm(r(j,:,i));
        r_betrag(i,1,j) = r_betrag(j,1,i);
    end
    F_neu(:,:,i) = bsxfun(@times, (24*E./(r_betrag(:,:,i).^2).*(sigma./r_betrag(:,:,i)).^6.*(1-2*(sigma./r_betrag(:,:,i)).^6)), r(:,:,i));
    F = F_neu;
    F(isnan(F_neu)) = 0;
end
F_total(:,:) = sum(F, 3);
end

function E_pot_total = LJ_Pot(xyz, sigma, E, a, d)
n = length(xyz);
E_pot = zeros(n,1,n);
r = zeros(n,3,n);
r_betrag = zeros(n,1,n);
E_pot_total = zeros(n,1);

for i=1:n
    r(:,:,i) = bsxfun(@minus, xyz(i,:), xyz);
    r(:,:,i) = r(:,:,i) - d.*round(r(:,:,i)./d);
    for j=i:n
        r_betrag(j,1,i) = norm(r(j,:,i));
        r_betrag(i,1,j) = r_betrag(j,1,i);
    end
    E_pot(:,:,i) = a*E*((sigma./r_betrag(:,:,i)).^6-(sigma./r_betrag(:,:,i)).^12);
    E_korr = E_pot;
    E_korr(isnan(E_pot)) = 0; % Korrektur für i = j
end
E_pot_total(:,:) = sum(E_korr, 3);
end


function [xyz_all, v, E_kin_all, E_pot_all, E_tot, T_all] = Zeitintegration(xyz, t, delta_t, F_0, n_steps, v, m, sigma, E, d, np_Teilchen, a, tau)
% Velocity-Verlet
xyz_all = zeros(length(xyz), 3, n_steps);
E_kin_all = zeros(n_steps,1);
E_pot_all = zeros(n_steps,1);
T_all = zeros(n_steps,1);

iteration = 0;
v_Betrag = zeros(length(xyz), 1);
T_0 = 50; % Zieltemperatur
k_B = 3.1651e-06; % Boltzmann-Konstante in a. u.
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
    F_neu = LJ_Kraft(xyz, sigma, E, d);
    E_pot = sum(LJ_Pot(xyz, sigma, E, a, d));

    % Temperaturkontrolle über Skalierung der Geschwindigkeiten mithilfe
    % des Skalierungsfaktors lambda
    v = v + bsxfun(@rdivide, (F_0+F_neu), 2*m)*delta_t;
    for i = 1:size(xyz,1)
        v_Betrag(i,1) = norm(v(i,:));
    end
    E_kin = sum(.5*m*(v_Betrag.^2));
    T = 2*E_kin./(3*np_Teilchen*k_B);
    lambda = sqrt(1+delta_t/tau*((T_0./T)-1));
    v = v.*lambda;

    xyz_all(:, :, iteration) = xyz(:, :);
    E_kin_all(iteration,1) = E_kin(:,:);
    E_pot_all(iteration,1) = E_pot(:,:);
    T_all(iteration,1) = T(:,:);

    F_0 = F_neu;
end
E_tot = E_pot_all + E_kin_all;
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
Atomanzahl = repmat({'Ar'},n,1);
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
