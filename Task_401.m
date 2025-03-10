%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementierung einer MD-Klasse %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% Initialisierung der Daten mit np Punkt-Teilchen in einer nd-D Box mit den
% Boxmaßen d_x,y,z = 30

% Parameter importieren
[a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter('Parameter_302.txt');
n_steps = 20;

% Teilchen in Box
[xyz, v, d] = Initialisierung_PBC(30, 30, 30, 10, sigma);

tau = 41.32*1e3;
T_0 = 50; % Zieltemperatur, 50 K
xyz_0 = Molekueldynamik(xyz, v, a, sigma, E, m, t, delta_t, n_steps, tau, T_0);

Zeitintegration = xyz_0.VelocityVerlet;


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


function [xyz, v, d] = Initialisierung_PBC(d_x, d_y, d_z, np_Teilchen, sigma)
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
    r(:,:,i) = r(:,:,i) - .5*d.* round(r(:,:,i)./d); % PBC: wenn Abstand zu Teilchen in nächster Box kürzer
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