%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementierung der Linked Cell Methode für     %
% Lennard-Jones-Wechselwirkungen                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Anfangsbedingungen: 2D-System aus N_cell = N_cell,x * N_cell,y (r_cut = 2.5*sigma, sigma = 1) mit N Partikel auf einem Gitter im Abstand 1.5*sigma und v_0 = 0

% [a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter('Datenimport.txt');
sigma = 1;
np_Teilchen = 5*5; 
coordinates_0 = Anfangspositionen(sigma, sqrt(np_Teilchen)-1);

% Linked Particle List:
% Array mit N_cell, Aufruf der einzelnen Zellen mit N_cell(x,y)
% unterschiedliche Teilchen:
m = repmat(1, np_Teilchen, 1);
sigma = repmat(1, np_Teilchen, 1);
E = repmat(1, np_Teilchen, 1);

% Zellen erstellen --> Array of linked lists
% Zellen einteilen/definieren: Simulationsdomäne ist das Teilchengitter
% (10*1.5*sigma x 10*1.5*sigma), cell structure via meshgrid

% [X,Y] = meshgrid(N_x,N_y);
% Cell_coordinates = [X(:), Y(:)];

% Datenstruktur für jedes Teilchen N --> aus Unterklasse, die dann verlinkt
% werden in Zellen
C = zeros(np_Teilchen, 2);
N_all = Molekueldynamik(coordinates_0, 4, 1, 1, 1, 0, 1e-3, 1e4, 1e3, 50);
for i = 1:np_Teilchen
    N(i) = Teilchen(coordinates_0(i,:), 4, sigma(i), E(i), m(i)); % Teilcheneigenschaften zuordnen
    PL(i) = dlnode(N(i).coordinates_0); % jedes Teilchen: Particle List erstellen

    C(i,:) = round(N(i).coordinates_0./(2.5*sigma(i))) + 1; % Position der Teilchen zuordnen
end

[D, iA, iD] = unique(C, 'rows');

for i = 1:np_Teilchen
    for j = i+1 : np_Teilchen
        if iD(j) == iD(i) % gleich wenn in einer Zelle
            PL(j).insertAfter(PL(i)); % dann: Verlinkung der Particle Lists
        end
    end

end

% Zellen: gleicher Wert in Vektor iD
n_cells = length(D); % Anzahl Zellen
CellArray = cell(length(D), 1);

for i = 1:length(iD)
    CellArray{iD(i),1} = [CellArray{iD(i),1}, PL(i)];
end
% Particle List - head einführen: jede Zelle mit M+1 Teilchen (erstes
% Teilchen kopieren, dann ab i = 2
for i = 1:length(D)
    CellArray{i,1} = [CellArray{i,1}(1,1),CellArray{i,1}];
end

% Kraftberechnung - z. B. LJ_Kraft oder Coulomb ...
F = LJ_Kraft(xyz, sigma, E, d);

% % Zeitintegration - z. B. Velocity-Verlet
% [xyz_all, v, E_kin_all, E_pot_all, E_tot, T_all] = Zeitintegration(xyz, t, delta_t, F, n_steps, v, m, sigma, E, d, np_Teilchen, a, delta_t*1e3);


%% Matlab code for the recursive movement through list :
for np = 1: Np
ctrP = ctrP +1;
% recursively until depth -> Np
CellArray ( nc2 , nc1 ) . PL (1 ,1) . Next . Data . x (1 ,:) = X ( ctrP ,:) ;
CellArray ( nc2 , nc1 ) . PL (1 ,1) = CellArray ( nc2 , nc1 ) . PL (1 ,1) . Next ;
end
% reset ' pointer ' to initial :
for np = 1: Np
CellArray ( nc2 , nc1 ) . PL (1 ,1) = CellArray ( nc2 , nc1 ) . PL (1 ,1) . Prev ;
end


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


function coordinates_0 = Anfangspositionen(sigma, np_Teilchen)
% Ausgangspositionen und Speicherzuordnung
x = 0 : 1.5*sigma : 1.5*sigma*np_Teilchen;
y = 0 : 1.5*sigma : 1.5*sigma*np_Teilchen;
[X,Y] = meshgrid(x,y);
coordinates_0 = [X(:), Y(:)];
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
    E_pot = sum(LJ_Pot(xyz, sigma, E, a, d))*0.5; % Korrektur: *0.5, da i~=j

    % Temperaturkontrolle über Skalierung der Geschwindigkeiten mithilfe
    % des Skalierungsfaktors lambda
    v = v + bsxfun(@rdivide, (F_0+F_neu), 2*m)*delta_t;
    for i = 1:size(xyz,1)
        v_Betrag(i,1) = norm(v(i,:));
    end
    E_kin = sum(.5*m*(v_Betrag.^2))*0.5; % Korrektur: *0.5, da i~=j
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
