%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementierung der Linked Cell Methode für     %
% Lennard-Jones-Wechselwirkungen in 3D            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Anfangsbedingungen: 3D-System, Positionen aus 3.3 importiert
% N_cell = N_cell,x * N_cell,y (r_cut = 2.5*sigma)

% Parameter importieren
[a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter('Parameter_302.txt');
n_steps = 1e5;

% Teilchen in Box
[coordinates_0, v, d, np_Teilchen] = Initialisierung_PBC(sigma*2.5*3, sigma*2.5*3, sigma*2.5*3, 10, sigma);

% Zellen erstellen --> Array of linked lists
% Zellen einteilen/definieren: Simulationsdomäne ist das Teilchengitter

% Datenstruktur für jedes Teilchen N --> aus Unterklasse, die dann verlinkt
% werden in Zellen
C = zeros(np_Teilchen, 3);
tau = delta_t*1e3;
T_0 = 50;
N_all = Molekueldynamik(coordinates_0, a, sigma, E, m, t, delta_t, n_steps, tau, T_0);
velocities_0 = zeros(np_Teilchen, 2);
forces_0 = zeros(np_Teilchen, 2);

for i = 1:np_Teilchen
    N(i) = Teilchen(coordinates_0(i,:), 4, sigma, E, m); % Teilcheneigenschaften zuordnen
    PL(i) = dlnode(coordinates_0(i,:), velocities_0(i,:), forces_0(i,:)); % jedes Teilchen: Particle List erstellen

    C(i,:) = floor(N(i).coordinates_0./(2.5*sigma)); % Position der Teilchen in Zelle zuordnen     
end

% für Zellindizierung später - sonst negativ oder 0
while sum(bsxfun(@ge, zeros(length(C), 3), C), 'all') > 0
    C = C + 1;
end


%%
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
CellArray_linear = cell(length(D),1);

for i = 1:length(iD)
    CellArray_linear{iD(i),1} = [CellArray_linear{iD(i),1}, PL(i)];
end

% Particle List - head einführen: jede Zelle mit M+1 Teilchen (erstes
% Teilchen kopieren
% Zellen neu anordnen

for i = 1:length(D)
    CellArray{D(i,1), D(i,2), D(i,3)} = [CellArray_linear{i,1}(1,1),CellArray_linear{i,:}];
end


% Kraftberechnung - z. B. LJ_Kraft oder Coulomb ...
% Linked-Cell: Kraft Fi auf N(i) in Zelle ic aus kc Zellen (Nachbarzellen)
% F = LJ_Kraft(xyz, sigma, E, d);
%%

% benachbarte Zellen finden:
% alle relevanten Positionen für dieses Teilchen

% [a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter('Parameter_302.txt');

tau = delta_t*1e3;
T_0 = 50;

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);
n_cells_z = size(CellArray, 3);
n_particle = zeros(n_cells_x, n_cells_y, n_cells_z);

xyz = [];
F_0 = [];

% alle Variationen (ab R2024b mit combinations möglich)
d_ijk = [-1,-1,-1; -1,-1,0; -1,-1,1; -1,0,-1; -1,0,0; -1,0,1; -1,1,-1; -1,1,0; -1,1,1;
    0,-1,-1; 0,-1,0; 0,-1,1; 0,0,-1; 0,0,0; 0,0,1; 0,1,-1; 0,1,0; 0,1,1;
    1,-1,-1; 1,-1,0; 1,-1,1; 1,0,-1; 1,0,0; 1,0,1; 1,1,-1; 1,1,0; 1,1,1;];

x = 1:n_cells_x;
y = 1:n_cells_y;
z = 1:n_cells_z;

[X,Y,Z] = meshgrid(x,y,z);
cells = [X(:),Y(:),Z(:)]; % alle Zellkombinationen um loops zu sparen
np_cell = zeros(length(cells), 1);

for i=1:length(cells)
    np_cell(i) = size(CellArray{cells(i,1),cells(i,2),cells(i,3)},2);
end

for i= 1:n_cells_x % alle möglichen Zellen
    for j = 1:n_cells_y
        for k = 1:n_cells_z
            % alle relevanten Positionen für dieses Teilchen
            n_particle(i,j,k) = size(CellArray{i,j,k}, 2); % Anzahl Teilchen pro Zelle {i,j}
            for l = 2: n_particle(i,j,k) % Schleife über alle Teilchen in Zelle
                current_cell(l-1,:,i,j,k) = CellArray{i,j,k}(1,l).coordinates;
            end
        end
    end
end
for i=1:n_cells_x
    for j=1:n_cells_y
        for k = 1:n_cells_z
            xyz = []; % clear xyz
            for m = 1:length(d_ijk)
                if i+d_ijk(m,1) > 0 && j+d_ijk(m,2) > 0 && k+d_ijk(m,3) > 0 && i+d_ijk(m,1) <= n_cells_x && j+d_ijk(m,2) <= n_cells_y && k+d_ijk(m,3) <= n_cells_z
                    xyz = [xyz; current_cell(:,:,i+d_ijk(m,1),j+d_ijk(m,2), k+d_ijk(m,3))]; % alle relevanten Koordinaten für Teilchen i in Zelle k
                end
            end
            for l = 2:n_particle(i,j,k)
                if n_particle(i,j,k) > 0 % nicht leere Zelle
                    xyz(all(xyz == 0, 2), :) = [];
                    n = size(xyz, 1);
                    F_neu = zeros(n,3);
                    r = zeros(n,3);
                    for n = 1:n
                        r(n,:) = xyz(n,:) - CellArray{i,j,k}(1,l).coordinates; % Abstand zu jedem Teilchen k in Zelle ij (aus CellArray)
                        F_neu(n,:) = 24*E/sum(r(n,:).^2, 2) * sigma/((sum(r(n,:).^2, 2))^3) * (1-2*sigma^6/(sum(r(n,:).^2,2)^3)) * r(n,:);
                        F = F_neu;
                        F(isnan(F_neu)) = 0;
                    end
                    F(l-1,:) = sum(F, 1);
                    F_total(l-1,:,i,j,k) = F(l-1,:);
                end
            end
            F_0 = [F_0; F_total(:,:,i,j,k)];
        end
    end
end

F_0(all(F_0 == 0, 2), :) = [];



% % Zeitintegration - z. B. Velocity-Verlet
% [xyz_all, v, E_kin_all, E_pot_all, E_tot, T_all] = Zeitintegration(xyz, t, delta_t, F, n_steps, v, m, sigma, E, d, np_Teilchen, a, delta_t*1e3);


%% Matlab code for the recursive movement through list :
for np = 1: Np
ctrP = ctrP +1;
% recursively until depth -> Np
CellArray_linear ( nc2 , nc1 ) . PL (1 ,1) . Next . Data . x (1 ,:) = X ( ctrP ,:) ;
CellArray_linear ( nc2 , nc1 ) . PL (1 ,1) = CellArray_linear ( nc2 , nc1 ) . PL (1 ,1) . Next ;
end
% reset ' pointer ' to initial :
for np = 1: Np
CellArray_linear ( nc2 , nc1 ) . PL (1 ,1) = CellArray_linear ( nc2 , nc1 ) . PL (1 ,1) . Prev ;
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


function F_total = LJ_Kraft(xyz, sigma, E) %d
% Berechnung der Kräfte
n = length(xyz);
F_neu = zeros(n,2,n);
r = zeros(n,2,n);
F_total = zeros(n,2);
r_betrag = zeros(n,1,n);

for i=1:n
    r(:,:,i) = bsxfun(@minus, xyz(i,:), xyz);
%     r(:,:,i) = r(:,:,i) - d.*round(r(:,:,i)./d);
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
