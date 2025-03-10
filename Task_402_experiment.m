%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementierung der Linked Cell Methode für     %
% Lennard-Jones-Wechselwirkungen                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Anfangsbedingungen: 2D-System aus N_cell = N_cell,x * N_cell,y (r_cut = 2.5*sigma, sigma = 1) mit N Partikel auf einem Gitter im Abstand 1.5*sigma und v_0 = 0

% [a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter('Datenimport.txt');
sigma = 1;
np_Teilchen = 100*100; 
coordinates_0 = Anfangspositionen(sigma, sqrt(np_Teilchen)-1);

% Linked Particle List:
% Array mit N_cell, Aufruf der einzelnen Zellen mit N_cell(x,y)
% unterschiedliche Teilchen:
m = 1;
sigma = 1;
E = sigma;

% Zellen erstellen --> Array of linked lists
% Zellen einteilen/definieren: Simulationsdomäne ist das Teilchengitter
% (10*1.5*sigma x 10*1.5*sigma), cell structure via meshgrid

% [X,Y] = meshgrid(N_x,N_y);
% Cell_coordinates = [X(:), Y(:)];

% Datenstruktur für jedes Teilchen N --> aus Unterklasse, die dann verlinkt
% werden in Zellen
C = zeros(np_Teilchen, 2);
N_all = Molekueldynamik(coordinates_0, 4, 1, 1, 1, 0, 1e-3, 1e4, 1e3, 50);
velocities_0 = zeros(np_Teilchen, 2);
forces_0 = zeros(np_Teilchen, 2);
for i = 1:np_Teilchen
%     N(i) = Teilchen(coordinates_0(i,:), 4, sigma(i), E(i), m(i)); % Teilcheneigenschaften zuordnen
%     PL(i) = dlnode(N(i).coordinates_0); % jedes Teilchen: Particle List erstellen
    PL(i) = dlnode(coordinates_0(i,:), velocities_0(i,:), forces_0(i,:));

    C(i,:) = round(PL(i).coordinates./(2.5*sigma)) + 1; % Position der Teilchen zuordnen
end

[D, iA, iD] = unique(C, 'rows');

for i = 1:np_Teilchen
    for j = i+1 : np_Teilchen
        if iD(j) == iD(i) % gleich wenn in einer Zelle
            PL(j).insertAfter(PL(i)); % dann: Verlinkung der Particle Lists
%             r(j,:,i) = PL(i).Data - PL(j).Data; % Abstände Teilchen in einer Zelle
%             r(i,:,j) = -r(j,:,i);
        end
    end

end

% Zellen: gleicher Wert in Vektor iD
n_cells = length(D); % Anzahl Zellen
CellArray_linear = cell(length(D),1);

for i = 1:length(iD)
    CellArray_linear{iD(i),1} = [CellArray_linear{iD(i),1}, PL(i)];
end

CellArray_linear{1,1}(1,1).coordinates = [1e-5,1e-5]; % hilfe für später --> später beheben


% Particle List - head einführen: jede Zelle mit M+1 Teilchen (erstes
% Teilchen kopieren
% Zellen neu anordnen
% CellArray = cell(length(D));
for i = 1:length(D)
    CellArray{D(i,1), D(i,2)} = [CellArray_linear{i,1}(1,1),CellArray_linear{i,1}];
end


% Kraftberechnung - z. B. LJ_Kraft oder Coulomb ...
% Linked-Cell: Kraft Fi auf N(i) in Zelle ic aus kc Zellen (Nachbarzellen)
% F = LJ_Kraft(xyz, sigma, E, d);

% benachbarte Zellen finden:
% alle relevanten Positionen für dieses Teilchen

% [a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter('Parameter_302.txt');

% tau = delta_t*1e3;
% T_0 = 50;
 

%% Kraftberechnung

F_0 = LJ_Kraft(CellArray, sigma, E);

%% Zeitintegration

n_steps = 100;
t = 0;
v = zeros(np_Teilchen,2);
delta_t = 1e-3;
T_0 = 50; % Zieltemperatur
tau = delta_t*1e3;

[xyz_all, v_all, F_all] = Zeitintegration(coordinates_0, t, delta_t, F_0, n_steps, v, m, sigma, E, np_Teilchen, tau, CellArray);


% xyz = coordinates_0;
% E_kin_all = zeros(n_steps,1);
% E_pot_all = zeros(n_steps,1);
% T_all = zeros(n_steps,1);
% 
% iteration = 0;
% v_Betrag = zeros(length(xyz), 1);
% k_B = 3.1651e-06; % Boltzmann-Konstante in a. u.
% 
% xyz_all = zeros(np_Teilchen, 2, n_steps);
% v_all = zeros(np_Teilchen, 2, n_steps);
% F_all = zeros(np_Teilchen, 2, n_steps);
% while t < delta_t*n_steps
%     t = t + delta_t;
%     iteration = iteration + 1;
%     xyz = xyz + delta_t*(v + F_0.*delta_t*0.5/m); % neue Positionen --> einordnen in Zellen, dann neue Kraftberechnung
%     for i = 1:np_Teilchen
%         PL(i).coordinates = xyz(i,:); % Aktualisierung
%         PL(i).forces = F_0(i,:);
%     end
%     F_old = F_0; 
%     F_0 = LJ_Kraft(CellArray, sigma, E);
%     v = v + bsxfun(@rdivide, (F_0+F_old), 2*m)*delta_t;
% 
% 
%     xyz_all(:, :, iteration) = xyz(:, :);
%     v_all(:,:,iteration) = v(:,:);
%     F_all(:,:,iteration) = F_0(:,:);
% end

%%
% Visualisierung(xyz_all, 'Film.avi')
xyz_z = zeros(np_Teilchen, 1, n_steps);
xyz_all = [xyz_all, xyz_z];
Generate_xyz(xyz_all, '402_experiment.xyz')



%% Matlab code for the recursive movement through list :
for np = 1: np_Teilchen
ctrP = ctrP +1;
% recursively until depth -> Np
CellArray( nc2 , nc1 ) . PL (1 ,1) . Next . Data . x (1 ,:) = X ( ctrP ,:) ;
CellArray( nc2 , nc1 ) . PL (1 ,1) = CellArray( nc2 , nc1 ) . PL (1 ,1) . Next ;
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


function coordinates_0 = Anfangspositionen(sigma, np_Teilchen)
% Ausgangspositionen und Speicherzuordnung
x = 0 : 1.5*sigma : 1.5*sigma*np_Teilchen;
y = 0 : 1.5*sigma : 1.5*sigma*np_Teilchen;
[X,Y] = meshgrid(x,y);
coordinates_0 = [X(:), Y(:)];
end

function F_0 = LJ_Kraft(CellArray, sigma, E) %d
% Berechnung der Kräfte

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);

n_particle = zeros(n_cells_x, n_cells_y); % Teilchen pro Zelle
xyz = [];
F_0 = [];

for i= 1:n_cells_x % alle möglichen Zellen
    for j = 1:n_cells_y
        % alle relevanten Positionen für dieses Teilchen
        n_particle(i,j) = size(CellArray{i,j}, 2); % Anzahl Teilchen pro Zelle {i,j}
        for k = 2: n_particle(i,j) % Schleife über alle Teilchen in Zelle
            current_cell(k-1,:,i,j) = CellArray{i,j}(1,k).coordinates;
        end
    end
end
for i=1:n_cells_x
    for j=1:n_cells_y
        for k = 2:n_particle(i,j)
            xyz = []; % clear xyz
             for l = -1:1 %  nur benachbarte Zellen, in 2D: 9 Zellen
                for m = -1:1
                    if i+l > 0 && j+m > 0 && i+l <= n_cells_x && j+m <= n_cells_y
                        xyz = [xyz; current_cell(:,:,i+l,j+m)]; % alle relevanten Koordinaten für Teilchen i in Zelle k 
                    end
                end
             end
            xyz(all(xyz == 0, 2), :) = [];
            n = size(xyz, 1);
            F_neu = zeros(n,2);
            r = zeros(n,2);
            for n = 1:n 
                r(n,:) = xyz(n,:) - CellArray{i,j}(1,k).coordinates; % Abstand zu jedem Teilchen k in Zelle ij (aus CellArray)
                F_neu(n,:) = 24*E/sum(r(n,:).^2, 2) * sigma/((sum(r(n,:).^2, 2))^3) * (1-2*sigma^6/(sum(r(n,:).^2,2)^3)) * r(n,:);
                F = F_neu;
                F(isnan(F_neu)) = 0;
            end
            F(k-1,:) = sum(F, 1);
            F_total(k-1,:,i,j) = F(k-1,:);
        end
        F_0 = [F_0; F_total(:,:,i,j)];
    end
end

F_0(all(F_0 == 0, 2), :) = [];
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


function [xyz_all, v, F_all] = Zeitintegration(xyz, t, delta_t, F_0, n_steps, v, m, sigma, E,np_Teilchen, tau, CellArray)
% Velocity-Verlet

iteration = 0;
k_B = 3.1651e-06; % Boltzmann-Konstante in a. u.

xyz_all = zeros(np_Teilchen, 2, n_steps);
v_all = zeros(np_Teilchen, 2, n_steps);
F_all = zeros(np_Teilchen, 2, n_steps);
while t < delta_t*n_steps
    t = t + delta_t;
    iteration = iteration + 1;
    xyz = xyz + delta_t*(v + F_0.*delta_t*0.5/m); % neue Positionen --> einordnen in Zellen, dann neue Kraftberechnung
    for i = 1:np_Teilchen
        PL(i).coordinates = xyz(i,:); % Aktualisierung
        PL(i).forces = F_0(i,:);
    end
    F_old = F_0; 
    F_0 = LJ_Kraft(CellArray, sigma, E);
    v = v + bsxfun(@rdivide, (F_0+F_old), 2*m)*delta_t;


    xyz_all(:, :, iteration) = xyz(:, :);
    v_all(:,:,iteration) = v(:,:);
    F_all(:,:,iteration) = F_0(:,:);
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
    plot(xyz_all(:, 1, i), xyz_all(:, 2, i), 'o', 'MarkerSize', 6);
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
        fprintf(fileID, '%s %f %f %f\n', Atomanzahl{j}, M(j,1,i), M(j,2,i), M(j,3,i));
    end
end

fclose(fileID);
end

