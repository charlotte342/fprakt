%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementierung der Linked Cell Methode für     %
% Lennard-Jones-Wechselwirkungen                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Anfangsbedingungen: 2D-System aus N_cell = N_cell,x * N_cell,y (r_cut = 2.5*sigma, sigma = 1) mit N Partikel auf einem Gitter im Abstand 1.5*sigma und v_0 = 0

sigma = 1;
np_Teilchen = 8*8; 
coordinates_0 = Anfangspositionen(sigma, sqrt(np_Teilchen)-1);

% Linked Particle List:
% Array mit N_cell, Aufruf der einzelnen Zellen mit N_cell(x,y)F_0(all(F_0 == 0, 2), :) = [];

% unterschiedliche Teilchen:
m = 1;
sigma = 1;
E = sigma;
T_0 = 50;
velocities_0 = zeros(np_Teilchen, 2);
forces_0 = zeros(np_Teilchen, 2);

% Zellen erstellen --> Array of linked lists
% Datenstruktur für jedes Teilchen N --> aus Unterklasse, die dann verlinkt
% werden in Zellen
C_0 = zeros(np_Teilchen, 2);
N_all = Molekueldynamik(coordinates_0, velocities_0, 4, 1, 1, 1, 0, 1e-3, 1e4, 1, T_0);

for i = 1:np_Teilchen
    PL(i) = dlnode(coordinates_0(i,:), velocities_0(i,:), forces_0(i,:));
    C_0(i,:) = floor(PL(i).coordinates./(2.5*sigma)); % Position der Teilchen zuordnen
end

C = C_0;
% für Zellindizierung später - sonst negativ oder 0
while sum(bsxfun(@ge, zeros(length(C), 2), C), 'all') > 0
    C = C + 1;
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
CellArray_linear = cell(length(D),1);

for i = 1:length(iD)
    CellArray_linear{iD(i),1} = [CellArray_linear{iD(i),1}, PL(i)];
end

% Particle List - head einführen: jede Zelle mit M+1 Teilchen (erstes
% Teilchen kopieren
% Zellen neu anordnen
CellArray = cell(max(D(:,1)), max(D(:,2)));
for i = 1:length(D)
    particle_head(D(i,1), D(i,2)) = CellArray_linear{i,1}(1,1).deepCopy;
    particle_head(D(i,1), D(i,2)).insertBefore(CellArray_linear{i,1}(1,1));
    CellArray{D(i,1), D(i,2)} = [particle_head(D(i,1),D(i,2)), CellArray_linear{i,1}]; % deepCopy damit unveränderlich
end

% Ausklamüsern

% Berechnung der Kräfte

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);

F_0 = [];

dij = [-1,-1; -1,0; -1,1; 0,-1; 0,0; 0,1; 1,-1; 1,0; 1,1]; % alle möglichen Variationen in 2D

x = 1:n_cells_x;
y = 1:n_cells_y;

[X,Y] = meshgrid(x,y);
cells = [X(:),Y(:)]; % alle Zellkombinationen um loops zu sparen

% Ausklamüsern

% Berechnung der Kräfte

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);

x = 1:n_cells_x;
y = 1:n_cells_y;

[X,Y] = meshgrid(x,y);
cells = [X(:),Y(:)]; % alle Zellkombinationen um loops zu sparen


F_0 = zeros(np_Teilchen,2);

dij = [-1,-1; -1,0; -1,1; 0,-1; 0,0; 0,1; 1,-1; 1,0; 1,1]; % alle möglichen Variationen in 2D
xyz = [];

for h = 1:np_Teilchen
    for l = 1:length(dij)
        if D(iD(h),1)+dij(l,1) > 0 && D(iD(h),2)+dij(l,2) > 0 && D(iD(h),1)+dij(l,1)<= n_cells_x && D(iD(h),2)+dij(l,2) <= n_cells_y
            while ~isempty(CellArray{D(iD(h),1)+dij(l,1),  D(iD(h),2)+dij(l,2)}(1,1).Next)
                xyz =[xyz; CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Next.coordinates];
                CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1) = CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Next;
            end
        end
    end
    n = size(xyz, 1);
    F_neu = zeros(n,2);
    r = zeros(n,2);
    for n = 1:n
        r(n,:) = xyz(n,:)-coordinates_0(h,:); % Abstand zu jedem Teilchen k in Zelle ij (aus CellArray)
        F_neu(n,:) = 24*E/sum(r(n,:).^2, 2) * sigma/((sum(r(n,:).^2, 2))^3) * (1-2*sigma^6/(sum(r(n,:).^2,2)^3)) * r(n,:);
        F = F_neu;
        F(isnan(F_neu)) = 0;
    end
    F_0(h,:) = sum(F, 1);
    xyz = [];
    if D(iD(h),1)+dij(l,1) > 0 && D(iD(h),2)+dij(l,2) > 0 && D(iD(h),1)+dij(l,1)<= n_cells_x && D(iD(h),2)+dij(l,2) <= n_cells_y
        while ~isempty(CellArray{D(iD(h),1)+dij(l,1),  D(iD(h),2)+dij(l,2)}(1,1).Prev)
            CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1) = CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Prev;
        end
    end
end

%% Velocity-Verlet
tic

coordinates = coordinates_0;
n_steps = 1;
t = 0;
v = zeros(np_Teilchen,2);
delta_t = 1e-3;
tau = delta_t*10;
a = 4;

E_kin_all = zeros(n_steps,1);
E_pot_all = zeros(n_steps,1);
T_all = zeros(n_steps,1);

iteration = 0;
v_Betrag = zeros(length(xyz), 1);
T_0 = 50; % Zieltemperatur
k_B = 3.1651e-06; % Boltzmann-Konstante in a. u.

xyz_all = zeros(np_Teilchen, 2, n_steps);
v_all = zeros(np_Teilchen, 2, n_steps);
F_all = zeros(np_Teilchen, 2, n_steps);
while t < delta_t*n_steps
    t = t + delta_t;
    iteration = iteration + 1;
    coordinates = coordinates + delta_t*(v + F_0.*delta_t*0.5/m); % neue Positionen --> einordnen in Zellen, dann neue Kraftberechnung
    
    C_0 = zeros(np_Teilchen, 2);
    for i = 1:np_Teilchen
        PL(i).coordinates = coordinates(i,:); % Aktualisierung
        PL(i).forces = F_0(i,:);
        C_0(i,:) = round(coordinates(i,:)/(2.5*sigma)); % Position der Teilchen zuordnen
    end
% 
    while sum(bsxfun(@ge, zeros(length(C_0), 2), C_0), 'all') > 0
        C_0 = C_0 + 1;
    end
    [D, iA, iD] = unique(C_0, 'rows');

    % Entstehung neuer Zellen
    for i = 1:np_Teilchen
        if bsxfun(@gt, C_0(i,:), size(particle_head))
            cells = [cells; C_0(i,:)];
            PL(i).removeNode
            particle_head(C_0(i,1), C_0(i,2)) = PL(i).deepCopy;
        end
    end

    for i = 1:np_Teilchen
        if C_0(i,:)/ C(i,:) ~= 1
            PL(i).removeNode;
            PL(i).insertAfter(CellArray{C_0(i,1),C_0(i,2)}(1,1)); % neue Verknüpfung
        end
    end

    F_old = F_0;
%     F_0 = LJ_Kraft(CellArray, sigma, E, np_Teilchen, D, iD, coordinates);


    for h = 1:np_Teilchen
        xyz = [];
        if D(iD(h),1)+dij(l,1) > 0 && D(iD(h),2)+dij(l,2) > 0 && D(iD(h),1)+dij(l,1)<= n_cells_x && D(iD(h),2)+dij(l,2) <= n_cells_y
            while ~isempty(CellArray{D(iD(h),1)+dij(l,1),  D(iD(h),2)+dij(l,2)}(1,1).Prev)
                CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1) = CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Prev; % Anfang der Liste, wenn CellArray-Indizierung nicht mehr stimmt
            end
        end
        for l = 1:length(dij)
            if D(iD(h),1)+dij(l,1) > 0 && D(iD(h),2)+dij(l,2) > 0 && D(iD(h),1)+dij(l,1)<= n_cells_x && D(iD(h),2)+dij(l,2) <= n_cells_y
                while ~isempty(CellArray{D(iD(h),1)+dij(l,1),  D(iD(h),2)+dij(l,2)}(1,1).Next)
                    xyz =[xyz; CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Next.coordinates];
                    CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1) = CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Next;
                end
            end
        end
        n = size(xyz, 1);
        F_neu = zeros(n,2);
        r = zeros(n,2);
        for n = 1:n
            r(n,:) = xyz(n,:) - coordinates_0(h,:); % Abstand zu jedem Teilchen k in Zelle ij (aus CellArray)
            F_neu(n,:) = 24*E/sum(r(n,:).^2, 2) * sigma/((sum(r(n,:).^2, 2))^3) * (1-2*sigma^6/(sum(r(n,:).^2,2)^3)) * r(n,:);
            F = F_neu;
            F(isnan(F_neu)) = 0;
        end
        F_0(h,:) = sum(F, 1);
    end

    for h= 1:np_Teilchen
        while ~isempty(CellArray{D(iD(h),1)+dij(l,1),  D(iD(h),2)+dij(l,2)}(1,1).Prev)
            CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1) = CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Prev;
        end
    end



    v = v + bsxfun(@rdivide, (F_0+F_old), 2*m)*delta_t;
    E_pot = sum(LJ_Pot(coordinates, sigma, E, a))*0.5; % Korrektur: *0.5, für i~=j

    % Temperaturkontrolle über Skalierung der Geschwindigkeiten mithilfe
    % des Skalierungsfaktors lambda
    for i = 1:size(coordinates,1)
        v_Betrag(i,1) = norm(v(i,:));
    end


    E_kin = sum(.5*m*(v_Betrag.^2))*0.5; % Korrektur: *0.5, da i~=j
    T = 2*E_kin./(3*np_Teilchen*k_B);
    lambda = sqrt(1+delta_t/tau*((T_0./T)-1));
    v = v.*lambda;

    for i = 1:np_Teilchen
        PL(i).velocities = v(i,:);
    end

    xyz_all(:, :, iteration) = coordinates(:, :);
    v_all(:,:,iteration) = v(:,:);
    F_all(:,:,iteration) = F_0(:,:);
    E_pot_all(iteration,:) = E_pot;
    E_kin_all(iteration,:) = E_kin;
    T_all(iteration,:) = T;
end


toc
%
figure;
plot(T_all, '-b', 'LineWidth', 2);
title('Temperaturverlauf');
xlabel('Zeit'); ylabel('Temperatur / K'); grid on;
figure;
plot(E_pot_all+E_kin_all, '-g', 'LineWidth', 2);
hold on
plot(E_kin_all, '-b','LineWidth',2);
hold on
plot(E_pot_all, '-r', 'LineWidth', 2);
title('Energie des Systems');
xlabel('Zeit'); ylabel('Energie'); grid on;
legend('Gesamtenergie', 'kinetische Energie', 'potentielle Energie');






%% Kraftberechnung
F_0 = LJ_Kraft(CellArray, sigma, E);

%% Zeitintegration
n_steps = 200;
t = 0;
v = zeros(np_Teilchen,2);
delta_t = 1e-2;
tau = delta_t*100;
a = 4;

[xyz_all, v_all] = Zeitintegration(CellArray, coordinates_0, t, delta_t, F_0, n_steps, v, m, sigma, E, np_Teilchen, a, tau);

%%
% Visualisierung(xyz_all, 'Film.avi')
Generate_xyz(xyz_all, '402.xyz')

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

function F_0 = LJ_Kraft(CellArray, sigma, E, np_Teilchen, D, iD, coordinates)
% Berechnung der Kräfte

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);

xyz = [];
F_0 = [];

dij = [-1,-1; -1,0; -1,1; 0,-1; 0,0; 0,1; 1,-1; 1,0; 1,1]; % alle möglichen Variationen in 2D
F = [];
for h = 1:np_Teilchen
    for l = 1:length(dij)
        if D(iD(h),1)+dij(l,1) > 0 && D(iD(h),2)+dij(l,2) > 0 && D(iD(h),1)+dij(l,1)<= n_cells_x && D(iD(h),2)+dij(l,2) <= n_cells_y
            while ~isempty(CellArray{D(iD(h),1)+dij(l,1),  D(iD(h),2)+dij(l,2)}(1,1).Next)
                xyz =[xyz; CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Next.coordinates];
                CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1) = CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Next;
            end
        end
    end
    n = size(xyz, 1);
    F_neu = zeros(n,2);
    r = zeros(n,2);
    for n = 1:n
        r(n,:) = xyz(n,:) - coordinates(h,:); % Abstand zu jedem Teilchen k in Zelle ij (aus CellArray)
        F_neu(n,:) = 24*E/sum(r(n,:).^2, 2) * sigma/((sum(r(n,:).^2, 2))^3) * (1-2*sigma^6/(sum(r(n,:).^2,2)^3)) * r(n,:);
        F = F_neu;
        F(isnan(F_neu)) = 0;
    end
    F_0(h,:) = sum(F, 1);
    xyz = [];
    % reset
    if D(iD(h),1)+dij(l,1) > 0 && D(iD(h),2)+dij(l,2) > 0 && D(iD(h),1)+dij(l,1)<= n_cells_x && D(iD(h),2)+dij(l,2) <= n_cells_y
        while ~isempty(CellArray{D(iD(h),1)+dij(l,1),  D(iD(h),2)+dij(l,2)}(1,1).Prev)
            CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1) = CellArray{D(iD(h),1) + dij(l,1) , D(iD(h),2) + dij(l,2)}(1,1).Prev;
        end
    end
end
end

function E_pot_total = LJ_Pot(xyz, sigma, E, a)
n = length(xyz);
E_pot = zeros(n,1,n);
r = zeros(n,2,n);
E_pot_total = zeros(n,1);

for i=1:n
    r(:,:,i) = bsxfun(@minus, xyz(i,:), xyz);
    % neu: Energien auch direkt aus r berechnet
    E_pot(:,:,i) = a*E*(((sigma.^6)./(sum(r(:,:,i).^2,2).^3))-((sigma.^12)./((sum(r(:,:,i).^2,2).^6))));
    E_korr = E_pot;
    E_korr(isnan(E_pot)) = 0; % Korrektur für i = j
end
E_pot_total(:,:) = sum(E_korr, 3);
end


function [xyz_all, v] = Zeitintegration(CellArray, xyz, t, delta_t, F_0, n_steps, v, m, sigma, E, np_Teilchen, a, tau)
% Velocity-Verlet
E_kin_all = zeros(n_steps,1);
E_pot_all = zeros(n_steps,1);
T_all = zeros(n_steps,1);

iteration = 0;
v_Betrag = zeros(length(xyz), 1);
T_0 = 50; % Zieltemperatur
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
    for np = 1:n_particles
        ctrp = ctrP +1;
        CellArray{i,j}(1,1).Next.coordinates = current_cell(ctrP,:); % Koordinaten der Zelle
        CellArray{i,j}(1,1) = CellArray{i,j}(1,1).Next;
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

        F_old = F_0;
        F_0 = LJ_Kraft(CellArray, sigma, E);
        v = v + bsxfun(@rdivide, (F_0+F_old), 2*m)*delta_t;
        E_pot = sum(LJ_Pot(xyz, sigma, E, a))*0.5; % Korrektur: *0.5, für i~=j

        % Temperaturkontrolle über Skalierung der Geschwindigkeiten mithilfe
        % des Skalierungsfaktors lambda
        v = v + bsxfun(@rdivide, (F_0+F_old), 2*m)*delta_t;
        for i = 1:size(xyz,1)
            v_Betrag(i,1) = norm(v(i,:));
        end

        E_kin = sum(.5*m*(v_Betrag.^2))*0.5; % Korrektur: *0.5, da i~=j
        T = 2*E_kin./(3*np_Teilchen*k_B);
        lambda = sqrt(1+delta_t/tau*((T_0./T)-1));
        v = v.*lambda;

        xyz_all(:, :, iteration) = xyz(:, :);
        v_all(:,:,iteration) = v(:,:);
        F_all(:,:,iteration) = F_0(:,:);
        E_pot_all(iteration,:) = E_pot;
        E_kin_all(iteration,:) = E_kin;
        T_all(iteration,:) = T;
    end
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
% z_positions = zeros(n,1); 


for i=1:size(M,3)
    fprintf(fileID, '%d\n', n);
    fprintf(fileID, 'Moleküldynamik_3D\n');
    for j = 1:n
        fprintf(fileID, '%s %f %f %f\n', Atomanzahl{j}, M(j,1,i), M(j,2,i), 0);
    end
end

fclose(fileID);
end
