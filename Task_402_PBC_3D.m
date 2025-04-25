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

% Linked Particle List:
% Array mit N_cell, Aufruf der einzelnen Zellen mit N_cell(x,y)
% unterschiedliche Teilchen:

T_0 = 50;
tau = delta_t*1e3;
velocities_0 = zeros(np_Teilchen, 3);
forces_0 = zeros(np_Teilchen, 3);

% Zellen erstellen --> Array of linked lists
% Datenstruktur für jedes Teilchen N --> aus Unterklasse, die dann verlinkt
% werden in Zellen

np_Teilchen = size(coordinates_0,1);
C_0 = zeros(np_Teilchen, 3);
N_all = Molekueldynamik(coordinates_0, velocities_0, a, sigma, E, m, t, delta_t, n_steps, tau, T_0);

for i = 1:np_Teilchen
    PL(i) = dlnode(coordinates_0(i,:), velocities_0(i,:), forces_0(i,:));
    C_0(i,:) = floor(PL(i).coordinates./(2.5*sigma)); % Position der Teilchen zuordnen
end

C = C_0;
% für Zellindizierung später - sonst negativ oder 0
while sum(bsxfun(@ge, zeros(length(C), 3), C), 'all') > 0
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
    particle_head(D(i,1), D(i,2), D(i,3)) = CellArray_linear{i,1}(1,1).deepCopy;
    particle_head(D(i,1), D(i,2), D(i,3)).insertBefore(CellArray_linear{i,1}(1,1));
    CellArray{D(i,1), D(i,2), D(i,3)} = [particle_head(D(i,1),D(i,2),D(i,3)), CellArray_linear{i,1}]; % deepCopy damit unveränderlich
end

% Ausklamüsern

%% Berechnung der Kräfte

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);
n_cells_z = size(CellArray, 3);

x = -1:1;
y = -1:1;
z = -1:1;

[X,Y,Z] = meshgrid(x,y,z);
dxyz = [X(:), Y(:), Z(:)];


x = 1:n_cells_x;
y = 1:n_cells_y;
z = 1:n_cells_z;

[X,Y,Z] = meshgrid(x,y,z);
cells = [X(:),Y(:),Z(:)]; % alle Zellkombinationen um loops zu sparen


for i = 1:length(cells)
    k=1;
    if ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)})
        while ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next)
            currcell(k,:,cells(i,1),cells(i,2),cells(i,3)) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next.coordinates;
            CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next;
            k=k+1;
        end
    end
end
% reset ' pointer ' to initial :
for i = 1:length(cells)
    if ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)})
        while ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Prev)
            CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Prev;
        end
    end
end


coordinates = [];
F_0 = [];
for i= 1:length(cells) % alle möglichen Zellen
    k = 1;
    F = [];
    if ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)})
        while ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next) % zum Zählen der Elemente pro linked list (pro Zelle, ungleich CellArray)
            CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next;
            coordinates = [coordinates; currcell(k,:,cells(i,1),cells(i,2),cells(i,3))];
            xyz = []; % clear xyz
            for l = 1:length(dxyz)
                % PBC!!
                for j = 1:3
                    if cells(i,j) + dxyz(l,j) <= 0
                        dxyz(l,j) = dxyz(l,j) + size(CellArray,j);
                    elseif cells(i,j) + dxyz(l,j) > size(CellArray,j)
                        dxyz(l,j) = dxyz(l,j)-size(CellArray,j);
                    end
                end
                xyz = [xyz; currcell(:,:,cells(i,1)+dxyz(l,1),cells(i,2)+dxyz(l,2), cells(i,3)+dxyz(l,3))]; % alle relevanten Koordinaten für Teilchen k in Zelle i
            end
            xyz(all(xyz == 0, 2),:) = [];
            if size(xyz,1) == 1
                F_total(k,:,cells(i,1),cells(i,2),cells(i,3)) = zeros(1,3);
                break
            end
            n=size(xyz,1);
            F_neu = zeros(n,3);
            r = zeros(n,3);
            for n = 1:size(xyz,1)
                r(n,:) = xyz(n,:) - currcell(k,:,cells(i,1),cells(i,2), cells(i,3)); % Abstand zu jedem Teilchen k in Zelle ij (aus CellArray)
                r(n,:) = r(n,:) - d.* round(r(n,:)./d); % PBC: wenn Abstand zu Teilchen in nächster Box kürzer

                if norm(r(n,:)) > 2.5*sigma
                    r(n,:) = zeros(1,3);
                end
                F_neu(n,:) = 24*E/sum(r(n,:).^2, 2) * sigma/((sum(r(n,:).^2, 2))^3) * (1-2*sigma^6/(sum(r(n,:).^2,2)^3)) * r(n,:);
                F_neu(isnan(F_neu)) = 0;
            end
            F(k,:) = sum(F_neu, 1);
            F_total(k,:,cells(i,1),cells(i,2),cells(i,3)) = F(k,:);
            k = k+1;
        end
        F_0 = [F_0; F_total(:,:,cells(i,1),cells(i,2),cells(i,3))];
    end
end

% F_0(all(F_0 == 0, 3), :) = [];

% reset ' pointer ' to initial :
for i = 1:length(cells)
    if ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)})
        while ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Prev)
            CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Prev;
        end
    end
end





%% Velocity-Verlet
tic


n_steps = 1;
t = 0;
% for i = 1:np_Teilchen
%     v(i,:) = PL(i).velocities;
% end
v = zeros(np_Teilchen,3);
delta_t = 1e-2;
tau = delta_t*1e3;
a = 4;

E_kin_all = zeros(n_steps,1);
E_pot_all = zeros(n_steps,1);
T_all = zeros(n_steps,1);

iteration = 0;
v_Betrag = zeros(length(coordinates), 1);
T_0 = 50; % Zieltemperatur
k_B = 3.1651e-06; % Boltzmann-Konstante in a. u.

xyz_all = zeros(np_Teilchen, 3, n_steps);
v_all = zeros(np_Teilchen, 3, n_steps);
F_all = zeros(np_Teilchen, 3, n_steps);
while t < delta_t*n_steps
    t = t + delta_t;
    iteration = iteration + 1;
    coordinates = coordinates + delta_t*(v + F_0.*delta_t*0.5/m); % neue Positionen --> einordnen in Zellen, dann neue Kraftberechnung
    for i = 1:size(coordinates, 1)
        for j = 1:size(coordinates, 2)
            if coordinates(i,j) >= d(j)
                coordinates(i,j) = coordinates(i,j) - d(j);
            elseif coordinates(i,j) < 0
                coordinates(i,j) = coordinates(i,j) + d(j);
            end
        end
    end

    C_0 = zeros(np_Teilchen, 3);
    for i = 1:np_Teilchen
        PL(i).coordinates = coordinates(i,:); % Aktualisierung
        PL(i).forces = F_0(i,:);
        C_0(i,:) = ceil(PL(i).coordinates./(2.5*sigma)); % Position der Teilchen zuordnen
    end
 
    while sum(bsxfun(@ge, zeros(length(C_0), 2), C_0), 'all') > 0
        C_0 = C_0 + 1;
    end

    for i = 1:np_Teilchen
        if C_0(i,:)/ C(i,:) ~= 1
            PL(i).removeNode;
            PL(i).insertAfter(particle_head(C_0(i,1), C_0(i,2), C_0(i,3))); % neue Verknüpfung
        end
    end
    F_0 = LJ_Kraft(CellArray, sigma, E, d);


    % reset ' pointer ' to initial :
    for i = 1:length(cells)
        while ~isempty(particle_head(cells(i,1),cells(i,2)).Prev)
            particle_head(cells(i,1),cells(i,2)) = particle_head(cells(i,1),cells(i,2)).Prev;
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


%% Kraftberechnung
F_ok = LJ_Kraft(CellArray, sigma, E, d);

%% Zeitintegration
n_steps = 100;
t = 0;
v = zeros(np_Teilchen,2);
delta_t = 1e-2;
tau = delta_t*1e3;
a = 4;

[xyz_all, v_all] = Zeitintegration(CellArray, coordinates_0, t, delta_t, F_0, n_steps, v, m, sigma, E, np_Teilchen, a, tau);

%%
Visualisierung(xyz_all, 'Film.avi')

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
            r(:,:,i) = r(:,:,i) - d.* round(r(:,:,i)./d);
            r_betrag(j,1,i) = norm(r(j,:,i));
            r_betrag(i,1,j) = r_betrag(j,1,i);
        end
    end
end     
for i = 1:n
    for j = i:n
        while r_betrag(j,1,i) < sigma
            xyz_0(j,:) = rand(1,3);
            xyz(j,:) = bsxfun(@minus, bsxfun(@times, xyz_0(j,:), d), 0.5*d);
            r(:,:,i) = bsxfun(@minus, xyz(i,:), xyz);
            r(:,:,i) = r(:,:,i) - d.* round(r(:,:,i)./d);
            r_betrag(j,1,i) = norm(r(j,:,i));
        end
    end
end
v = zeros(n,3);
end

function F_0 = LJ_Kraft(CellArray, sigma, E, d) %d
% Berechnung der Kräfte

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);
n_cells_z = size(CellArray, 3);

x = -1:1;
y = -1:1;
z = -1:1;

[X,Y,Z] = meshgrid(x,y,z);
dxyz = [X(:), Y(:), Z(:)];


x = 1:n_cells_x;
y = 1:n_cells_y;
z = 1:n_cells_z;

[X,Y,Z] = meshgrid(x,y,z);
cells = [X(:),Y(:),Z(:)]; % alle Zellkombinationen um loops zu sparen

currcell = [];
for i = 1:length(cells)
    k=1;
    if ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)})
        while ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next)
            currcell(k,:,cells(i,1),cells(i,2),cells(i,3)) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next.coordinates;
            CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next;
            k=k+1;
        end
    end
end
% reset ' pointer ' to initial :
for i = 1:length(cells)
    if ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)})
        while ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Prev)
            CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Prev;
        end
    end
end


coordinates = [];
F_0 = [];
for i= 1:length(cells) % alle möglichen Zellen
    k = 1;
    F = [];
    if ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)})
        while ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next) % zum Zählen der Elemente pro linked list (pro Zelle, ungleich CellArray)
            CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Next;
            coordinates = [coordinates; currcell(k,:,cells(i,1),cells(i,2),cells(i,3))];
            xyz = []; % clear xyz
            for l = 1:length(dxyz)
                % PBC!!
                for j = 1:3
                    if cells(i,j) + dxyz(l,j) <= 0
                        dxyz(l,j) = dxyz(l,j) + size(CellArray,j);
                    elseif cells(i,j) + dxyz(l,j) > size(CellArray,j)
                        dxyz(l,j) = dxyz(l,j)-size(CellArray,j);
                    end
                end
                xyz = [xyz; currcell(:,:,cells(i,1)+dxyz(l,1),cells(i,2)+dxyz(l,2), cells(i,3)+dxyz(l,3))]; % alle relevanten Koordinaten für Teilchen k in Zelle i
            end
            xyz(all(xyz == 0, 2),:) = [];
            if size(xyz,1) == 1
                F_total(k,:,cells(i,1),cells(i,2),cells(i,3)) = zeros(1,3);
                break
            end
            n=size(xyz,1);
            F_neu = zeros(n,3);
            r = zeros(n,3);
            for n = 1:size(xyz,1)
                r(n,:) = xyz(n,:) - currcell(k,:,cells(i,1),cells(i,2), cells(i,3)); % Abstand zu jedem Teilchen k in Zelle ij (aus CellArray)
                r(n,:) = r(n,:) - d.* round(r(n,:)./d); % PBC: wenn Abstand zu Teilchen in nächster Box kürzer

                if norm(r(n,:)) > 2.5*sigma
                    r(n,:) = zeros(1,3);
                end
                F_neu(n,:) = 24*E/sum(r(n,:).^2, 2) * sigma/((sum(r(n,:).^2, 2))^3) * (1-2*sigma^6/(sum(r(n,:).^2,2)^3)) * r(n,:);
                F_neu(isnan(F_neu)) = 0;
            end
            F(k,:) = sum(F_neu, 1);
            F_total(k,:,cells(i,1),cells(i,2),cells(i,3)) = F(k,:);
            k = k+1;
        end
        F_0 = [F_0; F_total(:,:,cells(i,1),cells(i,2),cells(i,3))];
    end
end

% reset ' pointer ' to initial :
for i = 1:length(cells)
    if ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)})
        while ~isempty(CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Prev)
            CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1) = CellArray{cells(i,1),cells(i,2),cells(i,3)}(1,1).Prev;
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

