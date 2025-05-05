%%%%%%%%%%
% Wasser - Ewald Summe %
%%%%%%%%%%
%Real-space force: use your short-range Coulomb code (as above) with cutoff rcrc​
%Assign charges to a 3D grid (e.g., using B-splines or Gaussians)
% 3D FFT the grid to get ρ(k)ρ(k)
% Multiply by the Ewald Green’s function in Fourier space
% Inverse FFT to get the potential field
% Interpolate forces back to particle positions 
% FFTN for the 3D Fourier transforms
% ndgrid or meshgrid for reciprocal vectors
% erfc, exp, erfcinv for accurate real/reciprocal splitting



%% Anfangsbedingungen: 3D-System, Positionen aus 3.3 importiert
% N_cell = N_cell,x * N_cell,y (r_cut = 2.5*sigma)

% Parameter importieren, https://docs.lammps.org/Howto_tip3p.html
mass_O = 15.9994; % amu
mass_H = 1.008;
q_O = -0.834 * 1.602e-19; % e
q_H = 0.417 * 1.602e-19;

N_A = 6.022e23;
epsilon_OO = 0.1521 * 1/N_A * 4.184; % kcal/mole to kJ
sigma_OO = 3.1507; % Angström
epsilon_HH = 0.0 * 1/N_A * 4.184;
sigma_HH = 1.0;
epsilon_OH = 0.0 * 1/N_A * 4.184;
sigma_OH = 1.0;
K_OH = 450; % kcal/mole/A², bond, Federkonstante
r0_OH = 0.9572; % A
K_HOH = 55.0; % kcal/mole, angle
theta0_HOH = 104.52;

epsilon_0 =8.86e-22*1e3; % C²/(N*Angstrom²)
skal = 1/(4*pi*epsilon_0);
G = 0.35; % Parameter der Gaußschen Glockenkurve / A⁻¹
a = 4;

Atome_dim = 3;

% Teilchen in Box
x = 0 : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim;
y = 0 : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim;
z = 0 : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim;
[X,Y,Z] = meshgrid(x,y,z);

coordinates_O = [X(:), Y(:), Z(:)];
np_oxygen = length(coordinates_O);
coordinates_O = [repmat(q_O, np_oxygen, 1), coordinates_O];

u = 0 - r0_OH : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim;
v = 0 : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim;
w = 0 : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim;
[U,V,W] = meshgrid(u,v,w);

q = 0 - cos(theta0_HOH - 90)*r0_OH : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim - cos(theta0_HOH - 90)*r0_OH;
r = 0 + sin(theta0_HOH - 90)*r0_OH : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim + sin(theta0_HOH - 90)*r0_OH;
s = 0 : 1.5*sigma_OO : 1.5*sigma_OO*Atome_dim;
[Q,R,S] = meshgrid(q,r,s);

coordinates_H = [U(:), V(:), W(:); Q(:), R(:), S(:)];
np_hydrogen = length(coordinates_H);
% velocities_0 = zeros(3,np_Teilchen);

coordinates_0 = [coordinates_O; repmat(q_H, np_hydrogen, 1), coordinates_H];
np_Teilchen = length(coordinates_0);
%%
Generate_xyz(coordinates_0, '501.xyz', q_O, q_H)

%% Linked Particle List:
% Array mit N_cell, Aufruf der einzelnen Zellen mit N_cell(x,y)
% unterschiedliche Teilchen:

T_0 = 50;
delta_t = 1e-3;
tau = delta_t*1e3;
velocities_0 = zeros(np_Teilchen, 3);
forces_0 = zeros(np_Teilchen, 3);

% Zellen erstellen --> Array of linked lists
% Datenstruktur für jedes Teilchen N --> aus Unterklasse, die dann verlinkt
% werden in Zellen

np_Teilchen = size(coordinates_0,1);
C_0 = zeros(np_Teilchen, 3);
% N_all = Molekueldynamik(coordinates_0, velocities_0, a, sigma, E, m, t,
% delta_t, n_steps, tau, T_0); später alles als Klasse

for i = 1:np_Teilchen
    PL(i) = dlnode(coordinates_0(i,:), velocities_0(i,:), forces_0(i,:));
    C_0(i,:) = floor(PL(i).coordinates(2:4)./(2.5*sigma_OO)); % Position der Teilchen zuordnen
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

coordinates = coordinates_0;

%% Energieberechnung

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);
n_cells_z = size(CellArray, 3);

x = -1:1;
y = -1:1;
z = -1:1;

[X,Y,Z] = meshgrid(x,y,z);
E_kr = zeros(np_Teilchen, 1);
F_kurz = zeros(np_Teilchen, 3);

for h = 1:np_Teilchen
    xyz = [];
    dxyz = [X(:), Y(:), Z(:)];
    for l = 1:length(dxyz)
        if D(iD(h),1)+dxyz(l,1) < 1
            dxyz(l,1) = size(CellArray,1) - D(iD(h),1);
        end
        if D(iD(h),2)+dxyz(l,2) < 1
            dxyz(l,2) = size(CellArray,2) - D(iD(h),2);
        end
        if D(iD(h),3)+dxyz(l,3) < 1
            dxyz(l,3) = size(CellArray,3) - D(iD(h),3);
        end
        if D(iD(h),1)+dxyz(l,1) > n_cells_x
            dxyz(l,1) = 1 - D(iD(h),1);
        end
        if D(iD(h),2)+dxyz(l,2) > n_cells_y
            dxyz(l,2) = 1 - D(iD(h),2);
        end
        if D(iD(h),3)+dxyz(l,3) > n_cells_z
            dxyz(l,3) = 1 - D(iD(h),2);
        end
        if ~isempty(CellArray{D(iD(h),1)+dxyz(l,1),  D(iD(h),2)+dxyz(l,2), D(iD(h),3)+dxyz(l,3)}) % Zelle gefüllt oder nicht
            while ~isempty(CellArray{D(iD(h),1),  D(iD(h),2), D(iD(h),3)}(1,1).Prev) % zurücksetzen
                CellArray{D(iD(h),1), D(iD(h),2), D(iD(h),3)}(1,1) = CellArray{D(iD(h),1), D(iD(h),2), D(iD(h),3)}(1,1).Prev;
            end
            while ~isempty(CellArray{D(iD(h),1)+dxyz(l,1),  D(iD(h),2)+dxyz(l,2), D(iD(h),3)+dxyz(l,3)}(1,1).Next)
                xyz =[xyz; CellArray{D(iD(h),1) + dxyz(l,1) , D(iD(h),2) + dxyz(l,2), D(iD(h),3) + dxyz(l,3)}(1,1).Next.coordinates];
                CellArray{D(iD(h),1) + dxyz(l,1) , D(iD(h),2) + dxyz(l,2), D(iD(h),3) + dxyz(l,3)}(1,1) = CellArray{D(iD(h),1) + dxyz(l,1) , D(iD(h),2) + dxyz(l,2), D(iD(h),3) + dxyz(l,3)}(1,1).Next;
            end
        end
    end
    n = size(xyz, 1);
    r = zeros(n,4);
    r_betrag = zeros(n,1);
    E_neu = zeros(n,1);
    for i = 1:n
        r(i,:) = xyz(i,:)-xyz(h,:); % Abstand zu jedem Teilchen k in Zelle aus CellArray
        r_betrag(i) = norm(r(i,2:4));
    end
    for i = 1:n
        if r_betrag(i) > 2.5*sigma_OO
            E_neu(i) = 0;
            continue
        end
        if r_betrag(i) == 0
            E_neu(i) = 0;
        elseif xyz(h,1) == q_O && r(i,1) == 0
            E_neu(i) = a*epsilon_OO*(((sigma_OO.^6)./(sum(r(i,2:4).^2,2).^3))-((sigma_OO.^12)./((sum(r(i, 2:4).^2,2).^6)))) + .5*skal*q_O^2*erfc(G*r_betrag(i))/r_betrag(i);
        elseif xyz(h,1) == q_O && r(i,1) ~= 0 || xyz(h,1) == q_H && r(i,1) ~= 0
            E_neu(i) = a*epsilon_OH*(((sigma_OH.^6)./(sum(r(i,2:4).^2,2).^3))-((sigma_OH.^12)./((sum(r(i, 2:4).^2,2).^6)))) + .5*skal*q_O*q_H*erfc(G*r_betrag(i))/r_betrag(i);
        elseif xyz(h,1) == q_H && r(i,1) == 0
            E_neu(i) = a*epsilon_HH*(((sigma_HH.^6)./(sum(r(i,2:4).^2,2).^3))-((sigma_HH.^12)./((sum(r(i, 2:4).^2,2).^6)))) + .5*skal*q_H^2*erfc(G*r_betrag(i))/r_betrag(i);
        end
    end
    E_kr(h,:) = sum(E_neu, 1);
    for i = 1:np_Teilchen
        while ~isempty(CellArray{D(iD(i),1),  D(iD(i),2), D(iD(i),3)}(1,1).Prev)
            CellArray{D(iD(i),1), D(iD(i),2), D(iD(i),3)}(1,1) = CellArray{D(iD(i),1), D(iD(i),2), D(iD(i),3)}(1,1).Prev;
        end
    end
end
E_pot = .5*sum(E_kr,1);

%% Berechnung der Kräfte

F_LJ = LJ_Kraft(CellArray, sigma_HH, sigma_OH, sigma_OO, epsilon_HH, epsilon_OH, epsilon_OO, np_Teilchen, D, iD, coordinates, q_O, q_H);
F_Coulomb = Coulomb_Kraft(coordinates, q_O, q_H, np_Teilchen, sigma_OO, Atome_dim);
F_0 = F_LJ + F_Coulomb;

%% Velocity-Verlet
tic

n_steps = 10;
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
 
    while sum(bsxfun(@ge, zeros(length(C_0), 3), C_0), 'all') > 0
        C_0 = C_0 + 1;
    end

    for i = 1:np_Teilchen
        if C_0(i,:)/ C(i,:) ~= 1
            PL(i).removeNode;
            PL(i).insertAfter(particle_head(C_0(i,1), C_0(i,2), C_0(i,3))); % neue Verknüpfung
        end
    end

    F_old = F_0;
    F_LJ = LJ_Kraft(CellArray, sigma_HH, sigma_OH, sigma_OO, epsilon_HH, epsilon_OH, epsilon_OO, np_Teilchen, D, iD, coordinates, q_O, q_H);
    F_Coulomb = Coulomb_Kraft(coordinates, q_O, q_H, np_Teilchen, sigma_OO, Atome_dim);
    F_
    F_0 = F_LJ + F_Coulomb;

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

function F_0 = LJ_Kraft(CellArray, sigma_HH, sigma_OH, sigma_OO, epsilon_HH, epsilon_OH, epsilon_OO, np_Teilchen, D, iD, coordinates, q_O, q_H) %d
% Berechnung der Kräfte

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);
n_cells_z = size(CellArray, 3);

x = -1:1;
y = -1:1;
z = -1:1;

[X,Y,Z] = meshgrid(x,y,z);

for h = 1:np_Teilchen
    xyz = [];
    dxyz = [X(:), Y(:), Z(:)];
    for l = 1:length(dxyz)
        if D(iD(h),1)+dxyz(l,1) < 1
            dxyz(l,1) = size(CellArray,1) - D(iD(h),1);
        end
        if D(iD(h),2)+dxyz(l,2) < 1
            dxyz(l,2) = size(CellArray,2) - D(iD(h),2);
        end
        if D(iD(h),3)+dxyz(l,3) < 1
            dxyz(l,3) = size(CellArray,3) - D(iD(h),3);
        end
        if D(iD(h),1)+dxyz(l,1) > n_cells_x
            dxyz(l,1) = 1 - D(iD(h),1);
        end
        if D(iD(h),2)+dxyz(l,2) > n_cells_y
            dxyz(l,2) = 1 - D(iD(h),2);
        end
        if D(iD(h),3)+dxyz(l,3) > n_cells_z
            dxyz(l,3) = 1 - D(iD(h),2);
        end
        if ~isempty(CellArray{D(iD(h),1)+dxyz(l,1),  D(iD(h),2)+dxyz(l,2), D(iD(h),3)+dxyz(l,3)}) % Zelle gefüllt oder nicht
            while ~isempty(CellArray{D(iD(h),1),  D(iD(h),2), D(iD(h),3)}(1,1).Prev) % zurücksetzen
                CellArray{D(iD(h),1), D(iD(h),2), D(iD(h),3)}(1,1) = CellArray{D(iD(h),1), D(iD(h),2), D(iD(h),3)}(1,1).Prev;
            end
            while ~isempty(CellArray{D(iD(h),1)+dxyz(l,1),  D(iD(h),2)+dxyz(l,2), D(iD(h),3)+dxyz(l,3)}(1,1).Next)
                xyz =[xyz; CellArray{D(iD(h),1) + dxyz(l,1) , D(iD(h),2) + dxyz(l,2), D(iD(h),3) + dxyz(l,3)}(1,1).Next.coordinates];
                CellArray{D(iD(h),1) + dxyz(l,1) , D(iD(h),2) + dxyz(l,2), D(iD(h),3) + dxyz(l,3)}(1,1) = CellArray{D(iD(h),1) + dxyz(l,1) , D(iD(h),2) + dxyz(l,2), D(iD(h),3) + dxyz(l,3)}(1,1).Next;
            end
        end
    end
    n = size(xyz, 1);
    F_neu = zeros(n,3);
    F = zeros(n,3);
    r = zeros(n,4);
    for n = 1:n
        r(n,:) = xyz(n,:)-coordinates(h,:); % Abstand zu jedem Teilchen k in Zelle aus CellArray
        if norm(r(n,:)) > 2.5*sigma_OO
            r(n,:) = zeros(1,4);
        end
        if coordinates(h,1) == q_O && r(n,1) == 0
            F_neu(n,:) = 24*epsilon_OO/sum(r(n,2:4).^2, 2) * sigma_OO/((sum(r(n,2:4).^2, 2))^3) * (1-2*sigma_OO^6/(sum(r(n,2:4).^2,2)^3)) * r(n,2:4);
        elseif coordinates(h,1) == q_O && r(n,1) ~= 0 || coordinates(h,1) == q_H && r(n,1) ~= 0
            F_neu(n,:) = 24*epsilon_OH/sum(r(n,2:4).^2, 2) * sigma_OH/((sum(r(n,2:4).^2, 2))^3) * (1-2*sigma_OH^6/(sum(r(n,2:4).^2,2)^3)) * r(n,2:4);
        elseif coordinates(h,1) == q_H && r(n,1) == 0
            F_neu(n,:) = 24*epsilon_HH/sum(r(n,2:4).^2, 2) * sigma_HH/((sum(r(n,2:4).^2, 2))^3) * (1-2*sigma_HH^6/(sum(r(n,2:4).^2,2)^3)) * r(n,2:4);
        end
        F_neu(isnan(F_neu)) = 0;
    end
    F_0(h,:) = sum(F_neu, 1);
    for i = 1:np_Teilchen
        while ~isempty(CellArray{D(iD(i),1),  D(iD(i),2), D(iD(i),3)}(1,1).Prev)
            CellArray{D(iD(i),1), D(iD(i),2), D(iD(i),3)}(1,1) = CellArray{D(iD(i),1), D(iD(i),2), D(iD(i),3)}(1,1).Prev;
        end
    end
end
end

function F_Coulomb = Coulomb_Kraft(coordinates, q_O, q_H, np_Teilchen, sigma_OO, Atome_dim)
F_Coulomb = zeros(np_Teilchen,3);
r = zeros(np_Teilchen, 4);

d = [1.5*sigma_OO*Atome_dim,1.5*sigma_OO*Atome_dim,1.5*sigma_OO*Atome_dim];
epsilon_0 = 1;

for i=1:np_Teilchen
    F_neu = zeros(np_Teilchen,3, np_Teilchen);
    r(:,:,i) = bsxfun(@minus, coordinates(i,:), coordinates);
    r(:,2:4,i) = r(:,2:4,i) - d.*round(r(:,2:4,i)./d);
    if coordinates(i,1) == q_O && r(i,1) == 0
        F_neu(:,:,i) = 1/(4*pi*epsilon_0).*q_O^2./(sum(r(:,2:4,i).^2, 2).^1.5) .* r(:,2:4,i);
    elseif coordinates(i,1) == q_O && r(i,1) ~= 0 || coordinates(i,1) == q_H && r(i,1) ~= 0
        F_neu(:,:,i) = 1/(4*pi*epsilon_0).*q_O*q_H./(sum(r(:,2:4,i).^2, 2).^1.5) .* r(:,2:4,i);
    elseif coordinates(i,1) == q_H && r(i,1) == 0
        F_neu(:,:,i) = 1/(4*pi*epsilon_0).*q_H^2./(sum(r(:,2:4,i).^2, 2).^1.5) .* r(:,2:4,i);
    end
    F_neu(isinf(F_neu)) = 0;
    F_neu(isnan(F_neu)) = 0;
end
F_Coulomb(:,:) = sum(F_neu, 3);
end

function E_pot_total = LJ_Pot(xyz, sigma, E, a)
n = length(xyz);
E_pot = zeros(n,1,n);
r = zeros(n,3,n);
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

function [xyz_all, v_all, F_all, E_pot_all, E_kin_all, T_all] = Zeitintegration(CellArray, coordinates, v, t, delta_t, F_0, n_steps, m, sigma, E, np_Teilchen, a, tau,d,C);

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
 
    while sum(bsxfun(@ge, zeros(length(C_0), 3), C_0), 'all') > 0
        C_0 = C_0 + 1;
    end

    for i = 1:np_Teilchen
        if C_0(i,:)/ C(i,:) ~= 1
            PL(i).removeNode;
            PL(i).insertAfter(particle_head(C_0(i,1), C_0(i,2), C_0(i,3))); % neue Verknüpfung
        end
    end

    F_old = F_0;
    F_0 = LJ_Kraft(CellArray, sigma, E, np_Teilchen, D, iD, coordinates);

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

function Generate_xyz(M, Dateiname, q_O, q_H)
% Output als xyz Datei
n = size(M, 1);
% Atomanzahl = repmat({'Ar'},n,1);
fileID = fopen(Dateiname, 'w'); 
Element = strings(n,1);

for i = 1:n
    if M(i,1) == q_O
        Element(i) = 'O';
    elseif M(i,1) == q_H
        Element(i) = 'H';
    end
end

for i=1:size(M,3)
    fprintf(fileID, '%d\n', n);
    fprintf(fileID, 'Moleküldynamik_3D\n');
    for j = 1:n
        fprintf(fileID, '%s %f %f %f\n', Element{j}, M(j,2,i), M(j,3,i), M(j,4,i));
    end
end

fclose(fileID);
end
