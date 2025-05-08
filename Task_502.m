%%%%%%%%%%%%%%%%%%%%%%%%
% Wasser - Ewald Summe %
%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Teilchen in Box
Atome_dim = 2;

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

%% Energieberechnung - kurzreichweitig

n_cells_x = size(CellArray, 1);
n_cells_y = size(CellArray, 2);
n_cells_z = size(CellArray, 3);

x = -1:1;
y = -1:1;
z = -1:1;

[X,Y,Z] = meshgrid(x,y,z);
E_LJ = zeros(np_Teilchen, 1);
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
    E_neu = zeros(n,1); % Hilfe für E_LJ
    E = zeros(n,1); % Hilfe für E_kr

    for i = 1:n
        r(i,:) = xyz(i,:)-coordinates_0(h,:); % Abstand zu jedem Teilchen k in Zelle aus CellArray
        r_betrag(i) = norm(r(i,2:4));
    end
    for i = 1:n
        if r_betrag(i) > 2.5*sigma_OO
            E_neu(i) = 0;
            continue
        end
        if r_betrag(i) == 0
            E_neu(i) = 0;
            E(i) = 0;
        elseif coordinates_0(h,1) == q_O && r(i,1) == 0
            E_neu(i) = a*epsilon_OO*((sigma_OO.^12)./((sum(r(i, 2:4).^2,2).^6))-(sigma_OO.^6)./(sum(r(i,2:4).^2,2).^3));
            E(i) = .5*skal*q_O^2*erfc(G*r_betrag(i))/r_betrag(i);
        elseif coordinates_0(h,1) == q_O && r(i,1) ~= 0 || coordinates_0(h,1) == q_H && r(i,1) ~= 0
            E_neu(i) = a*epsilon_OH*((-sigma_OH.^6)./(sum(r(i,2:4).^2,2).^3)+(sigma_OH.^12)./((sum(r(i, 2:4).^2,2).^6)));
            E(i)= .5*skal*q_O*q_H*erfc(G*r_betrag(i))/r_betrag(i);
        elseif coordinates_0(h,1) == q_H && r(i,1) == 0
            E_neu(i) = a*epsilon_HH*((-sigma_HH.^6)./(sum(r(i,2:4).^2,2).^3)+(sigma_HH.^12)./((sum(r(i, 2:4).^2,2).^6)));
            E(i) = .5*skal*q_H^2*erfc(G*r_betrag(i))/r_betrag(i);
        end
    end
    E_LJ(h,:) = sum(E_neu, 1);
    E_kr(h,:) = sum(E,1);
    for i = 1:np_Teilchen
        while ~isempty(CellArray{D(iD(i),1),  D(iD(i),2), D(iD(i),3)}(1,1).Prev)
            CellArray{D(iD(i),1), D(iD(i),2), D(iD(i),3)}(1,1) = CellArray{D(iD(i),1), D(iD(i),2), D(iD(i),3)}(1,1).Prev;
        end
    end
end
V_LJ = .5*sum(E_LJ,1);
V_kr = .5*sum(E,1);

%% langreichweitiger Term

E = zeros(np_Teilchen, 3);
L = 4*1.5*sigma_OO; % Boxlänge
n = -1*round(.5*np_Teilchen); % n Element Z
m = 0; % Laufindex
while n < round(.5*np_Teilchen) % nr = nk
    n = n+1;
    k = (2*pi)/L * n;
    m = m+1;
    if k == 0
        E(m,:) = zeros(1,3);
        continue
    end
    rho_x = zeros(i,1);
    rho_y = zeros(i,1);
    rho_z = zeros(i,1);
    for i = 1:np_Teilchen
        rho_x(i) = coordinates_0(i,1) * exp(1j*k*coordinates_0(i,2));
        rho_y(i) = coordinates_0(i,1) * exp(1j*k*coordinates_0(i,3));
        rho_z(i) = coordinates_0(i,1) * exp(1j*k*coordinates_0(i,4));
    end
    E(m,1) = 1/(epsilon_0*k^2) * norm(sum(rho_x,1))^2*exp(-k^2/(4*G^2));
    E(m,2) = 1/(epsilon_0*k^2) * norm(sum(rho_y,1))^2*exp(-k^2/(4*G^2));
    E(m,3) = 1/(epsilon_0*k^2) * norm(sum(rho_z,1))^2*exp(-k^2/(4*G^2));
end
V_lr = 1/(2*L^3) * sum((E(:,:)), "all");

%% Selbstwechselwirkung
E_self = zeros(np_Teilchen,1);
for i =1:np_Teilchen
    E_self(i) = coordinates_0(i,1)^2;
end
V_self = skal*G/sqrt(pi) * sum(E_self);

% Gesamtenergie

V_Coulomb = V_lr - V_self + V_kr;
V = V_Coulomb + V_LJ;

%% Funktionen

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


function [coordinates_0, v, d, np_Teilchen] = Initialisierung_PBC(d_x, d_y, d_z, np_Teilchen, sigma_OO, theta0_HOH, r0_OH, q_O, q_H)
% Initialisierung
% zentrierte Box
d = [d_x, d_y, d_z];
% np zufällig verteilte Punkt-Teilchen in den Grenzen der Box

r = zeros(np_Teilchen,3,np_Teilchen);
r_betrag = zeros(np_Teilchen,1,np_Teilchen);
% Initialisierte Zufallszahlen mit 1239465719 für Vergleichbarkeit mit rng
rng(1239465719);
xyz_0 = rand (np_Teilchen, 3); % rand: Zufallszahlen zwischen 0 und 1
np_oxygen = length(xyz_0);
coordinates_O = bsxfun(@minus, bsxfun(@times, xyz_0, d), 0.5*d); % Zufallszahlen im Intervall von d

n = length(xyz_0);

% r_min implementieren
for i = 1:n
    r(:,:,i) = bsxfun(@minus, coordinates_O, coordinates_O(i,:));
    r(:,:,i) = r(:,:,i) - d.* round(r(:,:,i)./d); % PBC: wenn Abstand zu Teilchen in nächster Box kürzer
end
for i = 1:n
    for j = i:n
        r_betrag(j,1,i) = norm(r(j,:,i));
        r_betrag(i,1,j) = r_betrag(j,1,i);
        r_betrag(j,1,j) = 2*sigma_OO;
        while r_betrag(j,1,i) < 2*sigma_OO
            xyz_0(j,:) = rand(1,3);
            coordinates_O(j,:) = bsxfun(@minus, bsxfun(@times, xyz_0(j,:), d), 0.5*d);
            r(:,:,i) = bsxfun(@minus, coordinates_O(i,:), coordinates_O);
            r(:,:,i) = r(:,:,i) - d.* round(r(:,:,i)./d);
            r_betrag(j,1,i) = norm(r(j,:,i));
            r_betrag(i,1,j) = r_betrag(j,1,i);
        end
    end
end     
for i = 1:np_oxygen
    coordinates_H(i,:) = [coordinates_O(i,1) - cos(theta0_HOH - 90) * r0_OH, coordinates_O(i,2) + sin(theta0_HOH - 90)*r0_OH, coordinates_O(i,3)];
    coordinates_H(np_oxygen+i,:) = [coordinates_O(i,1) - r0_OH, coordinates_O(i,2:3)];
end
coordinates_0 = [repmat(q_O, np_oxygen, 1), coordinates_O; repmat(q_H, 2*np_oxygen, 1), coordinates_H];
v = zeros(n,3);
end

function [V_LJ, V_Coulomb] = Energie(xyz, sigma_OO, sigma_OH, sigma_HH, epsilon_OO, epsilon_OH, epsilon_HH,a, skal, q_O, q_H)
n = length(xyz);
E_LJ = zeros(n,1);
E_Coulomb = zeros(n,1);

for h = 1:n
    r = zeros(n,4);
    r_betrag = zeros(n,1);
    E_neu = zeros(n,1);
    E = zeros(n,1);
    for i = 1:n
        r(i,:) = xyz(i,:)-xyz(h,:);
        r_betrag(i) = norm(r(i,2:4));
    end
    for i = 1:n
        if r_betrag(i) == 0
            E_neu(i) = 0;
            E(i) = 0;
        elseif xyz(h,1) == q_O && r(i,1) == 0 % O-O
            E_neu(i) = a*epsilon_OO*(((sigma_OO.^12)./(sum(r(i,2:4).^2,2).^6))-((sigma_OO.^6)./((sum(r(i, 2:4).^2,2).^3))));
            E(i) = skal*q_O^2/r_betrag(i);
        elseif xyz(h,1) == q_O && r(i,1) ~= 0 || xyz(h,1) == q_H && r(i,1) ~= 0 % O-H
            E_neu(i) = a*epsilon_OH*(((sigma_OH.^12)./(sum(r(i,2:4).^2,2).^6))-((sigma_OH.^6)./((sum(r(i, 2:4).^2,2).^3))));
            E(i) = skal*q_O*q_H/r_betrag(i);
        elseif xyz(h,1) == q_H && r(i,1) == 0 % H-H
            E_neu(i) = a*epsilon_HH*(((sigma_HH.^12)./(sum(r(i,2:4).^2,2).^6))-((sigma_HH.^6)./((sum(r(i, 2:4).^2,2).^3))));
            E(i) = skal*q_H^2/r_betrag(i);
        end
        if r_betrag(i) > 2.5*sigma_OO % Cut-off nur für LJ
            E_neu(i) = 0;
        end
    end
    E_LJ(h,:) = sum(E_neu, 1);
    E_Coulomb(h,:) = sum(E,1);
end
V_LJ = .5 * sum(E_LJ, 1);
V_Coulomb = .5*sum(E_Coulomb,1);

end
