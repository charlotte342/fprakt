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

% Particle List - head einführen: jede Zelle mit M+1 Teilchen (erstes
% Teilchen kopieren
% Zellen neu anordnen
CellArray = cell(length(D));
for i = 1:length(D)
    CellArray{D(i,1), D(i,2)} = [CellArray_linear{i,1}(1,1),CellArray_linear{i,1}];
end

% Kraftberechnung - z. B. LJ_Kraft oder Coulomb ...
% Linked-Cell: Kraft Fi auf N(i) in Zelle ic aus kc Zellen (Nachbarzellen)
% F = LJ_Kraft(xyz, sigma, E, d);


n_cells_x = 3;
n_cells_y = 3;
% benachbarte Zellen finden:
% alle relevanten Positionen für dieses Teilchen

% [a, sigma, E, m, t, delta_t, n_steps, t_end] = Parameter('Parameter_302.txt');

% tau = delta_t*1e3;
% T_0 = 50;

E = 1;
sigma = 1; 

%%
tic

F_total = zeros(np_Teilchen, 2);
n_particle = zeros(n_cells_x, n_cells_y); % maximal alle Teilchen in einer Zelle
for h = 1:np_Teilchen
    for i= 1:n_cells_x % alle möglichen Zellen
        for j = 1:n_cells_y
            % alle relevanten Positionen für dieses Teilchen
            n_particle(i,j) = size(CellArray{i,j}, 2); % Anzahl Teilchen pro Zelle {i,j}
            n_gesamt = 0;
            for k = -1:1 %  nur benachbarte Zellen, in 2D: 9 Zellen
                for l = -1:1
                    if i+k > 0 && j+l > 0 && i+k <= n_cells_x && j+l <= n_cells_y
                        n_gesamt = n_gesamt + (n_particle(i+k,j+l)-1); % Gesamtzahl relevanter Teilchen für die Kraftberechnung
%                         xyz = zeros(n_gesamt,2);
                        for m = 2:n_particle(i+k,j+l) % ab 2 wegen particle head
                            % Gesamtzahl an relevanten Teilchen
                            coordinate(m-1,:,i+k,j+l) = CellArray{i+k,j+l}(1,m).Data; % relevante Koordinaten für Teilchen in Zelle ij --> Kraftberechnung
%                             F(:,:,h) = LJ_Kraft(reshape(coordinate, size(coordinate,1) * size(coordinate, 3) * size(coordinate, 4), size(coordinate, 2)), 1, 1);
                            %                         neu(:,:) = [coordinate(m,:, k+2, l+2)];
                            xyz = reshape(coordinate, size(coordinate,1) * size(coordinate, 3) * size(coordinate, 4), size(coordinate, 2)); % neu angeordnet
%                             xyz = coordinate;


                        end
                    end
                end
            end

            n = size(xyz, 1);
            F_neu = zeros(n,2);
            r = zeros(n,2);
            r_betrag = zeros(n,1,n);

            for n = 1:n
                for o = 1:n
                    r(o,:) = xyz(n,:) - xyz(o,:); % Abstand n zu jedem Teilchen o
                    F_neu(o,:) = 24*E/sum(r(o,:).^2, 2) * sigma/((sum(r(o,:).^2, 2))^3) * (1-2*sigma^6/(sum(r(o,:).^2,2)^3)) * r(o,:);
                    %                                 r(:,:,n) = bsxfun(@minus, xyz(n,:), xyz);
                    %                                 F_neu(:,:,n) = bsxfun(@times, (24*E./sum(r(:,:,n).^2, 2).*sigma^6./((sum(r(:,:,n).^2, 2)).^3).*(1-2*(sigma^6./(sum(r(:,:,n).^2,2).^3)))), r(:,:,n));
                    F = F_neu;
                    F(isnan(F_neu)) = 0;
                end
                F_h(n,:) = sum(F,1);
            end
            F_total(h,:) = F_h(n,:); % Gesamtkraft auf Teilchen h
        end
    end
end

toc

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


function coordinates_0 = Anfangspositionen(sigma, np_Teilchen)
% Ausgangspositionen und Speicherzuordnung
x = 0 : 1.5*sigma : 1.5*sigma*np_Teilchen;
y = 0 : 1.5*sigma : 1.5*sigma*np_Teilchen;
[X,Y] = meshgrid(x,y);
coordinates_0 = [X(:), Y(:)];
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
