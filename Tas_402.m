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
% n_cells = 9; % in 2D
% N_cell_x = 3;
% N_cell_y = 3;
% Cell = zeros(n_cells, 1);
% 
N_x = 0 : 2.5*sigma : 2.5*sigma*3;
N_y = 0 : 2.5*sigma : 2.5*sigma*3;

[X,Y] = meshgrid(N_x,N_y);
Cell_coordinates = [X(:), Y(:)];

% When using the LC algorithm, using a 1D index representation as cell identifiers instead of
% coordinates can be beneficial. The purpose is to efficiently store and access elements in a
% multi-dimensional grid using a one-dimensional data structure like an array.
% Having a 3D Cartesian grid with dimensions dimx × dimy × dimz , where dimx, dimy,
% dimz represent the number of cells along the x, y, and z-axis, respectively, the bijective
% mapping f from the 3D grid coordinates x, y, z of a cell to a unique 1D index is given by:

% 1D index representation as cell identifiers: f = x ∗ (dimx)^2 + y ∗ dimy + z 
% Cell = Cell_coordinates(:,1) * 3^3 + Cell_coordinates(:,2) * 3;


% Datenstruktur für jedes Teilchen N --> aus Unterklasse, die dann verlinkt
% werden in Zellen
C = zeros(np_Teilchen, 2);
for i = 1:np_Teilchen
    N(i) = Molekueldynamik(coordinates_0(i,:), 4, sigma(i), E(i), m(i), 0, 1e-3, 1, 1, 50);
    C(i,:) = round(N(i).coordinates_0./(2.5*sigma(i))) + 1;
end


for i = 1:numParticles
    row = C(i, 1);
    col = C(i, 2);
    
    % Add the particle to the corresponding cell
    dlnode_ex{row, col} = dlnode_ex{row, col}(N(i));
end

% Desired number of rows per pair
n = 3;

% Initialize the cell array for storing the output structure
unique_pairs = unique(C, 'rows', 'stable'); % Get unique pairs
output_cell_array = cell(size(unique_pairs, 1), 1); % Initialize the output cell array

% Loop over each unique pair
for i = 1:size(unique_pairs, 1)
    % Find the rows corresponding to the current unique pair
    rows = find(ismember(C, unique_pairs(i, :), 'rows'));
    
    % Now create a new cell array for this unique pair
    temp_cell = cell(n, 1); % Create a cell array with n rows
    
    % Fill the temp cell array with row indices (or NaN if less than n rows)
    num_rows = numel(rows);
    for j = 1:n
        if j <= num_rows
            temp_cell{j} = rows(j); % Store the row index
        else
            temp_cell{j} = NaN; % Pad with NaN if less than n rows
        end
    end
    
    % Store the temp cell array in the output cell array
    output_cell_array{i} = temp_cell;
end

% Display the result
disp('Output cell array structure:');
for i = 1:numel(output_cell_array)
    fprintf('Unique pair: [%d, %d], Cell content:\n', unique_pairs(i, 1), unique_pairs(i, 2));
    disp(output_cell_array{i});
end


% Initialize a cell array to store the output matrices
unique_pairs = unique(C, 'rows', 'stable'); % Get unique pairs
output = cell(size(unique_pairs, 1), 1); % Initialize cell array for each unique pair

% Loop over each unique pair and store the row indices
for i = 1:size(unique_pairs, 1)
    % Find the rows corresponding to the current unique pair
    rows = find(ismember(C, unique_pairs(i, :), 'rows'));
    
    % Store the row indices in the output cell array
    output{i} = rows;
end

% Display the results
disp('Output structure:');
for i = 1:numel(output)
    fprintf('Unique pair: [%d, %d], Rows: ', unique_pairs(i, 1), unique_pairs(i, 2));
    disp(output{i}');
end
% 
% for i = 1: 3
%     for j = 1 : 3
%         if 
%             C(i,j) = {N(i)}
%         end
%     end
% end

%%
% C(i) = 
% C = zeros(3,3); % 9 Zellen insgesamt
% for i = 1:3
%     
%     for j = 1:3
%         if N(i).coordinates_0(i,:) <= Cell_coordinates(i,:)
%             C(i,j) = Molekueldynamik()
%         end
%     end
% end

% %% w/ example dlnode code
% sigma = 1;
% np_Teilchen = 10*10; 
% coordinates_0 = Anfangspositionen(sigma, sqrt(np_Teilchen)-1);
% 
% for i = 1:np_Teilchen
%     n(i) = dlnode_ex(coordinates_0(i));
% %  obj=dlnode(coordinates_0, a, sigma, E, m, t, delta_t, n_steps, tau, T_0)
% end
% for i = 1:n_cells
%     for j = 1:np_Teilchen
%         if coordinates
%         C(i) = {n(j)}
%     end
% end



%% Matlab code for the generation of linked lists via dlnode-class
n1 = dlnode(1);
n2 = dlnode(2);
n3 = dlnode(3);

% Build these nodes into a doubly linked list using the class methods:

n2.insertAfter(n1) % Insert n2 after n1
n3.insertAfter(n2) % Insert n3 after n2

% Navigate through the list:
n1.Next % Points to n2

n2.Next.Pref % Points back to n2

n1.Next.Next % Points to n3

n3.Prev.Prev % Points to n1


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
