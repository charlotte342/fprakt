%% Initialisierung
% Parameter
sigma = 1;
E = 1;
m = 1; % Masse für H-Atom, 1.00784 u

% Ausgangspositionen und Speicherzuordnung
Positions = [0,0,0; 1.5*sigma, 0, 0];
n = size(Positions, 1);

F = LJ_Kraft(Positions, sigma, E);

t = 0;
delta_t = 0.001;
t_end = delta_t*1000;
iteration = 0;
iteration_total = t_end/delta_t;

% Speicherzuordnung
v = zeros(n,3);
xyz_positions = zeros(n,3,int32(iteration_total));

% plot in matlab
% figure;
% for i=1:n
%     plot3(Positions(i,1), Positions(i,2), Positions(i,3), 'r*', 'LineWidth', 1);
%     hold on
% end
% hold on

%% Kraftberechnung und Moleküldynamik

tic
while t < t_end
    t = t + delta_t;
    iteration = iteration + 1;
    Positions = Positions + delta_t*(v + F.*delta_t*0.5/m);
    % falls plot in matlab
%     for i=1:n
%         plot3(Positions(i,1), Positions(i,2), Positions(i,3), 'bo','LineWidth', 0.7);
%         M(i) = getframe;
%         hold on
%     end
    F_neu = LJ_Kraft(Positions, sigma, E);
    v = v + bsxfun(@rdivide, (F+F_neu), 2*m)*delta_t;
    xyz_positions(:, :, iteration) = Positions(:, :);
    F = F_neu;
end
% movie(M)
toc

% Output als xyz-file 
atome = repmat({'H'},n,1);
fileID = fopen('202_2Atome.xyz', 'w'); 

for i=1:iteration
    fprintf(fileID, '%d\n', n);
    fprintf(fileID, 'Moleküldynamik\n');
    for j = 1:n
        fprintf(fileID, '%s %6.4f %6.4f %6.4f\n', atome{j}, xyz_positions(j,1,i), xyz_positions(j,2,i), xyz_positions(j,3,i));
    end
end

fclose(fileID);
%% Funktionen

function F_total = LJ_Kraft(x, sigma, E)
n = size(x,1);
F_neu = zeros(n,3,n);
r = zeros(n,3,n);
F_total = zeros(n,3);
r_betrag = zeros(n,1,n);
for i=1:n
    r(:,:,i) = bsxfun(@minus, x(i,:), x);
    for j=1:n
        r_betrag(j,1,i) = norm(r(j,:,i));
    end
    F_neu(:,:,i) = bsxfun(@times, (24*E./(r_betrag(:,:,i).^2).*(sigma./r_betrag(:,:,i)).^6.*(1-2*(sigma./r_betrag(:,:,i)).^6)), r(:,:,i));
    F = F_neu;
    F(isnan(F_neu)) = 0;
end
F_total(:,:) = sum(F, 3);
end


