%% Initialisierung
% Parameter
sigma = 1;
E = 1;
m = 1; % Masse für H-Atom, 1.00784 u

% Ausgangspositionen und Speicherzuordnung
x = 0 : 1.5*sigma : 10.5;
y = 0 : 1.5*sigma : 10.5;
[X,Y] = meshgrid(x,y);
Positions = [X(:), Y(:)];
n = length(Positions);

F = LJ_Kraft(Positions, sigma, E);

t = 0;
delta_t = 0.1;
t_end = delta_t*3;
iteration = 1;
iteration_total = t_end/delta_t;

% Speicherzuordnung
v = zeros(n,2);
xy_positions = zeros(n,2,iteration_total);

figure;
for i=1:n
    plot(Positions(i,1), Positions(i,2), 'ro', 'LineWidth', 0.7);
    hold on
end
hold on

% % Überprüfung des LJ_Pot
% F_neu = zeros(n,2,n);
% r = zeros(n,2,n);
% F_total = zeros(n,2);
% r_betrag = zeros(n,1,n);
% for i=1:n
%     r(:,:,i) = bsxfun(@minus, Positions(i,:), Positions);
%     for j=1:n
%         r_betrag(j,1,i) = norm(r(j,:,i));
%     end
%     F_neu(:,:,i) = bsxfun(@times, (24*E./(r_betrag(:,:,i).^2).*(sigma./r_betrag(:,:,i)).^6.*(1-2*(sigma./r_betrag(:,:,i)).^6)), r(:,:,i));
%     F = F_neu;
%     F(isnan(F_neu)) = 0;
% end
% F_total(:,:) = sum(F, 3);

%% Kraftberechnung und Moleküldynamik

while t < t_end
    t = t + delta_t;
    iteration = iteration + 1;
    Positions = Positions + delta_t*(v + F.*delta_t*0.5/m);
    for i=1:n
        plot(Positions(i,1), Positions(i,2),'bo','LineWidth', 0.7);
        M(i) = getframe;
        hold on
    end
    F_neu = LJ_Kraft(Positions, sigma, E);
    v = v + bsxfun(@rdivide, (F+F_neu), 2*m)*delta_t;
    xy_positions(:, :, iteration) = Positions(:, :);
    F = F_neu;
end
movie(M)

% Output als xyz-file
writematrix('202.xyz', [xy_positions(:,1), xy_positions(:,2), 0], 'Delimeter', '\t')
%% Funktionen

function F_total = LJ_Kraft(x, sigma, E)
n = length(x);
F_neu = zeros(n,2,n);
r = zeros(n,2,n);
F_total = zeros(n,2);
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

% function F_total = LJ_Kraft(x, sigma, E)
% n = length(x);
% F_neu = zeros(n,2,n);
% r = zeros(n,2,n);
% F_total = zeros(n,2);
% r_betrag = zeros(n,2);
% for i=1:n
%     r(:,:,i) = bsxfun(@minus, x(i,:), x);
%     for j=1:n
%         r_betrag(j,:) = norm(r(j,:,i));
%     end
%     F_neu(:,:,i) = bsxfun(@times, (24*E./(r_betrag.^2).*(sigma./r_betrag).^6.*(1-2*(sigma./r_betrag).^6)), r(:,:,i));
%     F = F_neu;
%     F(isnan(F_neu)) = 0;
% end
% F_total(:,:) = sum(F, 3);
% end

