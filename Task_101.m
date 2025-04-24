tic
% Anfangsbedingungen (2D)

x_N = [0, 0; 0, 1; 0, 5.36; 34.75, 0]; % Ort
v_N = [0, 0; -1, 0; -0.425, 0; 0, 0.0296]; % Geschwindigkeit
m_N = [1; 3e-6; 9.55e-4; 1e-14]; % Masse

n = length(x_N);



% Berechnung der Ausgangskräfte
F = Grav_Pot_vektorisiert(x_N,m_N);

figure;
C = {'k*:','b*:','r*:','g*:'};
for i=1:n
    plot(x_N(i,1), x_N(i,2),C{i}, 'LineWidth', 1);
    hold on
end
hold on

delta_t = 0.1;
t_end = delta_t*10;
j = t_end/delta_t;
t = 0;
x_positions = zeros(n,2,j);
iteration = 0;

while t < t_end
    t = t + delta_t;
    iteration = iteration + 1;
    x_N = x_N + delta_t*(v_N + F.*delta_t*0.5./m_N);
    F_neu = Grav_Pot_vektorisiert(x_N, m_N);
    v_N  = v_N + bsxfun(@rdivide, (F+F_neu), 2*m_N)*delta_t;
    x_positions(:, :, iteration) = x_N;
    F = F_neu;
end

for s = 1:iteration
    for i = 1:n
        plot(x_positions(i,1,s), x_positions(i,2,s), C{i}, 'LineWidth', 1);
        hold on
    end
end
title('Planetenbewegung');
legend('s = Sonne', 'b = Erde', 'r = Jupiter', 'g = Haley');
grid on;
hold off

toc

%% Section X: Funktionen
function F_total = Grav_Pot_vektorisiert(x, m)
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
    F_neu(:,:,i) = bsxfun(@times, [1*m(i)*m./(r_betrag(:,:,i).^3)], r(:,:,i));
    F = F_neu;
    F(isnan(F_neu)) = 0;
end
F_total(:,:) = sum(F, 3);
end

