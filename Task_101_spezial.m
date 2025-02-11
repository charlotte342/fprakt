% Anfangsbedingungen (2D)

x_N = [0, 0; 0, 1; 0, 5.36; 34.75, 0]; % Ort
v_N = [0, 0; -1, 0; -0.425, 0; 0, 0.0296]; % Geschwindigkeit
m_N = [1; 3e-6; 9.55e-4; 1e-14]; % Masse

% Berechnung der Ausgangskräfte
n = length(x_N);

F = zeros(n,2);

for i=1:n
    for j=1:n
        if j~=i
            F(i,:) = F(i,:) + Grav_Pot(x_N(i,:), x_N(j,:), m_N(i), m_N(j));  
        end
    end
end

figure;
C = {'k*:','b*:','r*:','g*:'};
for i=1:n
    plot(x_N(i,1), x_N(i,2),C{i}, 'LineWidth', 1);
    hold on
end
hold on

h = 0.015; % Zeitschritt;
t_end = h*10000;
t = 0;
x_Nneu = zeros(n,2);
v_Nneu = zeros(n,2);
x_positions = zeros(n,2,t_end/h);
iteration = 0;

k_1 = zeros(n,2);
k_2 = zeros(n,2);
k_3 = zeros(n,2);
k_4 = zeros(n,2);
M = zeros(n,2);

tic
while t < t_end
    t = t + h;
    iteration = iteration + 1;  
    %     Berechnung der vier Zwischenschritte
    k1_v = h*F./m_N;
    k1_x = h*v_N;
    F1 = zeros(n,2);
    for i=1:n
        for j=1:n
            if j~=i
                F1(i,:) = F1(i,:) + Grav_Pot((x_N(i,:) + 0.5*k1_x(i,:)), (x_N(j,:) + 0.5*k1_x(j,:)), m_N(i), m_N(j));
            end
        end
    end

    k2_v = h*F1./m_N;
    k2_x = h*(v_N + 0.5*k1_v);

    F2 = zeros(n,2);
    for i=1:n
        for j=1:n
            if j~=i
                F2(i,:) = F2(i,:) + Grav_Pot((x_N(i,:) + 0.5*k2_x(i,:)), (x_N(j,:) + 0.5*k2_x(j,:)), m_N(i), m_N(j));
            end
        end
    end

    k3_v = h* F2./m_N;
    k3_x = h*(v_N + 0.5*k2_v);
    
    F3 = zeros(n,2);
    for i=1:n
        for j=1:n
            if j~=i
                F3(i,:) = F3(i,:) + Grav_Pot((x_N(i,:) + k3_x(i,:)), (x_N(j,:) + k3_x(j,:)), m_N(i), m_N(j));
            end
        end
    end

    k4_v = h*F3./m_N;
    k4_x = h*(v_N + k3_v);
    
    % Aktualisierte Positionen und Geschwindigkeiten
    x_Nneu = x_N + (1/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x);
    v_Nneu = v_N + (1/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v);

    for i = 1:n
        x_positions(i, :, iteration) = x_Nneu(i, :);
    end
        F_neu = zeros(n,2);
    for i=1:n
        for j=1:n
            if j~=i
                F_neu(i,:) = F_neu(i,:) + Grav_Pot(x_Nneu(i,:), x_Nneu(j,:), m_N(i), m_N(j));
            end
        end
    end

    v_Nneu(:,:) = v_N(:,:) + bsxfun(@rdivide, (F+F_neu), 2*m_N)*h;

    F = F_neu;
    x_N = x_Nneu;
    v_N = v_Nneu;
end

for i = 1:iteration
    for j = 1:n
        plot(x_positions(j,1,i), x_positions(j,2,i), C{j}, 'LineWidth', 1);
    end
end
title('Planetenbewegung');
legend('s = Sonne', 'b = Erde', 'r = Jupiter', 'g = Haley');
grid on;

toc

