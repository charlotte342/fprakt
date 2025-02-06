%Anfangsbedingungen (2D)

x_N = [0, 0; 0, 1; 0, 5.36; 34.75, 0]; % Ort
v_N = [0, 0; -1, 0; -0.425, 0; 0, 0.0296]; % Geschwindigkeit
m_N = [1; 3e-6; 9.55e-4; 1e-14]; % Masse

% Berechnung der Ausgangskräfte
n = length(x_N);

F = zeros(n,2);
tic
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

delta_t = 0.1;
t_end = delta_t*1e4;
t = 0;
x_positions = zeros(n,2,t_end/delta_t);
x_Nneu = zeros(n,2);
v_Nneu = zeros(n,2);
iteration = 0;

while t < t_end
    t_start = tic;
    t = t + delta_t;
    iteration = iteration + 1;
    for i = 1:n
        x_Nneu(i,:) = x_N(i,:) + delta_t*(v_N(i,:) + F(i,:)*delta_t*0.5/m_N(i));
    end
    F_neu = zeros(n,2);
    for i=1:n
        for j=1:n
            if j~=i
                F_neu(i,:) = F_neu(i,:) + Grav_Pot(x_Nneu(i,:), x_Nneu(j,:), m_N(i), m_N(j));
            end
        end
    end
    for i=1:n
        v_Nneu(i,:) = v_N(i,:) + 0.5/m_N(i)*(F(i,:) + F_neu(i,:))*delta_t;
    end
    for i = 1:n
        x_positions(i, :, iteration) = x_Nneu(i, :);
    end
    F = F_neu;
    x_N = x_Nneu;
    v_N = v_Nneu;

end
figure;
hold on;
for s = 1:iteration
    for i = 1:n
        plot(x_positions(i,1,s), x_positions(i,2,s), C{i}, 'LineWidth', 1);
    end
end
title('Planetenbewegung');
legend('s = Sonne', 'b = Erde', 'r = Jupiter', 'g = Haley');
grid on;

toc

