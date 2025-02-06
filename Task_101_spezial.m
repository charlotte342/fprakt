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
t_end = h*1e5;
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
    for i = 1:n
        k_1(i,:) = (x_N(i,:) + h*(v_N(i,:) + F(i,:)*h*0.5/m_N(i))) - x_N(i,:); % k_1 gleich wie x_N, k_1 = f(x,t)
        k_2(i,:) = ((x_N(i,:) + 0.5*h*k_1(i,:)) + (0.5*h)*(v_N(i,:) + F(i,:)*(0.5*h)*0.5/m_N(i))) - x_N(i,:); % k_2 = f(x+0.5*h*k_1, t+0.5*h), Steigung zur Position bei halber Schrittweite
        k_3(i,:) = ((x_N(i,:) + 0.5*h*k_2(i,:)) + (0.5*h)*(v_N(i,:) + F(i,:)*(0.5*h)*0.5/m_N(i))) - x_N(i,:); % k_3 = f(x+0.5*h*k_2, t+0.5*h)
        k_4(i,:) = ((x_N(i,:) + h*k_3(i,:)) + h*(v_N(i,:) + F(i,:)*h*0.5/m_N(i))) - x_N(i,:); % k_4 = f(x+h*k_3, t+h)
        x_Nneu(i,:) = x_N(i,:) + (1/6)*(k_1(i,:)+2*k_2(i,:)+2*k_3(i,:)+k_4(i,:));
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
        v_Nneu(i,:) = v_N(i,:) + 0.5/m_N(i)*(F(i,:) + F_neu(i,:))*h;
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

