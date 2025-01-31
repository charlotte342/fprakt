%Anfangsbedingungen (2D)

x_N = [0, 0; 0, 1; 0, 5.36; 34.75, 0]; % Ort
v_N = [0, 0; -1, 0; -0.425, 0; 0, 0.0296]; % Geschwindigkeit
m_N = [1; 3e-6; 9.55e-4; 1e-14]; % Masse

% für mich für die manuelle Überprüfung meiner Werte
x_Sonne = x_N(1,1:2);
x_Erde = x_N(2,1:2);
x_Jupiter = x_N(3,1:2);
x_Haley = x_N(4,1:2);

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

delta_t = 0.15;
t_end = delta_t*1e3;
t = 0;
x_Nneu = zeros(n,2);
v_Nneu = zeros(n,2);

while t < t_end
    t = t + delta_t
    for i = 1:n
        x_Nneu(i,:) = x_N(i,:) + delta_t*(v_N(i,:) + F(i,:)*delta_t*0.5/m_N(i));
        plot(x_Nneu(i,1), x_Nneu(i,2), C{i}, 'LineWidth', 1);
        title(['Planetenbewegung']);
        legend('s = Sonne', 'b = Erde', 'r = Jupiter', 'g = Haley');
        grid on;
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
    F = F_neu;
    x_N = x_Nneu;
    v_N = v_Nneu;
%     disp(x_N)
%     disp(v_N)
%     disp(F_0)
    M(i) = getframe;
end

hold off
% disp(v_N);
% disp(x_N);
movie(M);
% % 
