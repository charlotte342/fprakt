% Anfangsbedingungen (2D)

x_N = [0, 0; 0, 1; 0, 5.36; 0, 34.75]; % Ort
v_N = [0, 0; -1, 0; -0.425, 0; 0, 0.0296]; % Geschwindigkeit
m_N = [1; 3e-6; 9.55e-4; 1e-14]; % Masse

% für mich für die manuelle Überprüfung meiner Werte
x_Sonne = x_N(1,1:2);
x_Erde = x_N(2,1:2);
x_Jupiter = x_N(3,1:2);
x_Haley = x_N(4,1:2);

% Berechnung der Ausgangskräfte
n = length(x_N);
r = zeros(n);
F = zeros(n);
for i=1:n
    for j=1:n
        if i==j
            F_x(i,j) = 0;
            F_y(i,j) = 0;
        else
            r_x(i,j) = x_N(i,1) - x_N(j,1);
            r_y(i,j) = x_N(i,2) - x_N(j,2);
            F_x(i,j) = Grav_Pot((r_x(i,j)), m_N(i), m_N(j));
            F_y(i,j) = Grav_Pot((r_y(i,j)), m_N(i), m_N(j));

            if r_x == 0
                F_x = 0;
            end
            if r_y == 0
                F_y = 0;
            end
        end
    end
end
F_0 = [sum(F_x,2), sum(F_y,2)]
% disp(r(1,2))
% disp(F(1,2))
% F_man = Grav_Pot(r_y(2,3), m_N(2), m_N(3))

% Velocity Verlet Propagator: neue x, v und F
delta_t = 5; % Zeitschritt
t_end = 5;
t = 0; % Startzeit
F_neu = zeros(n,2);

figure;
C = {'k*:','b*:','r*:','g*:'};
for i=1:n
    plot(x_N(i,1), x_N(i,2),C{i}, 'LineWidth', 2);
    hold on
end
hold on
while t <= t_end
    t = t + delta_t;
    for i = 1:n
        x_N(i,:) = x_N(i,:) + delta_t*v_N(i,:) + F_0(i,:)*(delta_t)^2*0.5/m_N(i)
        plot(x_N(i,1), x_N(i,2), C{i}, 'LineWidth', 2);
        title(['Vergleich der numerischen Performance'],['von naivemult und optimizedmm']);
        legend('y1=naivemult', 'y2=optimizedmm');
        xlabel('Größe der Matrix'); ylabel('Zeit / s'); grid on;
    end
    for i=1:n
        for j=1:n
            if i==j
                F_x(i,j) = 0;
                F_y(i,j) = 0;
            else
                r_x(i,j) = x_N(i,1) - x_N(j,1);
                r_y(i,j) = x_N(i,2) - x_N(j,2);
                F_x(i,j) = Grav_Pot((r_x(i,j)), m_N(i), m_N(j));
                F_y(i,j) = Grav_Pot((r_y(i,j)), m_N(i), m_N(j));
            end
            if r_x == 0
                F_x = 0;
            end
            if r_y == 0
                F_y = 0;
            end
        end
    end
    F_neu = [sum(F_x,2), sum(F_y,2)]
    for i=1
        v_N(i,:) = v_N(i,:) + 0.5/m_N(i)*(F_0(i,:) + F_neu(i,:)*delta_t)
    end
    F_0 = F_neu;
    disp(v_N);
    disp(x_N);
    M(i) = getframe;
end
movie(M)
%Test
x_Haley_neu = x_Haley + delta_t*v_N(4,:) + F_0(4,:)*(delta_t)^2*0.5/m_N(4,:);
x_Haley_neuneu = x_Haley_neu + delta_t*v_N(4,:) + F_0(4,:)*(delta_t)^2*0.5/m_N(4,:);