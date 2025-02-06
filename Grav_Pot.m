function F = Grav_Pot(x1, x2, m_1, m_2);
% Diese Funktion berechnet das Gravitationspotential zweier Körper in
% Abhängigkeit von den Vektoren x1 und x2

G = 1;
% G = 6.6743015e-11/1.98892e30*1/149597870700^3; % Gravitationskonstante
% r = x_1-x_2
r = norm(x2-x1);
r_dir = x2 - x1;
F = G*m_1*m_2/r^3*r_dir; 
% F = [F_x, F_y];


