function F_LJ = LJ_Grad(E, sigma, x_1, x_2)

r_dir = x_1 - x_2;
r = norm(r_dir);
F_LJ = 24*E/r^2*(sigma/r)^6*(1-2*(sigma/r)^6)*r_dir;