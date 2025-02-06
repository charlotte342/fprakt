function V_LJ = LJ_Pot(a, E, sigma, x_1, x_2)

r_ij = norm(x_1 - x_2)
V_LJ(r_ij) = a*E*((sigma/r_ij)^n-(sigma/r_ij)function ^m);