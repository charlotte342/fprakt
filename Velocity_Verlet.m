while t < t_end
    t = t + delta_t;
    iteration = iteration + 1;
    Positions = Positions + delta_t*(v + F.*delta_t*0.5/m);
    F_neu = LJ_Kraft(Positions, sigma, E);
    v = v + bsxfun(@rdivide, (F+F_neu), 2*m)*delta_t;
    xyz_positions(:, :, iteration) = Positions(:, :);
    F = F_neu;
end