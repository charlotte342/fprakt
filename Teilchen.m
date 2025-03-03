classdef Teilchen
    properties
        coordinates_0
%         velocities_0
        a
        sigma
        E
        m
        t
        delta_t
        n_steps
        tau
        T_0 % Zieltemperatur
    end
    methods
        % initialize, propagate...
        function obj=Teilchen(coordinates_0, a, sigma, E, m, t, delta_t, n_steps, tau, T_0) % heißt wie Klasse
            % Input
            obj.a=a;
            obj.sigma=sigma;
            obj.E=E;
            obj.m=m;
            obj.t=t;
            obj.delta_t=delta_t;
            obj.n_steps=n_steps;
            obj.coordinates_0=coordinates_0;
            obj.tau = tau;
            obj.T_0 = T_0;
        end
    end
end
