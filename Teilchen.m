classdef Teilchen
    properties
        coordinates_0
%         velocities_0
        a
        sigma
        E
        m

    end
    methods
        % initialize, propagate...
        function obj=Teilchen(coordinates_0, a, sigma, E, m) % heißt wie Klasse
            % Input
            obj.a=a;
            obj.sigma=sigma;
            obj.E=E;
            obj.m=m;
            obj.coordinates_0=coordinates_0;
        end
    end
end
