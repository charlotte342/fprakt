classdef MD_Programm
    % coordinates, velocities...

    properties
        a
        sigma
        E
        d
        np_Teilchen
    end
    properties (Dependent)
        coordinates_0
        coordinates
        velocities
        temperature
        energy
        n
        velocities_0
    end
    properties (Constant, Hidden)
        k_B = 3.1651e-06;
    end
    methods
        % initialize, propagate...
        function obj=MD_Programm(a, sigma, E, d, np_Teilchen) % heißt wie Klasse
            % Input
            obj.a=a;
            obj.sigma=sigma;
            obj.E=E;
            obj.d=d;
            obj.np_Teilchen=np_Teilchen;
        end
%         function coordinates_0=set.coordinates_0(obj)
%             coordinates_0 = Initialisierung_PBC(obj.d, obj.np_Teilchen, obj.sigma);
% 
% %             r = zeros(obj.np_Teilchen,3,obj.np_Teilchen);
% %             r_betrag = zeros(obj.np_Teilchen,1,obj.np_Teilchen);
% %             % Initialisierte Zufallszahlen mit 1239465719 für Vergleichbarkeit mit rng
% %             rng(1239465719);
% %             xyz_0 = rand (obj.np_Teilchen, 3); % rand: Zufallszahlen zwischen 0 und 1
% %             coordinates_0 = bsxfun(@minus, bsxfun(@times, xyz_0, obj.d), 0.5*obj.d); % Zufallszahlen im Intervall von d
% % 
% %             % r_min implementieren
% %             for i = 1:obj.n
% %                 r(:,:,i) = bsxfun(@minus, coordinates_0, coordinates_0(i,:));
% %                 r(:,:,i) = r(:,:,i) - .5*obj.d.* round(r(:,:,i)./obj.d); % PBC: wenn Abstand zu Teilchen in nächster Box kürzer
% %             end
% %             for i = 1:obj.n
% %                 for j = i:obj.n
% %                     r_betrag(j,1,i) = norm(r(j,:,i));
% %                     r_betrag(i,1,j) = r_betrag(j,1,i);
% %                     r_betrag(j,1,j) = obj.sigma;
% %                     while r_betrag(j,1,i) < obj.sigma
% %                         xyz_0(j,:) = rand(1,3);
% %                         coordinates_0(j,:) = bsxfun(@minus, bsxfun(@times, xyz_0(j,:), obj.d), 0.5*obj.d);
% %                         r(:,:,i) = bsxfun(@minus, coordinates_0(i,:), coordinates_0);
% %                         r_betrag(j,1,i) = norm(r(j,:,i));
% %                     end
% %                 end
% %             end
%             obj.coordinates_0=coordinates_0;
%         end        
      
        function n=get.n(obj)
%             n=size(obj.coordinates_0,1);
            n = size(obj.d,2);
        end
        function velocities_0=get.velocities_0(obj)
            velocities_0=zeros(size(obj.coordinates_0,1), size(obj_coordinates_0,2));
        end
        function [xyz_all, v, E_kin_all, E_pot_all, E_tot, T_all] = Zeitintegration(xyz, t, delta_t, F_0, n_steps, v, m, sigma, E, d, np_Teilchen, a, tau)
            % Velocity-Verlet-Propagator
            xyz_all = zeros(length(xyz), 3, n_steps);
            E_kin_all = zeros(n_steps,1);
            E_pot_all = zeros(n_steps,1);
            T_all = zeros(n_steps,1);

            iteration = 0;
            v_Betrag = zeros(length(xyz), 1);
            T_0 = 50; % Zieltemperatur
            k_B = 3.1651e-06; % Boltzmann-Konstante in a. u.
            while t < delta_t*n_steps
                t = t + delta_t;
                iteration = iteration + 1;
                xyz = xyz + delta_t*(v + F_0.*delta_t*0.5/m);

                % neu: Periodische Randbedingungen
                for i = 1:size(xyz, 1)
                    for j = 1:size(xyz, 2)
                        if xyz(i,j) >= 0.5*d(j)
                            xyz(i,j) = xyz(i,j) - d(j);
                        elseif xyz(i,j) < -0.5*d(j)
                            xyz(i,j) = xyz(i,j) + d(j);
                        end
                    end
                end
                F_neu = LJ_Kraft(xyz, sigma, E, d);
                E_pot = sum(LJ_Pot(xyz, sigma, E, a, d));

                % Temperaturkontrolle über Skalierung der Geschwindigkeiten mithilfe
                % des Skalierungsfaktors lambda
                v = v + bsxfun(@rdivide, (F_0+F_neu), 2*m)*delta_t;
                for i = 1:size(xyz,1)
                    v_Betrag(i,1) = norm(v(i,:));
                end
                E_kin = sum(.5*m*(v_Betrag.^2));
                T = 2*E_kin./(3*np_Teilchen*k_B);
                lambda = sqrt(1+delta_t/tau*((T_0./T)-1));
                v = v.*lambda;

                xyz_all(:, :, iteration) = xyz(:, :);
                E_kin_all(iteration,1) = E_kin(:,:);
                E_pot_all(iteration,1) = E_pot(:,:);
                T_all(iteration,1) = T(:,:);

                F_0 = F_neu;
            end
            E_tot = E_pot_all + E_kin_all;
        end

    end
end


function xyz = Initialisierung_PBC(d, np_Teilchen, sigma)
% Initialisierung
% zentrierte Box
% np zufällig verteilte Punkt-Teilchen in den Grenzen der Box

r = zeros(np_Teilchen,3,np_Teilchen);
r_betrag = zeros(np_Teilchen,1,np_Teilchen);
% Initialisierte Zufallszahlen mit 1239465719 für Vergleichbarkeit mit rng
rng(1239465719);
xyz_0 = rand (np_Teilchen, 3); % rand: Zufallszahlen zwischen 0 und 1
n = length(xyz_0);
xyz = bsxfun(@minus, bsxfun(@times, xyz_0, d), 0.5*d); % Zufallszahlen im Intervall von d

% r_min implementieren
for i = 1:n
    r(:,:,i) = bsxfun(@minus, xyz, xyz(i,:));
    r(:,:,i) = r(:,:,i) - .5*d.* round(r(:,:,i)./d); % PBC: wenn Abstand zu Teilchen in nächster Box kürzer
end
for i = 1:n
    for j = i:n
        r_betrag(j,1,i) = norm(r(j,:,i));
        r_betrag(i,1,j) = r_betrag(j,1,i);
        r_betrag(j,1,j) = sigma;
        while r_betrag(j,1,i) < sigma
            xyz_0(j,:) = rand(1,3);
            xyz(j,:) = bsxfun(@minus, bsxfun(@times, xyz_0(j,:), d), 0.5*d);
            r(:,:,i) = bsxfun(@minus, xyz(i,:), xyz);
            r_betrag(j,1,i) = norm(r(j,:,i));
        end
    end
end     
v = zeros(size(xyz, 1), 3);
end