% Matlab Class code dlnode for the implementation of linked lists
classdef dlnode < handle
    % handle: superclass of all handle classes (object, that indirectly
    % references its data)
    % dlnode A class to represent a doubly - linked node .
    % Link multiple dlnode objects together to create linked lists .
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
    properties (Dependent)
        coordinates
        velocities
        temperature
        energy
        E_pot
        E_pot_all
        E_kin_all
        n
        t_end
        F % Force
        %         T_all % Gesamttemperatur
        %         E_kin_all % kinetische Energie, alle Schritte
    end
    properties (Constant, Hidden)
        k_B = 3.1651e-06;
    end
    properties ( SetAccess = private )
        Next = dlnode . empty
        Prev = dlnode . empty
    end
    methods
        % initialize, propagate...
        function obj=dlnode(coordinates_0, a, sigma, E, m, t, delta_t, n_steps, tau, T_0) % heißt wie Klasse
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

        function n=get.n(obj)
            %             n=size(obj.coordinates_0,1);
            n = size(obj.coordinates_0,1);
        end
        function t_end = get.t_end(obj)
            t_end = obj.delta_t*obj.n_steps;
        end
        %         function velocities_0=get.velocities_0(obj)
        %             velocities_0=zeros(size(obj.coordinates_0,1), size(obj.coordinates_0,2));
        %         end

        function F = get.F(obj)
            % Berechnung der Kräfte
            F_neu = zeros(obj.n,3,obj.n);
            r = zeros(obj.n,3,obj.n);
            F = zeros(obj.n,3);
            r_betrag = zeros(obj.n,1,obj.n);
            d = [30,30,30];

            for i=1:obj.n
                r(:,:,i) = bsxfun(@minus, obj.coordinates_0(i,:), obj.coordinates_0);
                r(:,:,i) = r(:,:,i) - d.*round(r(:,:,i)./d);
                for j=i:obj.n
                    r_betrag(j,1,i) = norm(r(j,:,i));
                    r_betrag(i,1,j) = r_betrag(j,1,i);
                end
                F_neu(:,:,i) = bsxfun(@times, (24*obj.E./(r_betrag(:,:,i).^2).*(obj.sigma./r_betrag(:,:,i)).^6.*(1-2*(obj.sigma./r_betrag(:,:,i)).^6)), r(:,:,i));
                F_korr = F_neu;
                F_korr(isnan(F_neu)) = 0;
            end
            F(:,:) = sum(F_korr, 3);
        end

        function E_pot = get.E_pot(obj)
            E_pot_iteration = zeros(obj.n,1,obj.n);
            r = zeros(obj.n,3,obj.n);
            r_betrag = zeros(obj.n,1,obj.n);
            E_pot= zeros(obj.n,1);
            d = [30,30,30];

            for i=1:obj.n
                r(:,:,i) = bsxfun(@minus, obj.coordinates_0(i,:), obj.coordinates_0);
                r(:,:,i) = r(:,:,i) - d.*round(r(:,:,i)./d);
                for j=i:obj.n
                    r_betrag(j,1,i) = norm(r(j,:,i));
                    r_betrag(i,1,j) = r_betrag(j,1,i);
                end
                E_pot_iteration(:,:,i) = obj.a*obj.E*((obj.sigma./r_betrag(:,:,i)).^6-(obj.sigma./r_betrag(:,:,i)).^12);
                E_korr = E_pot_iteration;
                E_korr(isnan(E_pot_iteration)) = 0; % Korrektur für i = j
            end
            E_pot(:,:) = sum(E_korr, 3);
        end


        function result = Zeitintegration(obj)

            % Velocity-Verlet-Propagator
            result.coordinates = zeros(length(obj.coordinates_0), 3, obj.n_steps);
            result.coordinates(:,:,1) = obj.coordinates_0;
            %             obj.E_all = zeros(obj.n_steps, 1);
            result.energy = zeros(obj.n_steps,1);
            result.E_kin_all = zeros(obj.n_steps,1);
            result.E_pot_all = zeros(obj.n_steps,1);
            result.temperature = zeros(obj.n_steps,1);
            velocities_0 = zeros(size(obj.coordinates_0,1), size(obj.coordinates_0,2));
            result.velocities = zeros(length(obj.coordinates_0), 3, obj.n_steps);
            d = [30,30,30];


            iteration = 0;
            v_Betrag = zeros(length(obj.coordinates_0), 1);
            while obj.t < obj.delta_t*obj.n_steps
                obj.t = obj.t + obj.delta_t;
                iteration = iteration + 1;
                F_0 = obj.F;
                xyz = obj.coordinates_0 + obj.delta_t*(velocities_0 + obj.F.*obj.delta_t*0.5/obj.m);

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

                % Temperaturkontrolle über Skalierung der Geschwindigkeiten mithilfe
                % des Skalierungsfaktors lambda
                new_velocities = velocities_0 + bsxfun(@rdivide, (F_0+obj.F), 2*obj.m)*obj.delta_t;
                for i = 1:size(xyz,1)
                    v_Betrag(i,1) = norm(new_velocities(i,:));
                end
                E_kin = sum(.5*obj.m*(v_Betrag.^2));
                T = 2*E_kin./(3*10*obj.k_B);
                lambda = sqrt(1+obj.delta_t/obj.tau*((obj.T_0./T)-1));
                new_velocities = new_velocities.*lambda;
                obj.coordinates_0 = xyz; % Koordinaten aktualisieren für die korrekte Berechnung der potentiellen Energie

                result.E_pot_all(iteration,1) = sum(obj.E_pot(:,:));
                result.E_kin_all(iteration,1) = E_kin;
                result.energy(iteration,1) = E_kin + sum(obj.E_pot(:,:));

                result.temperature(iteration,1) = T(:,:);
                result.velocities(:,:,iteration) = velocities_0;
                velocities_0 = new_velocities;
                result.coordinates(:,:,iteration) = xyz;
            end
        end
    end

    methods
        function node = dlnode ( Data )
            % Construct a dlnode object
            if nargin > 0
                node . Data = Data ;
            end
        end
        function insertAfter ( newNode , nodeBefore )
            % Insert newNode after nodeBefore .
            removeNode ( newNode ) ;
            newNode . Next = nodeBefore . Next ;
            newNode . Prev = nodeBefore ;
            if ~ isempty ( nodeBefore . Next )
                nodeBefore . Next . Prev = newNode ;
            end
            nodeBefore . Next = newNode ;
        end
        function insertBefore ( newNode , nodeAfter )
            % Insert newNode before nodeAfter .
            removeNode ( newNode ) ;
            newNode . Next = nodeAfter ;
            newNode . Prev = nodeAfter . Prev ;
            if ~ isempty ( nodeAfter . Prev )
                nodeAfter . Prev . Next = newNode ;
            end
            nodeAfter . Prev = newNode ;
        end
        function removeNode ( node )
            % Remove a node from a linked list .
            if ~ isscalar ( node )
                error ( ' Input must be scalar ')
                20
            end
            prevNode = node . Prev ;
            nextNode = node . Next ;
            if ~ isempty ( prevNode )
                prevNode . Next = nextNode ;
            end
            if ~ isempty ( nextNode )
                nextNode . Prev = prevNode ;
            end
            node . Next = dlnode . empty ;
            node . Prev = dlnode . empty ;
        end
        function clearList ( node )
            % Clear the list before
            % clearing list variable
            prev = node . Prev ;
            next = node . Next ;
            removeNode ( node )
            while ~ isempty ( next )
                node = next ;
                next = node . Next ;
                removeNode ( node ) ;
            end
            while ~ isempty ( prev )
                node = prev ;
                prev = node . Prev ;
                removeNode ( node )
            end
        end
    end
    methods ( Access = private )
        function delete ( node )
            clearList ( node )
        end
    end
end