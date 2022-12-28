function [t,y] = propagateOrbit(obj,tspan,options,grafica,r_segnato)
    %    [t,y] = propagateOrbit(obj,tspan,options,grafica,r_segnato)
    %
    % plots the orbit, given the parameters ORBIT_OBJ
    % using as initial and final times tspan(1) and tspan(end).
    %
    % INPUT:
    %   - mu              :     
    %   - y0[6x1]         :    initial state, equatorial cartesian coordinates
    %   - tspan           :    define it as [ti tf] or [ti,t2,t3,...,tf]
    %                          time vector
    %   - options         :  struct containing all possible options:
    %                        + (ode_options) ode solver opts struct, define it with odeset(...)]
    %                        + (propagationType) function to use for propagation (if not
    %                          inserted the default is 'odefun_2bp'):
    %                           - 'odefun_2bp' (no perturbation)
    %                           - 'odefun_2bp_perturbed' *
    %                           - 'motionEquation_RSW' *
    %                        + (mu) Note: this and other parameters are
    %                        required if a different propagationType from
    %                        'odefun_2bp' is used!! 
    %   - grafica         :    the graphical style of the function "plot" used to plot
    %                          the orbit
    %   - r_segnato[3x1]  :    optional: the function will add a star in the position
    %                          corresponding given by r_segnato
    %
    % OUTPUT
    %   - t[Nx1]          :     vector of times
    %   - y[Nx6]          :     output of ode integration WARNING! the
    %                           output can either be a set of keplerian
    %                           orbital parameters or the vector [r;v]
    %                           depending on the method of integration
    %                           defined by the parameter propagationType
    %
    % LOGCHANGE:
    % [Guglielmo Gomiero, 26/11/2022]:
    %   - language fully set to english
    %   - added input description
    % [Guglielmo Gomiero, 27/11/2022]:
    %   - commented out plotting
    %   - commented out third party Earth plotting function, as it causes problems
    %     with the variable y(declared 2 times, y is also the output)
    % [Carlo Zambaldo, 06/12/2022]:
    %   - added plotting only if nargout <= 1, fixed "y" bug
    % [Carlo Zambaldo, 07/12/2022]:
    %   - removed planet plotting: use instead ASTRO to define an ASTRO object
    %     and then call ASTRO.plotAstro method to draw the celestial body
    % [Carlo Zambaldo, 12/12/2022]:
    %   - changed the prototype, now the propagateOrbit call is simplier and more intuitive:
    %     it is similar to calling an ode solver.
    %   - time_vect is now deprecated. please use tspan to define (if
    %     needed) a time span where to compute the function
    %     i.e. tspan = [1:10] is equivalent to the old time_vect = [1 10],
    %     N_dt = 10
    % [Carlo Zambaldo, 24/12/2022]:
    %   - added 'options' parameter
    %   - unfortunately to optimise the funciton the prototype has been drastically
    %     changed, this may cause an incompatibility with older code.
    %   - implemented different perturbations


    %% ode solver setup and integration
    if (nargin < 3) || isempty(options) || ~isfield(options,'ode_options')
        ode_options = odeset('RelTol', 1e-13);
    else
        ode_options = options.ode_options;
    end

    if (nargin < 3) || isempty(options) || ~isfield(options,'param')
        param.mu = obj.mu;
    else
        param = options.param;
    end
    
    % selecting the solver once the requested perturbation is known
    if (nargin < 3) || isempty(options) || ~isfield(options,'propagationType') || isequal(lower(options.propagationType),'odefun_2bp')
        y_0 = [obj.r0; obj.v0]; % state vector for cartesian approach
        [t,y] = ode113( @(t,y) odefun_2bp(t, y, param), tspan, y_0, ode_options);
    elseif isfield(options,'param')
        param = options.param;
        if isequal(lower(options.propagationType),'odefun_2bp_perturbed')
            y_0 = [obj.r0; obj.v0]; % state vector for cartesian approach
            [t,y] = ode113( @(t,y) odefun_2bp_perturbed(t, y, param), tspan, y_0, ode_options);
        elseif isequal(lower(options.propagationType),'motionequation_rsw')
            y_0 = [obj.a; obj.e; obj.i; obj.O; obj.w; obj.theta]; % state vector for Gauss' method
            [t,y] = ode113( @(t,y) motionEquation_RSW(t,y, @(t,y)acc_pert_fun_rsw(t,y,param),param), tspan, y_0, ode_options);
        elseif isequal(lower(options.propagationType),'motionequation_tnh')
            y_0 = [obj.a; obj.e; obj.i; obj.O; obj.w; obj.theta]; % state vector for Gauss' method
            [t,y] = ode113( @(t,y) motionEquation_TNH(t,y, @(t,y)acc_pert_fun_tnh(t,y,param),param), tspan, y_0, ode_options);
        else
            error("An error has occurred. Could not propagate the orbit.");
        end
    else
        error("'param' field is not defined. Please define it and feed it into the groundTrack function through the options struct.")
    end

    %% disegno
    if nargout == 0 % plot the propagation if no output is required.
        if nargin < 4 % tratto del disegno impostato di default
            grafica = '-';
        end

        if ~isequal(lower(options.propagationType),'odefun_2bp') && ~isequal(lower(options.propagationType),'odefun_2bp_perturbed')
            warning("This method is not efficient for only plotting. Consider using propagationType = 'odefun_2bp_perturbed'.");
            for hh = 1:length(t)   % translate the state vector in cartesian coordinates
                [rr, vv] = ORBIT.kep2car(y(hh,1),y(hh,2),y(hh,3),y(hh,4),y(hh,5),y(hh,6),param.mu);
                y(hh,:) = [rr(:)', vv(:)'];
            end
        end

        hold on;
        plot3(y(:,1), y(:,2), y(:,3), grafica, 'LineWidth',1.3);
        axis image; 

        if nargin == 5 % significa che si è inserito anche r_segnato (metto un asterisco in r_segnato)
            hold on;
            plot3(r_segnato(1), r_segnato(2), r_segnato(3), 'r*', 'LineWidth',1.3);
        end
    end
end
