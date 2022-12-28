function [lon, lat, alpha, delta] = computeGroundTrack(obj,theta_G_0,omega_planet,N_orbits,options)
    % [lon, lat, alpha, delta] = computeGroundTrack(obj,theta_G_0,omega_planet,N_orbits,options)
    %
    %  INPUT:
    %   - obj             :  orbit object
    %   - theta_G_0       :  initial Greenwich Sidereal Time
    %   - omega_planet    :  angular velocity of the planet
    %   - N_orbits        :  number of orbits to draw
    %   - options         :  struct containing all possible options:
    %                         + (ode_options) ode solver opts struct, define it with odeset(...)]
    %                         + (tspan) note: this can override the 'N_orbits' parameter!
    %                         + (propagationType) function to use for propagation (if not
    %                           inserted the default is 'odefun_2bp'):
    %                            - 'odefun_2bp' (no perturbation)
    %                            - 'odefun_2bp_perturbed' *
    %                            - 'motionEquation_RSW' * 
    %                         + (Color, LineWidth ...) plotting optional field
    %
    %  OUTPUT:
    %   - lon             :  longitude
    %   - lat             :  latitude
    %   - alpha           :  right ascension
    %   - delta           :  declination
    %
    %  NOTE:  the function automatically plots the track if no output is requested,
    %         if not, to plot the groundtrack please use plotGroundTrack
    %
    %   example of how to use the options field:
    %      options.planet = 'Earth';
    %      options.LineWidth = 1.1;
    %      options.Color = 'g';
    %      options.ode_options = odeset('RelTol', 1e-10);
    %      options.propagationType = 'motionEquation_RSW'
    %
    %   NOTE:
    %    the propagating functions that model the perturbations use only
    %    the J2 (gravitational) perturbation and the aerodynamic drag (the
    %    latter uses CIRA atmosphere, use help miaCIRA to know more)
    %
    %   * NOTE: for these propagation types the options field MUST contain
    %   also the variable 'param' needed to model the various perturbances.
    %   type 'help motionEquation_RSW' to know more about the param struct
    %
    % CHANGELOG:
    %  [Guglielmo Gomiero 26/11/2022]:
    %    - lat and lon are now rows (since time instants should advance on
    %    rows)
    %  [Carlo Zambaldo 24/12/2022]:
    %    - added the possibility to choose the perturbations when
    %      propagating the orbit
    %    - tidied up options field
    %    - added a more detailed 'help'
    %    - moved the plot to plotGroundTrack

    %% setup
    theta_G = @(t) theta_G_0 + omega_planet .* t;

    %% ode solver setup and integration
    if (nargin < 5) || isempty(options) || ~isfield(options,'tspan')
        tspan = linspace(0, N_orbits*obj.T, 1e5);
    else
        tspan = options.tspan;
    end

    if (nargin < 5) || isempty(options) || ~isfield(options,'ode_options')
        options.ode_options = odeset('RelTol', 1e-10);
    end

    % selecting the solver once the requested perturbation is known
    if (nargin < 5) || isempty(options) || ~isfield(options,'propagationType') || isequal(lower(options.propagationType),'odefun_2bp')
        param.mu = obj.mu;
        y_0 = [obj.r0; obj.v0]; % state vector for cartesian approach
        [t_vect,y] = ode113( @(t,y) odefun_2bp(t, y, param), tspan, y_0, options.ode_options);
    elseif isfield(options,'param')
        param = options.param;
        if isequal(lower(options.propagationType),'odefun_2bp_perturbed')
            y_0 = [obj.r0; obj.v0]; % state vector for cartesian approach
            [t_vect,y] = ode113( @(t,y) odefun_2bp_perturbed(t, y, param), tspan, y_0, options.ode_options);
        elseif isequal(lower(options.propagationType),'motionequation_rsw')
            warning("This method is not efficient. Consider using propagationType = 'odefun_2bp_perturbed'.");
            y_0 = [obj.a; obj.e; obj.i; obj.O; obj.w; obj.theta]; % state vector for Gauss' method
            [t_vect,y] = ode113( @(t,y) motionEquation_RSW(t,state, @(t,y)acc_pert_fun_rsw(t,y,param),param), tspan, y_0, options.ode_options);
            for hh = 1:length(t_vect)   % translate the state vector in cartesian coordinates
                [rr, vv] = ORBIT.kep2car(y(hh,1),y(hh,2),y(hh,3),y(hh,4),y(hh,5),y(hh,6),param.mu);
                y(hh,:) = [rr(:)', vv(:)'];
            end
        else
            error("An error has occurred.");
        end
    else
        error("'param' field is not defined. Please define it and feed it into the groundTrack function through the options struct.")
    end

    %%
    r_vect_prop = y(:,1:3)'; % extract the position vectors for all the times
    rx = r_vect_prop(1,:);
    ry = r_vect_prop(2,:);
    rz = r_vect_prop(3,:);
    r_norm = sqrt(sum(r_vect_prop.^2,1)); % take the norm

    % angles
    delta = asin(rz./r_norm);
    alpha = atan2(ry,rx);

    % from here angles in degrees
    lon = rad2deg(wrapToPi((alpha-theta_G(t_vect'))));
    lat = rad2deg(wrapToPi(delta));


    %% PLOT
    if nargout == 0
        if (nargin < 5) || isempty(options)
            obj.plotGroundTrack(lon, lat);
        else
            obj.plotGroundTrack(lon, lat, options);
        end
    end
end






