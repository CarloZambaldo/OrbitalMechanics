classdef ORBIT
    % ORBIT is an object describing an orbit, the principal fields are:
    %    a        : semi-major axis [km]
    %    e        : eccentricity [-]
    %    i        : inclination [rad]
    %    O        : RAAN [rad]
    %    w        : anomaly of pericentre [rad]
    %    theta    : true anomaly (usually initial position) [rad]
    %    mu       : gravitational constant, for Earth use astroConstants(13)
    % 
    % in addition, the by defining the orbit using ORBIT(...) the following
    % fields are filled in:
    %    r0        : position on orbit at given theta
    %    v0        : velocity on orbit at given theta
    %    T         : period of the orbit [s]
    %    energy    : energy of the orbit
    %    p         : orbital parameter
    %    h_vect    : angular momentum vector (per unit mass)
    %    e_vect    : eccentricity vector
    %    rA        : radius of apoapsis
    %    rP        : radius of periapsis
    %    theta_inf : true anomaly asymptote
    %    delta     : deflection
    %    v_inf     : hyperbolic eccess speed
    %    type      : char [+, 0, -] respectively [prograde, polar, retrograde]
    %
    %
    % the orbit can be defined in parametrical form by inserting all 7
    % parameters like: obj = ORBIT(a,e,i,O,w,theta,mu)
    % if the orbit is given in Cartesian form use the function keplerian to
    % define the orbit.  WARNING: r_vect,v_vect and h are initialized to 0 vectors!
    %
    % - to define the orbit in degrees use ORBIT.define_in_deg(....) 
    %
    %
    % useful methods:
    %  - obj = ORBIT(a,e,i,O,w,theta,mu)
    %  - obj = define_in_deg(a,e,i_deg,O_deg,w_deg,theta_deg,mu)
    %  - [r_vect, v_vect] = cartesian(ORBIT_OBJ);
    %  - [ORBIT_OBJ] = keplerian(r_vect, v_vect, mu);
    %  - [r_prop_vect] = plotOrbit(ORBIT_OBJ,theta_i,theta_f,d_theta,grafica,theta_segnato);
    %  - [alpha, delta, lon, lat] = groundTrack(ORBIT_OBJ,theta_G_0,omega_planet,N,plot_options);
    %
    %
    % ------------------------------------------
    %   Carlo Zambaldo (info@carlozambaldo.it)
    %           last update: 30.11.22
    % ------------------------------------------
    %
    %   CHANGELOG:
    % [Guglielmo Gomiero 26/11/2022]:
    %   - various comments and typos
    %   - correct uses of obj variable and function calls
    %   - obj.p now is set correctly, after we made sure a for hyperbolas is <0
    %   - since we use them more than once and in groundTrack, added r0 and v0 as properties
    %   - modified h_vect and e_vect to account for this
    % [Carlo Zambaldo 30/11/2022]
    %   - properties "grafically" rearranged
    %   - fixed bug for obj.p in parabolas (p was overwritten during def)
    %   - to define r0 and v0 with kep2car use obj propreties to avoid possible bugs
    %   - added cartesian and keplerian methods (object-oriented kep2car
    %   and car2kep counterparts)
    % [Carlo Zambaldo 06/12/2022]
    %   - fixed bug: for hyperbolas, if a>0 then obj.a is set to negative

    properties
        a        {mustBeNumeric}                        % [km]
        e        {mustBeNumeric}                        % [-]
        i        {mustBeNumeric,mustBeNonnegative}      % [rad]
        O        {mustBeNumeric}                        % [rad]
        w        {mustBeNumeric}                        % [rad]
        theta    {mustBeNumeric}                        % [rad]
        mu       {mustBeNumeric}                        % [km^3/s^2]
        
        r0       (3,1){mustBeNumeric}
        v0       (3,1){mustBeNumeric}

        T        {mustBeNumeric}                        % [s]
        energy   {mustBeNumeric}
        p        {mustBeNumeric, mustBeNonnegative}     % orbital parameter
        h_vect   (3,1) {mustBeVector}
        e_vect   (3,1) {mustBeVector}

        rA       {mustBeNumeric}                        % radius of APOCENTRE
        rP       {mustBeNumeric}                        % radius of PERICENTRE

        theta_inf{mustBeNumeric}                        % true anomaly asymptote
        delta    {mustBeNumeric}                        % deflection
        v_inf    {mustBeNumeric}

        type     % ceph (circular, elliptical, parabolic, hyperbolic); +=- (prograde,polar,retrograde)
    end
    
    %% METHODS
    % these can be called simply as ORBITobject.methodname
    methods
        % Definition of ORBIT
        function obj = ORBIT(a,e,i,O,w,theta,mu)
            if nargin~=0
                obj.e = e;
                obj.i = wrapToPi(i);    % [0,pi]
                obj.O = wrapTo2Pi(O);   % [0,2pi]
                obj.w = wrapTo2Pi(w);   % [0,2pi]
                obj.theta = theta;
                obj.mu = mu;

                if e == 1        %%% PARABOLA    
                    obj.T = inf; % Infinity if parabola

                    if obj.a ~= Inf
                        fprintf("Interpreting the given value ""a"" as ""rP"" to fully define parabolic trajectory.\n");
                        obj.rP = a;
                        obj.a = Inf;
                        obj.p = 2*obj.rP; % define p
                    else
                        error("One value is missing to fully describe parabolic trajectory.")
                    end
                    
                else
                    if e>1             %%% HYPERBOLA    
                        if a>0
                            a = -a; warning("The semi-major axis has been set to negative.");
                        end
                        obj.T = NaN; % Not defined if hyperbolic trajectory
                        obj.theta_inf = acos(1./(-e)); % true anomaly asymptote
                        obj.delta = 2*asin(1/e);       % deflection
    
                    else                 %%% ELLIPSE/CIRCLE
                        obj.T = 2*pi*sqrt(a^3/mu);
                    end
                    % define here p to avoid replacing p of parabola
                    obj.p = a*(1-e^2);
                end
                obj.a = a;
                obj.energy = -mu/(2*a);

                % defining extra values (filling the "gaps")
                [obj.r0,obj.v0] = ORBIT.kep2car(obj.a,obj.e,obj.i,obj.O,obj.w,obj.theta,obj.mu);
                obj.h_vect = cross(obj.r0,obj.v0);
                obj.e_vect = (cross(obj.v0, obj.h_vect)./obj.mu) - (obj.r0 ./norm(obj.v0));
                obj.type = obj.typedef;

                % Compute other parameters:
                if e ~= 1 % Only for non-parabolic trajectories
                    obj.rP = a*(1-e);
                end
                obj.rA = a*(1+e);
                obj.v_inf = sqrt(mu/abs(a));
            end
        end

        % Keplerian angles in degrees
        function [i_deg,O_deg,w_deg,theta_deg] = get_angles_deg(obj)
            i_deg = rad2deg(obj.i);
            O_deg = rad2deg(obj.O);
            w_deg = rad2deg(obj.w);
            theta_deg = rad2deg(obj.theta);
        end

      % Type of orbit
        function type = typedef(obj)
            type = 'ERROR';

            if obj.h_vect(3)>0 && obj.i>=0 && obj.i<=pi/2
                type = '+'; % prograde orbit
            elseif obj.h_vect(3)==0 % do not compare i==pi/2 due to possible numerical error
                type = '='; % polar orbit
            elseif obj.h_vect(3)<0 && obj.i>=pi/2 && obj.i<=pi
                type = '-'; % retrograde orbit
            end

            if     obj.e==0
                type = strcat('c',type);
            elseif obj.e<1
                type = strcat('e',type);
            elseif obj.e==1
                type = strcat('p',type);
            elseif obj.e>1
                type = strcat('h',type);
            end
        end
        
        
        %%% Methods in other files
        % Plot orbit in a geocentric, equatorial frame
        [r_prop_vect] = plotOrbit(obj,theta_vect,d_theta,grafica,theta_segnato);

        % Propagate cartesian state with integration
        [t,y] = propagateOrbit(obj,y_0,time_vect,N_dt,grafica,r_segnato);

        % Plot ground track
        [alpha, delta, lon, lat] = groundTrack(obj,theta_G_0,omega_planet,N_orbits,plot_options);

        % Time of flight
        delta_time = TOF(orbit_obj, thetas);
    end

    %% static methods
    % These might be useful even when I don't want to declare a full ORBIT object
    % Must be called like ORBIT_object.methodname(params).
    % Do not require a ORBIT_object as input
    methods(Static)
        % allows to define orbit in degrees
        function obj = define_in_deg(a,e,i_deg,O_deg,w_deg,theta_deg,mu)
            i_rad = deg2rad(i_deg);
            O_rad = deg2rad(O_deg);
            w_rad = deg2rad(w_deg);
            theta_rad = deg2rad(theta_deg(:));
            obj = ORBIT(a,e,i_rad,O_rad,w_rad,theta_rad,mu);
        end
        
        %%% Static methods in other files
         [r_vect, v_vect] = kep2car(a,e,i,O,w,theta,mu);
         [a,e,i,O,w,theta] = car2kep(r_vect, v_vect, mu);
         ground_track(lat,lon,opts,planet);

        [r_vect, v_vect] = cartesian(orbit_obj);
        [orbit_obj]      = keplerian(r_vect, v_vect, mu);
      end
end