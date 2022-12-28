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
    %    type      : char [c,e,p,h] circular, elliptical, parabolic, hyperbolic
    %                 and [+,=,-] respectively prograde, polar, retrograde
    %
    %
    % ------------------ brief explanation of the class ------------------
    %  DEFINITION OF AN ORBIT OBJECT
    %    the orbit can be defined in parametrical form by inserting all 7
    %    parameters like: obj = ORBIT(a,e,i,O,w,theta,mu)
    %    if the orbit is given in cartesian form use the function keplerian to
    %    define the orbit like: obj = ORBIT.keplerian(r0, v0, mu).
    %    To define the orbit in degrees use ORBIT.define_in_deg(a,e,i,O,w,theta,mu).
    %    WARNING: r_vect, v_vect and h are initialized as 0 vectors, if
    %    there is a problem in the script these values are likely to remain
    %    zero. Moreover, check the 'type' proprety, in these cases it usually
    %    displays an 'error' message.
    %
    %  DEFINITION OF A SIMPLIFIED ORBIT OBJECT
    %    it is also possible to only store a,e,i,O,w,theta,mu values by
    %    calling the obj = simpleOrbit(a,e,i,O,w,theta,mu) method. If this type
    %    of orbit has to be filled, in a second moment, with all the other
    %    parameters: use obj.fillOrbit.
    %    
    %  PLOTTING OF AN ORBIT
    %    to plot a 3D orbit use the method obj.plotOrbit [if needed add the
    %    optional parameters], to propagate an orbit use instead
    %    obj.propagateOrbit, this last method uses ode113 to propagate the
    %    orbit given a tspan, it is possible to modify the propagation
    %    method (i.e. if perturbances have to be taken into account) using 
    %    option.propagationType. Write help propagateOrbit to know more.
    %
    %  PLOTTING OF GROUND TRACKS
    %    ground-tracks plotting can simply be performed by using
    %    computeGroundTrack method without requiring any output. If
    %    computeGroundTrack is called requiring the outputs use
    %    plotGroundTrack to also plot the required groundTrack.
    %     note: groundTrack function is deprecated and in a future release
    %     will no longer work.
    %
    %   CONVERT BETWEEN CARTESIAN AND KEPLERIAN ELEMENTS AND VICE VERSA
    %    two possible approaches are available: if the orbit is not
    %    hyperbolic nor parabolic the user is highly encouraged to use
    %    kep2car and car2kep methods. The second approach is
    %    object-oriented and requires the definition of an ORBIT object to
    %    work: use therefore cartesian or keplerian methods to convert
    %    between the two rapresentations.
    %
    %   EXTRAS
    %     to compute the time of flight between two points use TOF()
    %     to check if a position vector is on the orbit use isOnOrbit
    % 
    %
    % ------------------ list of methods ------------------
    % methods:
    %  - [t,y] = propagateOrbit(obj,tspan,options,grafica,r_segnato);
    %  - [r_prop_vect] = plotOrbit(obj,theta_vect,d_theta,grafica,theta_segnato);
    %  - [lon, lat, alpha, delta] = computeGroundTrack(obj,theta_G_0,omega_planet,N_orbits,options)
    %  - type = typedef(obj)
    %  - delta_time = TOF(obj, thetas);
    %  - delta_theta = AOF(obj, delta_time, toll)
    %  - v_esc = escapeVelocity(obj, position)
    %  - check = isOnOrbit(obj, r_vect)
    %  - obj = fillOrbit(obj)
    %
    % static methods:
    %  - obj = define_in_deg(a,e,i_deg,O_deg,w_deg,theta_deg,mu)
    %  - obj = simpleOrbit(a,e,i,O,w,theta,mu)
    %  - [r_vect, v_vect]  = kep2car(a,e,i,O,w,theta,mu);
    %  - [a,e,i,O,w,theta] = car2kep(r_vect, v_vect, mu);
    %  - [r_vect, v_vect] = cartesian(obj, theta);
    %  - [obj] = keplerian(r_vect, v_vect, mu);
    %  - [] = plotGroundTrack(lon, lat, options)
    %
    %
    %
    % EXTREME WARNING
    %   THE ORBIT OBJECT TO WORK PROPERLY HAS TO BE IN A FOLDER NAMED
    %   PRECISELY "@ORBIT" WITH ALL ITS METHODS, THIS FOLDER HAS INDEED 
    %   TO BE ADDED TO THE PATH. NO EXCEPTIONS CAN BE MADE.
    %
    % ------------------------------------------
    %   Carlo Zambaldo (info@carlozambaldo.it)
    %         first version: 20.10.2022
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
    % [Carlo Zambaldo 09/12/2022]
    %   - fixed bug for parabolas: "if a ~= Inf" instead of if "obj.a ~= Inf"
    %   - added "isOnOrbit" and "escapeVelocity" methods
    %   - fixed bug for parabolas: use cartesian method instead of kep2car
    %     to define r0 and v0 (this allows more control on orbit type.)
    %   - added "AOF" (angle of flight), which is the counterpart of TOF,
    %     this function allows to compute the angle between two times
    % [Carlo Zambaldo 14/12/2022]
    %   - added "simpleOrbit" and "fillOrbit" methods
    % [Carlo Zambaldo 28/12/2022]
    %   - written a more accurate help

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

        theta_inf{mustBeNumeric}                        % true anomaly asymptote (only for hyperbola)
        delta    {mustBeNumeric}                        % deflection             (only for hyperbola)
        v_inf    {mustBeNumeric}                        % asymptotic speed       (only for hyperbola)

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

                    if a ~= Inf
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
                %[obj.r0,obj.v0] = ORBIT.kep2car(obj.a,obj.e,obj.i,obj.O,obj.w,obj.theta,obj.mu);
                [obj.r0,obj.v0] = obj.cartesian(obj,theta); % use this (or change kep2car to account for a=Inf for parabolas)
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


        % Type of orbit
        function type = typedef(obj)
            type = 'ERROR';

            % this also allows to check if there is any inconsistency between h_vect and i
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

        % Plot ground track
        [lon, lat, alpha, delta] = computeGroundTrack(obj,theta_G_0,omega_planet,N_orbits,options)

        % to maintain the backward compatibility
        function [alpha, delta, lon, lat] = groundTrack(obj,theta_G_0,omega_planet,N_orbits,plot_options)
            warning("This function is deprecated. Please use computeGroundTrack and plotGroundTrack instead.");
            [lon, lat, alpha, delta] = computeGroundTrack(obj,theta_G_0,omega_planet,N_orbits,options)
        end

        % Time of flight
        delta_time = TOF(obj, thetas);

        % delta_theta between delta_t
        delta_theta = AOF(obj, delta_time, toll)

        % escape velocity of given point
        function v_esc = escapeVelocity(obj, position)
            % computes v_esc for given position on orbit object, if
            % position is [1x1] the function interprets it as an angle, if
            % position is [3x1] or [1x3] it is interpreted as a radius
            position = position(:);
            if size(position,1) == 1 % theta
                r = (norm(obj.h_vect)^2/obj.mu)/(1+obj.e * cos(position));
            elseif size(position,1) == 3 % r_vect
                if obj.isOnOrbit(position)
                    r = norm(position);
                end
            else
                error("Unknown position on orbit");
            end
            v_esc = sqrt(2*obj.mu/norm(r));
        end

        % isOnOrbit
        function check = isOnOrbit(obj, r_vect)
            % is on orbit checks if the given radius is on the obj
            toll = 1e-13;
            r_vect = r_vect(:);
            r = norm(r_vect);
            ang = acos((obj.h^2/(obj.mu*r)-1)/obj.e);
            if ~(dot(obj.h_vect,r_vect)<=toll) % first check if radius is on the right plane
                check = 0;
            elseif ~isreal(ang) % then check if angle is real 
                check = 0;
            else
                % at this stage checks if the vector is equal to the one
                % expected at that given ang
                [r_expected, ~] = obj.cartesian(obj,ang);
                if norm(r_expected-r_vect) <= toll
                    % oh, dear! that's the right vector. thanks God.
                    check = 1;
                else
                    check = 0;
                end
            end
        end
    
        % fills all available values for given obj
        function obj = fillOrbit(obj)
            obj = ORBIT(obj.a,obj.e,obj.i,obj.O,obj.w,obj.theta,obj.mu);
        end

        % Propagate cartesian state with integration
        [t,y] = propagateOrbit(obj,tspan,options,grafica,r_segnato);
    
    end

    %% static methods
    % These might be useful even when I don't want to declare a full ORBIT object
    % Must be called like ORBIT_object.methodname(params).
    % Do not require a ORBIT_object as input
    methods(Static)
        
        function obj = define_in_deg(a,e,i_deg,O_deg,w_deg,theta_deg,mu)
            % allows to define orbit in degrees
            i_rad = deg2rad(i_deg);
            O_rad = deg2rad(O_deg);
            w_rad = deg2rad(w_deg);
            theta_rad = deg2rad(theta_deg(:));
            obj = ORBIT(a,e,i_rad,O_rad,w_rad,theta_rad,mu);
        end

        function obj = simpleOrbit(a,e,i,O,w,theta,mu)
            % uses obj ORBIT to store the geometrical values that describe the
            % orbit (useful to reduce the computation time i.e. for iterative
            % plotting procedure), to fill all the missing values call
            % obj.fillOrbit method.
            % NOTE: this method does NOT verify the type of the orbit.
            
                obj = ORBIT;
                obj.a = a;
                obj.e = e;
                obj.i = i;
                obj.O = O;
                obj.w = w;
                obj.theta = theta;
                obj.mu = mu;
                obj.type = 'simplified';
            if e >= 1
                if e == 1
                    warning("The defined simplified orbit is a parabola, " + ...
                        "there could be some problems with the definition of ""a"" " + ...
                        "and ""p"". Please use ORBIT() method to define a consistend parabolic orbit.");
                else
                    warning("The defined simplified orbit is an hyperbola.");
                end
            else
                obj.p = a*(1-e^2);
            end
        end

        %%% Static methods in other files
        [r_vect, v_vect]  = kep2car(a,e,i,O,w,theta,mu);
        [a,e,i,O,w,theta] = car2kep(r_vect, v_vect, mu);

        [r_vect, v_vect] = cartesian(obj, theta);
        [obj] = keplerian(r_vect, v_vect, mu);

        % plot ground track
        [] = plotGroundTrack(lon, lat, options)
      end
end