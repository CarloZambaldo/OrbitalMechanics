classdef ORBIT
    % the orbit can be defined in Parametrical form by inserting all 7 parameters
    % if the orbit is given in Cartesian form just use the function
    % car2par.
    
    properties
        a  {mustBeNumeric}      % [km]
        e  {mustBeNumeric}      % [-]
        i  {mustBeNumeric}      % [rad]
        O  {mustBeNumeric}      % [rad]
        w  {mustBeNumeric}      % [rad]
        theta {mustBeNumeric}  % [rad]
        mu {mustBeNumeric}      % [km^3/s^2]
        r0 (3,1) {mustBeVector}       % [km]
        v0 (3,1) {mustBeVector}       % [km]
        energy {mustBeNumeric} 
        h (3,1) {mustBeVector}
        type                    % +0- (prograde,polar,retrograde)
    end
    
    methods
        %% definition of ORBIT
        function obj = ORBIT(a,e,i,O,w,theta,mu)
            if nargin~=0
                obj.a = a;
                obj.e = e;
                obj.i = i;
                obj.O = O;
                obj.w = w;
                obj.theta = theta;
                obj.mu = mu;

                [r0,v0] = par2car(a,e,i,O,w,theta,mu);
                obj.r0 = r0;
                obj.v0 = v0;
            end
        end
    
        %% functions
        [ORBIT_OBJ] = car2kep(r0, v0);
        [r_vect, v_vect] = kep2car(ORBIT_OBJ);

        [alpha, delta, lon, lat] = groundTrack(ORBIT_OBJ,theta_G_0,omega_planet,N);

        
        function [i_deg,O_deg,w_deg,theta_deg] = orbitAngles2deg(ORBIT_OBJ)
            i_deg = ORBIT_OBJ.i * 180/pi;
            O_deg = ORBIT_OBJ.O * 180/pi;
            w_deg = ORBIT_OBJ.w * 180/pi;
            theta_deg = ORBIT_OBJ.theta * 180/pi;
        end
    end
end