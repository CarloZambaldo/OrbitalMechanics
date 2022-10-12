classdef ORBIT
    % the orbit can be defined in Parametrical form by inserting all 7 parameters
    % if the orbit is given in Cartesian form just use the function
    % car2par. WARNING: r0,v0 and h are initialized to 0 vectors!
    
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
        
        T {mustBeNonnegative, mustBeNumeric} % [s]
        energy {mustBeNumeric}
        p {mustBeNumeric}         % orbital parameter
        h_vect (3,1) {mustBeVector}
        e_vect (3,1) {mustBeVector}
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
                obj.energy = -mu/(2*a);
                obj.T = 2*pi*sqrt(a^3/mu);
                obj.p = a*(1-e^2); 
                
                if isequal(obj.r0(:),zeros(3,1)) % if not defined
                    [r0,v0] = kep2car(obj);
                    obj.r0 = r0(:);
                    obj.v0 = v0(:);
                    obj.h_vect = cross(r0, v0);
                    obj.e_vect = (cross(v0, h_vect)./mu) - (r0./norm(r0));
                end
            end
        end

        %% functions
        [r_vect, v_vect] = kep2car(ORBIT_OBJ);
  
        function [i_deg,O_deg,w_deg,theta_deg] = orbitAngles2deg(ORBIT_OBJ)
            i_deg = ORBIT_OBJ.i * 180/pi;
            O_deg = ORBIT_OBJ.O * 180/pi;
            w_deg = ORBIT_OBJ.w * 180/pi;
            theta_deg = ORBIT_OBJ.theta * 180/pi;
        end


    end
end