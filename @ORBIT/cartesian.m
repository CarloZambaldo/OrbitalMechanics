function [r_vect, v_vect] = cartesian(orbit_obj, theta)
    % [r_vect, v_vect] = cartesian(orbit_obj)
    %
    % the function outputs the corresponding cartesian coordinates of a
    % given orbit
    %
    %  INPUT:
    %     - orbit_obj   :   orbit object, defined with ORBIT
    %     - theta       :   (optional) true anomaly -> if not inserted uses orbit_obj.theta
    %     
    %
    %  OUTPUT:
    %     - r_vect      :   vector containing the x,y,z components of position [km]
    %     - v_vect      :   vector containing the x,y,z components of the velocity [km/s]
    %
    %  note: if nargoud == 1 the function outputs the state vector:
    %     - y_vect      :   state vector such that y_vect = [r_vect(:); v_vect(:)]


    %% initialisation
    a = orbit_obj.a;
    e = orbit_obj.e;
    i = orbit_obj.i;
    O = orbit_obj.O;
    w = orbit_obj.w;
    mu = orbit_obj.mu;

    if nargin<2
        theta = orbit_obj.theta;
    end

    % rotation matrix between perifocal to geocentric frame
    R_PF_GE = [ cos(w)*cos(O)-sin(w)*cos(i)*sin(O),	-sin(O)*cos(i)*cos(w)-cos(O)*sin(w),	sin(i)*sin(O);
                cos(O)*cos(i)*sin(w)+sin(O)*cos(w), cos(O)*cos(i)*cos(w)-sin(O)*sin(w),    -cos(O)*sin(i);
                             sin(i)*sin(w),                        sin(i)*cos(w),                cos(i)  ];
    
    % parameters
    p = orbit_obj.p;        % orbital parameter
    r = p/(1+e*cos(theta)); % radius

    %% r_vect and v_vect (in the PERIFOCAL FRAME)
    r_vect = r * [cos(theta); sin(theta); 0];
    v_vect = sqrt(mu/abs(p)) * [-sin(theta); e+cos(theta); 0];

    %% transform r_vect and v_vect in the GEOCENTRIC FRAME
    r_vect = R_PF_GE * r_vect;
    v_vect = R_PF_GE * v_vect;

    if isequal(r_vect(:), zeros(3,1)) || isequal(v_vect(:), zeros(3,1))
        warning("There could be a problem in the definition of radius and velocity vectors.");
    end

    if nargout == 1 % build state vector
        r_vect = [r_vect(:); v_vect(:)];
    end
end