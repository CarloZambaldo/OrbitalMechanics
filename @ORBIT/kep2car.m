function [r_vect, v_vect] = kep2car(a,e,i,O,w,theta,mu)
    % [r_vect, v_vect] = kep2car(a,e,i,O,w,theta,mu)
    %
    % the function outputs the corresponding cartesian coordinates of a
    % given orbit
    %
    %  INPUT:
    %     a        : semi-major axis [km]
    %     e        : eccentricity [-]
    %     i        : inclination [rad]
    %     O        : RAAN [rad]
    %     w        : anomaly of pericentre [rad]
    %     theta    : true anomaly (usually initial position) [rad]
    %     mu       : gravitational constant, for Earth use astroConstants(13)
    %     
    %
    %  OUTPUT:
    %     - r_vect      :   vector containing the x,y,z components of position [km]
    %     - v_vect      :   vector containing the x,y,z components of the velocity [km/s]
    %
    %  note: if nargoud == 1 the function outputs the state vector:
    %     - y_vect      :   state vector such that y_vect = [r_vect(:); v_vect(:)]
    % CHANGELOG:
    % [Guglielmo Gomiero 26/11/2022]:
    % -inputs and outputs set to be arrays, do not use objects

    %% initialisation
    % Ensure theta is row (needed for following matrices)
    theta = theta(:)';

    % rotation matrix between perifocal to geocentric frame
    R_PF_GE = [ cos(w)*cos(O)-sin(w)*cos(i)*sin(O),	-sin(O)*cos(i)*cos(w)-cos(O)*sin(w),	sin(i)*sin(O);
                cos(O)*cos(i)*sin(w)+sin(O)*cos(w), cos(O)*cos(i)*cos(w)-sin(O)*sin(w),    -cos(O)*sin(i);
                             sin(i)*sin(w),                        sin(i)*cos(w),                cos(i)  ];
    
    % parameters
    p = a*(1-e^2);        % orbital parameter
    r = p./(1+e.*cos(theta)); % radius

    %% r_vect and v_vect (in the PERIFOCAL FRAME)
    r_vect = r .* [cos(theta); sin(theta); zeros(1,length(theta))];
    v_vect = sqrt(mu/abs(p)) .* [-sin(theta); e+cos(theta); zeros(1,length(theta))];

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