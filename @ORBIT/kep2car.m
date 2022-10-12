function [r_vect, v_vect] = kep2car(ORBIT_OBJ)
    % [r_vect, v_vect] = kep2car(ORBIT_OBJ)
    %
    % the function outputs the corresponding cartesian coordinates of a
    % given orbit
    %
    %  INPUT:
    %     - ORBIT_OBJ   :   orbit object, defined with ORBIT
    %
    %  OUTPUT:
    %     - r_vect      :   vector containing the x,y,z components of position [km]
    %     - v_vect      :   vector containing the x,y,z components of the velocity [km/s]


    %% initialisation
    a = ORBIT_OBJ.a;
    e = ORBIT_OBJ.e;
    i = ORBIT_OBJ.i;
    O = ORBIT_OBJ.O;
    w = ORBIT_OBJ.w;
    theta = ORBIT_OBJ.theta;
    mu = ORBIT_OBJ.mu;

    % rotation matrix between perifocal to geocentric frame
    R_PF_GE = [ cos(w)*cos(O)-sin(w)*cos(i)*sin(O),	-sin(O)*cos(i)*cos(w)-cos(O)*sin(w),	sin(i)*sin(O);
                cos(O)*cos(i)*sin(w)+sin(O)*cos(w), cos(O)*cos(i)*cos(w)-sin(O)*sin(w),    -cos(O)*sin(i);
                             sin(i)*sin(w),                        sin(i)*cos(w),                cos(i)  ];
    
    % parameters
    p = a*(1-e^2);          % orbital parameter
    r = p/(1+e*cos(theta)); % radious

    %% r_vect and v_vect (in the PERIFOCAL FRAME)
    r_vect = r * [cos(theta); sin(theta); 0];
    v_vect = sqrt(mu/p) * [-sin(theta); e+cos(theta); 0];

    %% transform r_vect and v_vect in the GEOCENTRIC FRAME
    r_vect = R_PF_GE * r_vect;
    v_vect = R_PF_GE * v_vect;
end