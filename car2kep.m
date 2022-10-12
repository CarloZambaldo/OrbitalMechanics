function [ORBIT_OBJ] = car2kep(r_vect, v_vect, mu)
    % [ORBIT_OBJ] = car2kep(r_vect, v_vect, mu)
    %
    % the function outputs the corresponding Keplerian coordinates of a
    % given orbit
    %
    %  INPUT:
    %     - r_vect      :   vector containing the x,y,z components of position [km]
    %     - v_vect      :   vector containing the x,y,z components of the velocity [km/s]
    %     - mu          :   planetary constant
    %
    %  OUTPUT:
    %     - ORBIT_OBJ   :   orbit object, defined with ORBIT

    % to be sure r and v are column vectors
    r_vect = r_vect(:);
    v_vect = v_vect(:);

    % calculate the norm of r and v
    r_norm = norm(r_vect);
    v_norm = norm(v_vect);

    % semimajor axis (from energy equation) [km]
    a = 1/(2/r_norm - v_norm^2/mu);
    
    % angular momentum 
    h_vect = cross(r_vect, v_vect);
    h_norm = norm(h_vect);

    % inclination [rad]
    i = acos(h_vect(3)/h_norm);

    type = 'N';
    if h_vect(3)>0 && i>=0 && i<=pi/2
        type = '+'; % prograde orbit
    elseif h_vect(3)==0 % do not compare i==pi/2 due to possible numerical error
        type = '0'; % polar orbit
    elseif h_vect(3)<0 && i>=pi/2 && i<=pi
        type = '-'; % retrograde orbit
    end
    if type == 'N' % check if the type of orbit has correctly been determined
        warning("There could be a problem in car2kep.m function\n");
    end

    % node line 
    n_vect = [-h_vect(2); h_vect(1); 0];
    n_vect = n_vect./norm(n_vect);
    
    % RAAN [rad]
    if i == 0
        O = 0; % by convention (RAAN not defined)
    else
        if n_vect(2) >= 0
            O = acos(n_vect(1));
        else
            O = 2*pi - acos(n_vect(1));
        end
    end
    
    % eccentricity [-]
    e_vect = (cross(v_vect, h_vect)./mu) - (r_vect./r_norm);
    e_norm = norm(e_vect);

    % argument of pericentre [rad]
    if isequal(e_vect(:),zeros(3,1))
        w = 0; % by convention
    else
        if e_vect(3) >=0
            w = acos(dot(n_vect, e_vect)/e_norm);
        else
            w = 2*pi - acos(dot(n_vect, e_vect)/e_norm);
        end
    end


    % true anomaly [rad]
    if dot(r_vect, v_vect) >=0
        theta = acos(dot(r_vect./r_norm, e_vect./e_norm));
    else
        theta = 2*pi - acos(dot(r_vect./r_norm, e_vect./e_norm));
    end


    %% store the result in ORBIT object
    ORBIT_OBJ = ORBIT(a,e_norm,i,O,w,theta,mu);
    ORBIT_OBJ.h_vect = h_vect;
    ORBIT_OBJ.r0 = r_vect;
    ORBIT_OBJ.v0 = v_vect;
    ORBIT_OBJ.e_vect = e_vect;
    ORBIT_OBJ.type = type;
end