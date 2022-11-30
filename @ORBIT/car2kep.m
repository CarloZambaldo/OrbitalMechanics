function [a,e_norm,i,O,w,theta] = car2kep(r_vect, v_vect, mu)
    % [orbit_obj] = car2kep(r_vect, v_vect, mu)
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
    %     - a:              semi-major axis
    %       e_norm:         eccentricity
    %       i:              inclination
    %       O:              RAAN
    %       w:              argument of pericenter
    %       theta:          true anomaly
    %
    %
    % CHANGELOG:
    % [Guglielmo Gomiero 26/11/2022]:
    %   - modified to output Keplerian elements as an array
    % [Guglielmo Gomiero 27/11/2022]:
    %   - fixed w calculation when i = 0
    % [Carlo Zambaldo 30/11/2022]
    %   - added "=" to expression "elseif e_vect(3) > 0"

    % tollerance for numerical error
    toll = 1e-13;

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

    % node line (note: unitary vector)
    n_vect = [-h_vect(2); h_vect(1); 0]; % from cross prod: cross(k,h)
    n_vect = n_vect./norm(n_vect);
    
    % RAAN [rad]
    if i <= toll
        i = 0; % to "remove" numerical error
        O = 0; % by convention (RAAN not defined)
        n_vect = [1;0;0];
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
    if e_norm <= toll % toll due to possible numerical error
        e_norm = 0;
        e_vect = zeros(3,1);
        w = 0; % by convention
        e_direction = n_vect; % by convention
    else
        e_direction = e_vect./e_norm;
        if e_vect(3) == 0
            if e_vect(2) >= 0
                w = acos(dot(n_vect, e_direction));
            else
                w = 2*pi - acos(dot(n_vect, e_direction));
            end
        elseif e_vect(3) >= 0
            w = acos(dot(n_vect, e_direction));
        else
            w = 2*pi - acos(dot(n_vect, e_direction));
        end
    end

    % true anomaly [rad]
    if dot(r_vect, v_vect) >=0
        theta = acos(dot(r_vect./r_norm, e_direction));
    else
        theta = 2*pi - acos(dot(r_vect./r_norm, e_direction));
    end


    %% store the result in ORBIT object
    if e_norm==1
        a = .5*h_norm^2/mu;  % for parabolas use a to define rP
    end
end