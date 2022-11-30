function delta_time = TOF(obj, theta_vect)
    % delta_time = TOF(orbit, theta_vect)
    % Calculate time of flight of orbit given by ORBIT object, starting
    % from theta_vect(1) and ending on theta_vect(2)
    %
    % INPUT:
    % theta_vect = [theta_i, theta_f]
    %
    % OUTPUT:
    % delta_time
    %
    % CHANGELOG:
    % 


    %% initialise
    
    theta_i = theta_vect(1);
    theta_f = theta_vect(2);
    

    if theta_i > theta_f
        error("theta_i > theta_f");
    end
    
    if obj.e < 1
        while theta_i < 0
            % note that this function only works with positive angles.
            theta_i = wrapTo2Pi(theta_i);
            theta_f = wrapTo2Pi(theta_f);
        end
    end

    warning("Inserire logica se theta_i > theta_f");
    n_comp_rev = floor((theta_f-theta_i)/(2*pi));
    theta_f = theta_f - 2*pi*n_comp_rev;

    if obj.e>=1 && n_comp_rev~=0
        error("The given angles would imply multiple passages around the body, which is not consistent with the given trajectory");
    end
    
    %% orbit type 

    % CIRCUMFERENCE
    if obj.e==0
        t_1 = 0;
        t_2 = ((obj.mu*obj.a)^(3/2))/(obj.mu^2);
    
    % ELLIPSE 
    elseif obj.e<1
        E   = wrapTo2Pi(2*atan( sqrt((1-obj.e)/(1+obj.e)) * tan(theta_i/2) ))
        t_1 = obj.T/(2*pi) * ( E - obj.e*sin(E) )

        E   = wrapTo2Pi(2*atan( sqrt((1-obj.e)/(1+obj.e)) * tan(theta_f/2) ))
        t_2 = obj.T/(2*pi) * ( E - obj.e*sin(E) )
   
    % PARABOLA
    elseif obj.e==1
        D = tan(theta_i/2);
        t_1 = .5 * sqrt(obj.p^3/obj.mu) * (D + D^3/3);
    
        D = tan(theta_f/2);
        t_2 = .5 * sqrt(obj.p^3/obj.mu) * (D + D^3/3);

    % HYPERBOLA
    elseif obj.e>1
        F   = 2*atanh( sqrt(-(1-obj.e)/(1+obj.e)) * tan(theta_i/2) );
        t_1 = sqrt(abs(obj.a^3/obj.mu)) * (obj.e*sin(F) - F);

        F   = 2*atanh( sqrt(-(1-obj.e)/(1+obj.e)) * tan(theta_f/2) );
        t_2 = sqrt(abs(obj.a^3/obj.mu)) * (obj.e*sin(F) - F);
    end


    %% Compute delta_time
    if t_1 > t_2 || t_1 < 0 || t_2<0
        warning("There could be an issue in TOF function.");
        t_2 = t_2 + obj.T;
    end

    delta_time = t_2 - t_1;

    if theta_i > theta_f
        delta_time = delta_time + obj.T;
    end
    if obj.e<1
        delta_time = n_comp_rev*obj.T + delta_time;
    end
end