function [delta_time, intermTimes] = TOF(obj, theta_vect)
    % [delta_time, intermTimes] = TOF(obj, theta_vect)
    %
    % Calculate time of flight of orbit given by ORBIT object, starting
    % from theta_vect(1) and ending on theta_vect(2)
    %
    % INPUT:
    %   theta_vect = [theta_i, theta_f]  -> VECTOR IN RADIANS
    %
    % OUTPUT:
    %   delta_time = ToF between theta_i and theta_f on orbit "obj"
    %
    % CHANGELOG:
    %  [Carlo Zambaldo 09/12/2022]:
    %    - various comments
    %    - fixed ToF for hyperbolic trajectories (use sinh not sin!)
    %    - added intermediate times as output if nargout >1
    %    - added some checks and errors


    %% initialise
    theta_i = theta_vect(1);
    theta_f = theta_vect(2);
    
    if obj.e < 1 % for ellipses and circles
        %if theta_i < 0 || theta_f < 0 
        %    % note that this function only works with positive angles.
        %    theta_i = (theta_i);
        %    theta_f = (theta_f);
        %end
    elseif obj.e>1 && (abs(theta_i)>obj.theta_inf || abs(theta_f)>obj.theta_inf)
        error("The given thetas are out of range for this hyperbolic trajectory.");
    end

    % number of complete revolutions
    n_comp_rev = floor((theta_f-theta_i)/(2*pi));
    theta_f = theta_f - 2*pi*n_comp_rev; % remove the revolutions from theta_f

    if obj.e>=1 && n_comp_rev~=0 % multiple passages for hyperbolas and parabolas not allowed
        error("The given angles would imply multiple passages around the body, which is not consistent with the given trajectory");
    end
    
    %% orbit type 

    % CIRCUMFERENCE
    if obj.e==0
        t_1 = 0;
        t_2 = ((obj.mu*obj.a)^(3/2))/(obj.mu^2);
    
    % ELLIPSE 
    elseif obj.e<1
        E_1_unwrapped = 2*atan( sqrt((1-obj.e)/(1+obj.e)) * tan(theta_i/2) );
        E_1   = wrapTo2Pi(E_1_unwrapped);
        t_1 = obj.T/(2*pi) * ( E_1 - obj.e*sin(E_1) );
        if E_1_unwrapped ~= E_1
            t_1 = t_1 - obj.T;
        end

        E_2_unwrapped = wrapTo2Pi(2*atan( sqrt((1-obj.e)/(1+obj.e)) * tan(theta_f/2) ));
        E_2   = wrapTo2Pi(E_2_unwrapped);
        t_2 = obj.T/(2*pi) * ( E_2 - obj.e*sin(E_2) );
        if E_2_unwrapped ~= E_2
            t_2 = t_2 - obj.T;
        end
   
    % PARABOLA
    elseif obj.e==1
        D_1 = tan(theta_i/2);
        t_1 = .5 * sqrt(obj.p^3/obj.mu) * (D_1 + D_1^3/3);
    
        D_2 = tan(theta_f/2);
        t_2 = .5 * sqrt(obj.p^3/obj.mu) * (D_2 + D_2^3/3);

    % HYPERBOLA
    elseif obj.e>1
        F_1 = 2*atanh( sqrt(-(1-obj.e)/(1+obj.e)) * tan(theta_i/2) );
        t_1 = sqrt(abs(obj.a^3/obj.mu)) * (obj.e*sinh(F_1) - F_1);

        F_2 = 2*atanh( sqrt(-(1-obj.e)/(1+obj.e)) * tan(theta_f/2) );
        t_2 = sqrt(abs(obj.a^3/obj.mu)) * (obj.e*sinh(F_2) - F_2);
    end


    %% Compute delta_time
    if t_1 > t_2 || t_1 < 0 || t_2<0
        warning("There could be an issue in TOF function.");
        %t_2 = t_2 + obj.T;
    end

    delta_time = t_2 - t_1;

    if theta_i > theta_f
        warning("theta_i > theta_f");
        delta_time = delta_time + obj.T;
    end
    
    if obj.e<1
        delta_time = n_comp_rev*obj.T + delta_time;
    end

    if nargout>1
        intermTimes = [t_1; t_2];
    end
end