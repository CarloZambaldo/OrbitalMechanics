function [delta_theta, theta_vect] = AOF(obj, time_vect)
    % [delta_theta, theta_vect] = AOF(obj, time_vect)
    %
    % the function is not finished:
    % to do:
    %   - compute initial and final thetas,
    %   - check various possible combinations of errors (i.e. theta_f <
    %     theta_i....)
    %   - allow different types of orbit.
    
    error("This function has not been implemented yet");

    time_i = time_vect(1);
    time_f = time_vect(2);
    delta_time = time_f - time_i;

    f = @(E) E - obj.e*sin(E) - sqrt(obj.mu/obj.a^3) * delta_time;
   
    xvect = fzero(f,0);

    theta_i = 0;
    theta_f = 2*atan(sqrt((1+obj.e)/(1-obj.e))*tan(xvect(end)/2));

    delta_theta = theta_f - theta_i;

    theta_vect = [theta_i; theta_f];
    

end