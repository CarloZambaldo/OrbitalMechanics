function [delta_th] = delta_theta(obj, delta_time,toll)


    GM = 398600.4418; % graviational parameter in km^3/s^2
    
    % metodo di Newton
    f = @(E) E - e*sin(E) - sqrt(GM/a^3) * delta_t;
    df = @(E) 1 - e*cos(E);
    
    nmax = ceil(sqrt(1/toll));

    [xvect, ~] = metodo_newton(x0,nmax,toll,f,df,1);
    
    delta_th = 2*atan(sqrt((1+e)/(1-e))*tan(xvect(end)/2));