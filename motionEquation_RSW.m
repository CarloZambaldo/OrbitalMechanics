function dstate = motionEquation_RSW(t, state, acc_pert_fun_rsw, param)
    % dstate = motionEquation(t, state, acc_pert_fun, param)
    %
    % solves Gauss' Equation in TNH frame
    % outputs the derivative of the state [a,e,i,O,w,theta]
    %
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 20, 2022
    % --------------------------------------
    
    %% PARAMETERS
    % evaluating the perturbing accelerations
    acc_vect = acc_pert_fun_rsw(t, state);
    
    % extracts parameters from state:
    a = state(1);
    e = state(2);
    i = state(3);
    O = state(4);
    w = state(5);
    theta = state(6);
    mu = param.mu;

    % compute other parameters
    p = a*(1-e^2);
    r = p/(1+e*cos(theta));
    h = sqrt(p*mu);

    % compute the variation of the states
    da = 2*a^2/h * (e * sin(theta) * acc_vect(1) + p/r * acc_vect(2));
    de = 1/h * (p*sin(theta)*acc_vect(1) + ((p+r)*cos(theta)+r*e) * acc_vect(2));
    di = r * cos(theta+w)*acc_vect(3)/h;
    dO = r * sin(theta+w)/(h*sin(i)) * acc_vect(3);
    dw = 1/(h*e) * (-p * cos(theta)*acc_vect(1) + (p+r)*sin(theta)*acc_vect(2)) - r*sin(theta+w)*cos(i)/(h*sin(i))*acc_vect(3);
    dtheta = h/r^2 + 1/(e*h) * (p*cos(theta)*acc_vect(1) - (p+r)*sin(theta)*acc_vect(2));
    
    % build the dstate
    dstate = [da; de; di; dO; dw; dtheta];
end
