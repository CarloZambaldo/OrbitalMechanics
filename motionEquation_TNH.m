function dstate = motionEquation_TNH(t, state, acc_pert_fun_tnh, param)
    % dstate = motionEquation(t, state, acc_pert_fun, param)
    %
    %  NOTE: THE USE OF THIS FUNCTION IS HIGHLY NOT RECOMMENDED.
    %
    % solves Gauss' Equation in TNH frame. 
    % note: this method is more computational expensive, try using
    % motionEquation_RSW instead.
    %
    % outputs the derivative of the state [a,e,i,O,w,theta]
    %
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 20, 2022
    % --------------------------------------

    warning("NOTE: THE USE OF 'motionEquation_TNH' FUNCTION IS HIGHLY NOT RECOMMENDED. Use motionEquation_RSW instead.");

    %% PARAMETERS
    % evaluating the perturbing accelerations
    acc_vect = acc_pert_fun_tnh(t, state);
    a_t = acc_vect(1);
    a_n = acc_vect(2);
    a_h = acc_vect(3);

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
    v = sqrt(mu*(2/r-1/a));
    h = sqrt(p*mu);

    %% compute the variation of the states in TNH frame
    da     = 2*a^2*v/mu * a_t;
    de     = 1/v * (2 * (e+cos(theta) * a_t) - r/a * sin(theta) * a_n);
    di     = r * cos(theta+w)/h * a_h;
    dO     = r * sin(theta+w)/(h*sin(i)) * a_h;
    dw     = 2/(e*v) * (2*sin(theta) * a_t + (2*e + r/a * cos(theta)) * a_n) - r * sin(theta+w)*cos(i)/(h*sin(i)) * a_h;
    dtheta = h/r^2 - 1/(e*v) * (2*sin(theta)*a_t + (2*e + r/a * cos(theta)) * a_n);
   
    % build the dstate
    dstate = [da; de; di; dO; dw; dtheta];
end
