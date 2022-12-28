function acc_pert_vect = acc_pert_fun_tnh(t, state, param)
    % acc_pert_vect = acc_pert_fun_tnh(t, state, param)
    %
    % INPUT:
    %   - t     : actually useless
    %   - state : [a, e, i, O, w, theta]
    %   - param : struct containing the parameters needed i.e. param.J2, R, mu
    %
    % acc_pert_fun outputs the acceleration vector
    %
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 20, 2022
    % --------------------------------------

    %% PARAMETERS
    % extracts parameters from state:
    a = state(1);
    e = state(2);
    i = state(3);
    O = state(4);
    w = state(5);
    theta = state(6);

    % extracts parameters from param:
    J2 = param.J2; % constant value
    planetR  = param.R;
    mu = param.mu;
    cD = param.cD;
    A_M = param.A_M;

    % computes other parameters
    p = a*(1-e^2);
    h = sqrt(p*mu);
    r = p/(1+e*cos(theta));
    v = sqrt(mu*(2/r-1/a));
    height = r - planetR;

    R_PF_GE = [ cos(w)*cos(O)-sin(w)*cos(i)*sin(O),	-sin(O)*cos(i)*cos(w)-cos(O)*sin(w),	sin(i)*sin(O);
            cos(O)*cos(i)*sin(w)+sin(O)*cos(w), cos(O)*cos(i)*cos(w)-sin(O)*sin(w),    -cos(O)*sin(i);
                         sin(i)*sin(w),                        sin(i)*cos(w),                cos(i)  ];
    r_vect = R_PF_GE * (r .* [cos(theta); sin(theta); zeros(1,length(theta))]);
    v_vect = R_PF_GE * (sqrt(mu/abs(p)) .* [-sin(theta); e+cos(theta); zeros(1,length(theta))]);
    h_vect = cross(r_vect, v_vect);

    % air density estimation - data from CIRA
    rho = miaCIRA(height); % [kg/m^3]

    %% COMPUTE ACCELERATION
    
    % this is in RSW
    const = -3/2 * J2 * mu * planetR^2 / r^4;
    a_r = const .* ( 1 - 3*(sin(i))^2 * (sin(theta+w))^2 );
    a_s = const .* ( (sin(i))^2 * sin(2*(theta+w)) );
    a_w = const .* ( sin(2*i) * sin(theta+w) );
    acc_pert_vect_rsw = [a_r;a_s;a_w];

    % rotate in TNH
    ROT = h/(p*v) * [e*sin(theta), -(1+e*cos(theta)), 0; 1+e*cos(theta), e*sin(theta), 0; 0,0,1];
    acc_J2 = ROT' * acc_pert_vect_rsw;
    
    % adding the aerodynamic perturbation
    w_dir = h_vect./h;
    r_dir = r_vect./norm(r_vect);
    s_dir = cross(w_dir, r_dir);
    ROT_eci_srw = [r_dir(:), s_dir(:), w_dir(:)];
    v_rel_ECI = v_vect - cross([0 0 2*pi/86164.0905]', r_vect);
    v_rel = ROT' * (ROT_eci_srw' * v_rel_ECI);
    acc_AirDrag = -.5 * A_M * cD * rho * 1e3 * norm(v_rel)^2 * v_rel./norm(v_rel);

    %% total acceleration
    acc_pert_vect = acc_J2 + acc_AirDrag;

end