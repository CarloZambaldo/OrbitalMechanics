function fun = odefun_2bp_perturbed(~, y, param)
    % fun = odefun_2bp_perturbed(~, y, param)
    %
    % outputs a column vector.
    %
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 20, 2022
    % --------------------------------------
    
    % extract parameters
    J2 = param.J2; % constant value
    planetR  = param.R;
    mu = param.mu;
    cD = param.cD;
    A_M = param.A_M;

    % position and velocity
    r_vect = y(1:3);
    v_vect = y(4:6);
 
    % computes other parameters
    r_norm = norm(r_vect);
    height = r_norm - planetR;

    % air density estimation - data from CIRA
    rho = miaCIRA(height); % [kg/m^3]   

    %% COMPUTE ACCELERATION
    % J2 perturbation
    a_J2 = 3/2*J2*mu*planetR^2/r_norm^4 .* [1;1;1];
    a_J2(1) = a_J2(1) * r_vect(1)/r_norm * (5*r_vect(3)^2/r_norm^2 - 1);
    a_J2(2) = a_J2(2) * r_vect(2)/r_norm * (5*r_vect(3)^2/r_norm^2 - 1);
    a_J2(3) = a_J2(3) * r_vect(3)/r_norm * (5*r_vect(3)^2/r_norm^2 - 3);

    % air drag perturbation
    v_rel_vect = v_vect - cross([0 0 2*pi/86164.0905]', r_vect);
    v_rel = norm(v_rel_vect);
    acc_AirDrag = -.5 * A_M * cD * rho * 1e3 * v_rel^2 * v_rel_vect./v_rel; % [km/s^2]

    %% state vector
    fun = [v_vect(:); (-mu/(r_norm)^3)*r_vect(:) + a_J2 + acc_AirDrag(:)];
end