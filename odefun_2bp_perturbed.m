function fun = odefun_2bp_perturbed(~, y, mu)
% fun = odefun_2bp(~, y, mu)
%
% outputs a column vector.
%

    J2 = 0.00108263; % constant value
    R_E = astroConstants(23);

    % position and velocity
    r = y(1:3);
    v = y(4:6);
    
    r_norm = norm(r);
    a_J2 = 3/2*J2*mu*R_E^2/r_norm^4 .* [1;1;1];
    a_J2(1) = a_J2(1) * r(1)/r_norm * (5*r(3)^2/r_norm^2 - 1);
    a_J2(2) = a_J2(2) * r(2)/r_norm * (5*r(3)^2/r_norm^2 - 1);
    a_J2(3) = a_J2(3) * r(3)/r_norm * (5*r(3)^2/r_norm^2 - 3);

    % state vector
    fun = [v; (-mu/(r_norm)^3)*r + a_J2];
