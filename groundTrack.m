function [alpha, delta, lon, lat] = groundTrack(ORBIT_OBJ,theta_G_0,omega_planet,N,plot_options)
    %
    %
    %
    %
    % example of plot_orbit optional field:
    %    plot_options.LineWidth = 1.1;
    %    plot_options.Color = 'g';

    %% setup
    theta_G = @(t) theta_G_0 + omega_planet * t;
    
    mu = ORBIT_OBJ.mu;
    y_0 = [ORBIT_OBJ.r0; ORBIT_OBJ.v0]; % state vector
    
    %% ode solver setup and integration
    tspan = linspace(0, N*ORBIT_OBJ.T, 1e5);
    options = odeset('RelTol', 1e-7);
    
    [t_vect,y] = ode113( @(t,y) odefun_2bp(t, y, mu), tspan, y_0, options);



    %%
    r_vect_prop = y(:,1:3)'; % extract the position vectors for all the times
    rx = r_vect_prop(1,:);
    ry = r_vect_prop(2,:);
    rz = r_vect_prop(3,:);
    r_norm = sqrt(sum(r_vect_prop.^2,1)); % take the norm

    % angles
    delta = asin(rz./r_norm);
    alpha = atan2(ry,rx);

    % from here angles in degrees
    lon = wrapTo180((alpha(:) - theta_G(t_vect(:))) .* 180/pi);
    lat = wrapTo180(delta(:) .* 180/pi);

    if nargin < 5
        plot_options.LineWidth = 1.1;
        plot_options.Color = 'g';
        plot_options.planet = 'Earth';
    end
    ground_track(lat,lon,plot_options);
    title("a = "+num2str(ORBIT_OBJ.a)+"[km], e = "+num2str(ORBIT_OBJ.e)+"[-], i = "+num2str(ORBIT_OBJ.i*180/pi)+"[deg], \Omega = "+num2str(ORBIT_OBJ.O*180/pi)+"[deg], \omega = "+num2str(ORBIT_OBJ.w*180/pi)+"[deg]");






