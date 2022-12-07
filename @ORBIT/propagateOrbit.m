function [t,y] = propagateOrbit(orbit_obj,y_0,time_vect,N_dt,grafica,r_segnato)

%    [t,y] = propagateOrbit(orbit_obj,y_0,time_vect,N_dt,grafica,r_segnato)
%
% plots the orbit, given the parameters ORBIT_OBJ
% using as initial and final true anomaly: theta_i and theta_f
% with a step of d_theta
%
% INPUT:
% orbit_obj:            ORBIT class object
% y0[6x1]:              initial state, equatorial cartesian coordinates
% time_vect[1x2]:       tspan
% N_dt:                 discretization step for plotting
% grafica:              the graphical style of the function "plot" used to plot
%                       the orbit
% r_segnato[3x1]:       optional: the function will add a star in the position
%                       corresponding given by r_segnato
%
% OUTPUT:
% t[Nx1]:
% y[Nx6]:     output of ode integration
%
% LOGCHANGE:
% [Guglielmo Gomiero, 26/11/2022]:
%   - language fully set to english
%   - added input description
% [Guglielmo Gomiero, 27/11/2022]:
%   - commented out plotting
%   - commented out third party Earth plotting function, as it causes problems
%     with the variable y(declared 2 times, y is also the output)
% [Carlo Zambaldo, 06/12/2022]:
%   - added plotting only if nargout <= 1, fixed "y" bug

    y_0 = y_0(:);

    if nargin < 3
        t_1 = 0;
        t_2 = orbit_obj.T;
    else
        t_1 = time_vect(1);
        t_2 = time_vect(2);
    end
    if nargin < 5 % graphics set by default
        grafica = '-r';
    end
    if nargin < 4 % number of intervals set by default
        N_dt = 1e4;
    end


    %% propagate orbit
    tspan = linspace(t_1, t_2, N_dt);
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    [t,y] = ode113( @(t,y) odefun_2bp(t, y, orbit_obj.mu), tspan, y_0, options);
   
    %% disegno
    if nargout <= 1
        hold on;
        plot3(y(:,1), y(:,2), y(:,3), grafica, 'LineWidth',1.3);
        axis image; 
     
        if nargin == 5 % significa che si è inserito anche r_segnato (metto un asterisco in r_segnato)
            hold on;
            plot3(r_segnato(1), r_segnato(2), r_segnato(3), 'r*', 'LineWidth',1.3);
        end
    
    
        %% Textured 3D Earth example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Ryan Gray
        % 8 Sep 2004
        % Revised 9 March 2006, 31 Jan 2006, 16 Oct 2013
        
        %% Options
        npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
        alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
        GMST0 = []; % Don't set up rotatable globe (ECEF)
        %GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0
        
        % Earth texture image
        % Anything imread() will handle, but needs to be a 2:1 unprojected globe
        % image.
        
        image_file = '1024px-Land_ocean_ice_2048.jpg';
        
        % Mean spherical earth
        erad = 6371008.7714/1000; % equatorial radius (km)
        prad = 6371008.7714/1000; % polar radius (km)
        erot = 7.2921158553e-5; % earth rotation rate (radians/sec)
        
        %% Create figure
        hold on;
        %set(gca, 'NextPlot','add', 'Visible','off');
        grid on;
        
        % Set initial view
        view(30,30);
        
        
        %% Create wireframe globe
        % Create a 3D meshgrid of the sphere points using the ellipsoid function
        [xea, yea, zea] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
        
        globe = surf(xea, yea, -zea, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
        
        if ~isempty(GMST0)
            hgx = hgtransform;
            set(hgx,'Matrix', makehgtform('zrotate',GMST0));
            set(globe,'Parent',hgx);
        end
        
        %% Texturemap the globe
        % Load Earth image for texture map
        cdata = imread(image_file);
        
        % Set image as color data (cdata) property, and set face color to indicate
        % a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
        set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        
        
        %% legend on axis
        xlabel("X [km]");
        ylabel("Y [km]");
        zlabel("Z [km]");
    end
end
