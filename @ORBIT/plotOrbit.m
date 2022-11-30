function [r_prop_vect] = plotOrbit(orbit_obj,theta_vect,d_theta,grafica,theta_segnato)
%[r_prop_vect] = plotOrbit(orbit_obj,theta_vect,d_theta,grafica,theta_segnato)
%
% plots the orbit, given the parameters ORBIT_OBJ
% using as initial and final true anomaly: theta_i and theta_f
% with a step of d_theta
%
%
% INPUT:
% orbit_obj:            ORBIT class object
% theta_vect[1x2]:      array of initial and final true anomalies
% d_theta:              true anomalt step, for plot visualization
% theta_segnato:        optional: the function will add a star in the position
%                       corresponding to:
%                       if dim(theta_segnato) = 1 --> theta = theta_segnato
%                       if dim(theta_segnato) =3 --> r = theta_segnato
% grafica:              the graphical style of the function "plot" used to plot
%                       the orbit
%
% OUTPUT:
% r_prop_vect[3xN]:     (optional, not needed to draw orbit plot) r_prop_vect = [rx; ry; rz];
%
% LOGCHANGE:
% [Guglielmo Gomiero, 26/11/2022]:
% -language fully set to english
% -added input description

%% Initialization
    a = orbit_obj.a;
    e = orbit_obj.e;
    i = orbit_obj.i;
    O = orbit_obj.O;
    w = orbit_obj.w;

    if nargin < 4 % tratto del disegno impostato di default
        grafica = '-r';
    end

    if nargin < 3 % d_theta impostato di default
        %fprintf("Default d_theta set to 2*pi/1001\n");
        d_theta = 2*pi/1001;
    end

    if nargin < 2
        theta_i = 0;
        theta_f = 2*pi;
    else
        theta_i = theta_vect(1);
        theta_f = theta_vect(2);
    end
    
    
    % Only for hyperbolas, set the thetas inside limits of theta_inf
    if e > 1
        if nargin<2 || ( theta_f >= orbit_obj.theta_inf || theta_i <= -orbit_obj.theta_inf ) 
            fprintf(2,"theta out of range, setting automatically theta inside range [-%.3f, %.3f]. \n",orbit_obj.theta_inf,orbit_obj.theta_inf);
            theta_f = orbit_obj.theta_inf*5.5/6;
            theta_i = -orbit_obj.theta_inf*5.5/6;
        end
    end


    %% Initialize theta
    tt = theta_i:d_theta:theta_f;
    % Ensura that it is a row for next calculations
    if iscolumn(tt)
        tt = tt';
    end

    %% Rotation matrix pfc --> geo
    R_PF_GE = [ cos(w)*cos(O)-sin(w)*cos(i)*sin(O),	-sin(O)*cos(i)*cos(w)-cos(O)*sin(w),	sin(i)*sin(O);
                cos(O)*cos(i)*sin(w)+sin(O)*cos(w), cos(O)*cos(i)*cos(w)-sin(O)*sin(w),    -cos(O)*sin(i);
                         sin(i)*sin(w),                        sin(i)*cos(w),                  cos(i)     ];

    %% Parameters
    p = orbit_obj.p;
    R = @(theta) (p)./(1+e.*cos(theta));


%     rx = r(tt) .* (R_PF_GE(1,1).*cos(tt) + R_PF_GE(1,2).*sin(tt));
%     ry = r(tt) .* (R_PF_GE(2,1).*cos(tt) + R_PF_GE(2,2).*sin(tt));
%     rz = r(tt) .* (R_PF_GE(3,1).*cos(tt) + R_PF_GE(3,2).*sin(tt));

    % In perifocal frame
    r_pfc = R(tt).*[cos(tt); sin(tt); zeros(1,length(tt))];
    r_geo = R_PF_GE * r_pfc;

    if nargout == 1
        r_prop_vect = r_geo;
    else
        r_prop_vect = [];
    end

    %% Draw
    hold on;
    plot3(r_geo(1,:),r_geo(2,:),r_geo(3,:), grafica, 'LineWidth', 1.3);
    axis image; 
 
    if nargin == 5 %theta_segnato is given, mark position with star
        hold on;
        % Theta_segnato is a true anomaly
        if size(theta_segnato,1) == 1 && size(theta_segnato,2) == 1
            Punto = R(theta_segnato) .* R_PF_GE * [cos(theta_segnato); sin(theta_segnato); 0];
            plot3(Punto(1), Punto(2), Punto(3), 'r*', 'LineWidth',1.3);
        % Theta_segnato is a position vector
        else
            plot3(theta_segnato(1), theta_segnato(2), theta_segnato(3), 'r*', 'LineWidth',1.3);
        end
    end


%% Textured 3D Earth example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Ryan Gray
% 8 Sep 2004
% Revised 9 March 2006, 31 Jan 2006, 16 Oct 2013

%% Options
space_color = 'k';
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
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);

globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

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

