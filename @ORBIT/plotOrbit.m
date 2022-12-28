function [r_prop_vect] = plotOrbit(obj,theta_vect,d_theta,grafica,theta_segnato)
%[r_prop_vect] = plotOrbit(obj,theta_vect,d_theta,grafica,theta_segnato)
%
% plots the orbit, given the parameters ORBIT_OBJ
% using as initial and final true anomaly: theta_i and theta_f
% with a step of d_theta
%
%
% INPUT:
%   - obj                :  ORBIT class object
%   - theta_vect         :  [1x2] array of initial and final true anomalies
%   - d_theta            :  true anomalt step, for plot visualization
%   - theta_segnato      :  optional - the function will add a star in the position
%                           corresponding to:
%                           if theta_segnato = [1x1] --> theta = theta_segnato
%                           if theta_segnato = [1x3] --> r = theta_segnato
%   - grafica:              the graphical style of the function "plot" used to plot
%                           the orbit
%
% OUTPUT:
%   - r_prop_vect[3xN]   :  (optional) r_prop_vect = [rx; ry; rz];
%
% LOGCHANGE:
%  [Guglielmo Gomiero, 26/11/2022]:
%   - language fully set to english
%   - added input description
%  [Carlo Zambaldo, 07/12/2022]:
%   - removed planet plotting: use instead ASTRO to define an ASTRO object
%     and then call ASTRO.plotAstro method to draw the celestial body


%% Initialization
    a = obj.a;
    e = obj.e;
    i = obj.i;
    O = obj.O;
    w = obj.w;

    if nargin < 4 % tratto del disegno impostato di default
        grafica = '-';
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
        if nargin<2 || ( theta_f >= obj.theta_inf || theta_i <= -obj.theta_inf ) 
            fprintf(2,"theta out of range, setting automatically theta inside range [-%.3f, %.3f]. \n",obj.theta_inf,obj.theta_inf);
            theta_f = obj.theta_inf*5.5/6;
            theta_i = -obj.theta_inf*5.5/6;
        end
    end


    %% Initialize theta
    tt = theta_i:d_theta:theta_f;
    % Ensura that it is a row for next calculations
    tt = tt(:)';

    %% Rotation matrix pfc --> geo
   R_PF_GE = [ cos(w)*cos(O)-sin(w)*cos(i)*sin(O),	-sin(O)*cos(i)*cos(w)-cos(O)*sin(w),	sin(i)*sin(O);
               cos(O)*cos(i)*sin(w)+sin(O)*cos(w), cos(O)*cos(i)*cos(w)-sin(O)*sin(w),    -cos(O)*sin(i);
                        sin(i)*sin(w),                        sin(i)*cos(w),                  cos(i)     ];
 
    %% Parameters
    p = obj.p;
    R = @(theta) (p)./(1+e.*cos(theta));

    % In perifocal frame
    r_pfc = R(tt).*[cos(tt); sin(tt); zeros(1,length(tt))];
    r_geo = R_PF_GE * r_pfc;
    
    %     r_geo(1,:) = R(tt) .* (R_PF_GE(1,1).*cos(tt) + R_PF_GE(1,2).*sin(tt));
    %     r_geo(2,:) = R(tt) .* (R_PF_GE(2,1).*cos(tt) + R_PF_GE(2,2).*sin(tt));
    %     r_geo(3,:) = R(tt) .* (R_PF_GE(3,1).*cos(tt) + R_PF_GE(3,2).*sin(tt));

    if nargout == 1
        r_prop_vect = r_geo;
    else
        r_prop_vect = [];
    end

    %% Draw
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
        hold;
    end
end
