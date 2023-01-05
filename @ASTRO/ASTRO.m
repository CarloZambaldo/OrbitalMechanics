classdef ASTRO
    % ASTRO class defines a variable of a celestial body selected by its
    % name simply calling variable = ASTRO('planetName');
    %  the variable contains:
    %    - planetName   
    %    - mu           body constant (mu = mass * G) [km^3/s^2]:
    %    - R            radius of the body
    %    - pos          position of the Astro (useful for plotAstro)
    %
    %
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 20, 2022
    % --------------------------------------
    %
    %  TO DO:
    %   add a method "ephemeris" which uses uplanet or ephNEO to compute the
    %   ephemeris of the astro object

    properties
        name                        % name of the astro
        mu      {mustBeNumeric}     % body constant (mu = mass * G)
        R       {mustBeNumeric}     % radius of the body
        pos     (3,1){mustBeVector} % position of the planet (default [0;0;0])
        omega   {mustBeNumeric}     % rotation speed of planet on its axis [rad/s]
        atmosphere                  % [km] height of the atmosphere of the planet
    end

    methods
        function [obj] = ASTRO(name,position)
            if nargin ~= 0
                switch upper(name)
                    case 'EARTH'
                        obj.mu = astroConstants(13);
                        obj.R  = astroConstants(23);
                        obj.omega = 2*pi/(3600*23.9344696); % (sidereal day)
                        obj.atmosphere = 200;
                    case 'SUN'
                        obj.mu = astroConstants(4);
                        obj.R = astroConstants(3);
                        obj.atmosphere = 0;
                    case 'MERCURY'
                        obj.mu = astroConstants(11);
                        obj.R = astroConstants(21);
                    case 'VENUS'
                        obj.mu = astroConstants(12);
                        obj.R = astroConstants(22);
                    case 'MARS'
                        obj.mu = astroConstants(14);
                        obj.R = astroConstants(23);
                    case 'JUPITER'
                        obj.mu = astroConstants(15);
                        obj.R = astroConstants(25);
                    case 'SATURN'
                        obj.mu = astroConstants(16);
                        obj.R = astroConstants(26);
                    case 'URANUS'
                        obj.mu = astroConstants(17);
                        obj.R = astroConstants(27);
                    case 'NEPTUNE'
                        obj.mu = astroConstants(18);
                        obj.R = astroConstants(28);
                    case 'PLUTO'
                        obj.mu = astroConstants(19);
                        obj.R = astroConstants(29);
                    case 'MOON'
                        obj.mu = astroConstants(20);
                        obj.R = astroConstants(30);
                    otherwise
                        warning("Astro named """+name+""" not found using Earth as default.");
                        obj.mu = astroConstants(13);
                        obj.R  = astroConstants(23);
                        name = "noname";
                end
                obj.name = upper(name);

                if nargin > 1
                    obj.pos = position;
                end
            end
        end

        function [astro_plot] = plotAstro(obj, scale, position)
            % Options
            npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
            alpha   = 1;     % globe transparency level, 1 = opaque, through 0 = invisible
            body_pos = obj.pos;

            if nargin>1
                rho = scale * obj.R;
                if nargin>2
                    body_pos = position(:);
                end
            else
                rho = obj.R; % if no scale use default radius of planet
            end
            
            if obj.name == ""
                error("ASTRO object name not defined.");
            end

            % reads in image
            cdata = imread(strcat(lower(obj.name),'.png'));
            
            % set space
            hold on;
            grid on;
            axis image;
            view(30,30);
            
            % draw planet
            [xea, yea, zea] = ellipsoid(body_pos(1), body_pos(2), -body_pos(3), rho, rho, rho, npanels);
            astro_plot = surf(xea, yea, -zea, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
            set(astro_plot, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    
            % legend on axis
            xlabel("X [km]");
            ylabel("Y [km]");
            zlabel("Z [km]");
        end

        function [obj] = updatePosition(obj, pos)
            % right now this function is more useless than useful but who
            % knows... maybe it will be useful.
            obj.pos = pos(:);
        end
    end

end
