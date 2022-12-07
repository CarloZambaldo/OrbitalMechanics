classdef ASTRO
    % ASTRO class defines a variable of a celestial body selected by its
    % name simply calling variable = ASTRO('planetName');
    %  the variable contains:
    %    - planetName   
    %    - mu           body constant (mu = mass * G) [km^3/s^2]:
    %    - R            radius of the body
    %    - r_SOI        radius of sphere of influence

    properties
        name
        mu      {mustBeNumeric}
        R       {mustBeNumeric}
        r_SOI   {mustBeNumeric}
    end

    methods
        function [obj] = ASTRO(name)
            if nargin ~= 0
                switch upper(name)
                    case 'EARTH'
                        obj.mu = astroConstants(13);
                        obj.R  = astroConstants(23);
                    case 'SUN'
                        obj.mu = astroConstants(4);
                        obj.R = astroConstants(3);
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
                        error("Error occurred in the definition of type ASTRO variable");
                end
                obj.name = upper(name);
                obj.r_SOI = obj.R * (obj.mu/astroConstants(4)) ^ (2/5);
            end
        end
    end
end
