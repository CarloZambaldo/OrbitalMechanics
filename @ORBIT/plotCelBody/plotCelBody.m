function celbody = plotCelBody(in, R_p, R_c)

%  plotCelBody.m  - Returns the 3D plot of the selected celestial body
%
% INPUT
% in        number corresponding to the desired celestial body. 
%               0) Sun;
%               1) Mercury;
%               2) Venus;
%               3) Earth;
%               4) Mars;
%               5) Jupiter;
%               6) Saturn;
%               7) Uranus;
%               8) Neptune;
%               9) Pluton.
%  R_p      radius of the celestial body, if absent or equal to 'mean', 
%           automatically setted to the mean radius in km 
%           (from Horizon, see references).
%  R_c      3 x 1 or 1 x 3 vector conteing x, y and z coordinates of the 
%           celestial body centre, if absent, automatically setted to
%           [0 0 0].
%
% OUTPUT
%  celbody  name of the plotted object, not necessary.
%
% EXAMPLE
%  earth = plotCelBody(3)
%           return the 3D plot of Earth with its mean radius and centred
%           in [0 0 0].
%  plotCelBody(4, 1e3, [1250 -100 550])
%           return the 3D plot of Mars, with a radius of 1000 and centred 
%           in [1250 -100 550].
%  Jupiter = plotCelBody(6, 'mean')
%           return the 3D plot of Jupiter, with its mean radius  and 
%           centred in [0 0 0].
%
% REFERENCES
%  https://ssd.jpl.nasa.gov/horizons/app.html#/
%
% CALLED FUNCTION
%  none
%
% AUTHOR
%  Mauro Sgura, 2021, Matlab, plotCelBody.m

switch in
    case 0
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 695700;
        end
        cdata = imread('sun.jpg');
    case 1
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 2439.7;
        end
        cdata = imread('mercury.jpg');
    case 2
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 6051.8;
        end
        cdata = imread('venus.jpg');
    case 3
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 6371.01;
        end
        cdata = imread('earth.jpg');
    case 4
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 3389.92;
        end
        cdata = imread('mars.jpg');
    case 5
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 69911;
        end
        cdata = imread('jupiter.jpg');
    case 6
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 58232;
        end
        cdata = imread('saturn.jpg');
    case 7
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 25362;
        end
        cdata = imread('uranus.jpg');
    case 8
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 24624;
        end
        cdata = imread('neptune.jpg');
    case 9
        if nargin < 2 || isequal(R_p, 'mean')
            R_p = 1151;
        end
        cdata = imread('pluton.jpg');
end

if nargin < 3
    [X, Y, Z] = ellipsoid(0, 0, 0, R_p, R_p, R_p);
else
    [X, Y, Z] = ellipsoid(R_c(1), R_c(2), R_c(3), R_p, R_p, R_p);
end

celbody = surf(X, Y, -Z);
set(celbody, 'FaceColor', 'texturemap', 'CData', cdata,...
    'EdgeColor', 'none');
axis equal

end

       