function [rho, rho0, height0, scale_height] = miaCIRA(altitude)
    % [rho, rho0, height0, scale_height] = miaCIRA(altitude)
    % 
    % this function uses the given table (see on slides) to compute the air
    % density (CIRA model).
    % The exponential model is used.
    %
    % WARNING: the output density is in kg/m^3 !!!
    %
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 20, 2022
    % --------------------------------------


    if altitude <0
        error("altitude < 0. Required altitude: "+num2str(altitude));
    end
    
    legend_height = [0 25 30:10:150 180 200:50:500 600:100:1000 1000]';
    
    rho0 = [1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, 3.206e-4, 8.77e-5, 1.905e-5, 3.396e-6, 5.297e-7, 9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, ...
            2.070e-9, 5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, 9.158e-12, 3.725e-12, 1.585e-12, 6.967e-13, 1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15, 3.019e-15]';
    
    scale_height = [7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, 5.799, 5.382, 5.877, 7.263, 9.473, 12.636, 16.149, 22.523, 29.740, 37.105, 45.546, 53.628, 53.298, ...
                    58.515, 60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268.00]';
    
    index = 1;
    while index < length(legend_height)-1 && altitude > legend_height(index+1)
        index = index+1;
    end
    
    rho0  = rho0(index);
    height0 = legend_height(index);
    scale_height = scale_height(index);

    rho = rho0 * exp(-(altitude-height0)/scale_height);

end
