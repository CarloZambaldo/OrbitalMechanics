function [a_req, T_req] = repeating_groundTrack(k,m,omega_planet,mu)
    % [a_req, T_req] = repeating_groundTrack(k,m,omega_planet,mu)
    % 
    % INPUTs:
    %   - k             = number of satellite revolutions
    %   - m             = number of Earth revolutions
    %   - omega_planet  = angular velocity of the planet
    %   - mu            = astro body constant (mu = mass * G)
    % 
    % OUTPUTs:
    %   - a_req         = required semimajor axis
    %   - T_req         = required orbital periot
    %  to have a repeating groundtrack
    %
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 20, 2022
    % --------------------------------------
        
    T_req = m/k * 2*pi/omega_planet;
    n = omega_planet*k/m;
    a_req = (mu/n^2)^(1/3);
end