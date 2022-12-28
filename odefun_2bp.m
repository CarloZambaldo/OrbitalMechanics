function fun = odefun_2bp(~, y, param)
    % fun = odefun_2bp(~, y, param)
    %
    % outputs a column vector.
    % param is a struct containing the mu parameter
    %
    % --------------------------------------
    % Carlo Zambaldo (info@carlozambaldo.it)
    %       last update: Dec 20, 2022
    % --------------------------------------

    % position and velocity
    r = y(1:3);
    v = y(4:6);

    % extract parameters
    if isempty(param) || ~isfield(param,'mu')
         % note: this is to maintain backward compatibility with code 
         % programmed with an older version of this function, which once
         % required mu as input and not a struct containing mu (param.mu).
         mu = param;
    else
        mu = param.mu;
    end

    % state vector
    fun = [v; (-mu/(norm(r))^3)*r];
end