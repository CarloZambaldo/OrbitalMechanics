function fun = odefun_2bp(~, y, mu)
% fun = odefun_2bp(~, y, mu)
%
% outputs a column vector.
%

    % position and velocity
    r = y(1:3);
    v = y(4:6);

    % state vector
    fun = [v; (-mu/(norm(r))^3)*r];
end