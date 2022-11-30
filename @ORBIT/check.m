function [check] = check(given_ORBIT, toll)
    % checks whether the ORBIT object "obj" is defined properly
    % toll is optional, default values set to 1e-10
    %
            if nargin ~= 2
                toll = 1e-6; % tollerance set to half the machine epsilon
            end

            [test_ORBIT] = given_ORBIT.car2kep(given_ORBIT.r_vect, given_ORBIT.v_vect, given_ORBIT.mu);
            [r_vect, v_vect] = given_ORBIT.kep2car(given_ORBIT);

            test(1)  = abs(given_ORBIT.a - test_ORBIT.a) <= toll;
            test(2)  = abs(given_ORBIT.e - test_ORBIT.e) <= toll;
            test(3)  = abs(given_ORBIT.i - test_ORBIT.i) <= toll;
            test(4)  = abs(given_ORBIT.O - test_ORBIT.O) <= toll;
            test(5)  = abs(given_ORBIT.w - test_ORBIT.w) <= toll;
            test(6)  = abs(given_ORBIT.theta - test_ORBIT.theta) <= toll;
            test(7)  = abs(given_ORBIT.mu - test_ORBIT.mu) <= toll;
            if isnan(given_ORBIT.T) && isnan(test_ORBIT.T)
                test(8)=1;
            else
                test(8) = abs(given_ORBIT.T - test_ORBIT.T) <= toll;
            end
            test(9) = abs(given_ORBIT.energy - test_ORBIT.energy) <= toll;
            test(10) = abs(given_ORBIT.p - test_ORBIT.p) <= toll;
            test(11) = (norm(test_ORBIT.h_vect - given_ORBIT.h_vect) <= toll);
            test(12) = (norm(test_ORBIT.e_vect - given_ORBIT.e_vect) <= toll);
            test(13) = isequal(given_ORBIT.type, test_ORBIT.type);
           
            if isequal(test(:), ones(length(test),1))
                check = 1;
            else
                check = 0;
                fprintf("Test array: ");
                disp(test);
                warning("There could be a problem in ORBIT. Please check manually each field.");
                given_ORBIT
                test_ORBIT

                fprintf("If the warning persists try using a higher tollerance.\n");
            end
end