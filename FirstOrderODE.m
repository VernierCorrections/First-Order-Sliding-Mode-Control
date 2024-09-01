function [dSdt, n_c] = FirstOrderODE(t, S, N_prime, missile_saturation, T1)
    r_t = S(1:3);
    v_t = S(4:6);
    r_m = S(7:9);
    v_m = S(10:12);
    x0_m = S(13:15);
    n_t = [0; 0; 10];
    r_r = r_t - r_m;
    v_r = v_t - v_m;
    if norm(r_r) > 10^-12
        range = norm(r_r);
    else
        range = 10^-12;
    end
    closurerate = norm(v_r);
    unitlos = r_r / range;
    unitvel = v_r / closurerate;
    if norm(n_t) == 0
        unitacc = [0; 0; 0];
    else
        unitacc = n_t / norm(n_t);
    end
    t_go = -dot(v_r, r_r) / (closurerate^2);
    if t_go / T1 > (2^-6)
        x = t_go / T1;
    else
        x = 2^-6;
    end
    % n_c = N_prime * ((r_relative / t_go^2) + ((v_relative * abs(dot(unitlos, unitvel)) / t_go) + (0.5 * norm(cross(unitlos, unitacc)) * n_t)));
    N_prime = (6 * x^2 * (exp(-x) + x - 1)) / (2 * x^3 + 3 + 6 * x - 6 * x^2 - 12 * x * exp(-x) - 3 * exp(-2 * x));
    n_c = N_prime * ((r_r / t_go^2) + ((v_r * abs(dot(unitlos, unitvel)) / t_go) + (0.5 * norm(cross(unitlos, unitacc)) * n_t)) + ((-1 * x0_m * T1^2 * (exp(-x) + x - 1)) / t_go^2));
    % n_c = (N_prime / (t_go^2)) * ZEM_PLOS * unitperp;
    if norm(n_c) > missile_saturation
        n_c = (n_c / norm(n_c)) * missile_saturation;
    end
    x1_m = (-1 / T1) * x0_m + (1 / T1) * n_c; 
    dSdt = [v_t; n_t; v_m; x0_m; x1_m];
end
