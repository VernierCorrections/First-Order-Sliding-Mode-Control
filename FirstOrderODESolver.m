function [target_rmat, missile_rmat, time_mat, acceleration_mat, dv] = FirstOrderODESolver(t0, tf, S0, N_prime, missile_saturation, T1)
beta = 0.8;
Fmin = 0.125;
Fmax = 4.0;
attemptspertimestep = 12;
epsilon_relative = 10^(-9);
epsilon_absolute = 10^(-9);
min_error = 10^(-18);
t_guess = tf;
h_min = 10^(-16);
stepsizecontrol = [beta; Fmin; Fmax; attemptspertimestep; epsilon_relative; epsilon_absolute; min_error; t_guess; h_min];
h_initial = 10^(-9);
t = t0;
S = S0;
i = 1;
h = h_initial;
dv = 0;
target_rmat = zeros(3, 10^4);
missile_rmat = zeros(3, 10^4);
time_mat = zeros(1, 10^4);
acceleration_mat = zeros(1, 10^4);
target_rmat(:, 1) = S(1:3);
missile_rmat(:, i) = S(7:9);
while true
    [t, S, i, h, a] = FirstOrderRKF78(t, S, N_prime, missile_saturation, T1, i, h, stepsizecontrol);
    target_rmat(:, i) = S(1:3);
    missile_rmat(:, i) = S(7:9);
    time_mat(1, i) = t;
    acceleration_mat(1, i) = a;
    dv = dv + a * h;
    if t >= tf
        break
    end
    %if any(isnan(S)) == true
    %    break
    %end
    r_t = S(1:3);
    v_t = S(4:6);
    r_m = S(7:9);
    v_m = S(10:12);
    r_r = (r_t - r_m);
    v_r = (v_t - v_m);
    % if dot(r_r, v_r) > 0
    %     break
    % end
end
end


