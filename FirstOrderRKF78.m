function [t, S, i, h, a] = FirstOrderRKF78(t, S, N_prime, missile_saturation, T1, i, h, stepsizecontrol)
    beta = stepsizecontrol(1);
    Fmin = stepsizecontrol(2);
    Fmax = stepsizecontrol(3);
    attemptspertimestep = stepsizecontrol(4);
    epsilon_relative = stepsizecontrol(5);
    epsilon_absolute = stepsizecontrol(6);
    min_error = stepsizecontrol(7);
    t_guess = stepsizecontrol(8);
    h_min = stepsizecontrol(9);
    h_temp = h;
    j = 1;
    while true
        t0         = t;
        S0         = S;
        [k0, a0]   = FirstOrderODE(t0, S0, N_prime, missile_saturation, T1);
        t1         = t + 2/27 * h_temp;
        S1         = S + h_temp * (2/27 * k0);
        [k1, ~]    = FirstOrderODE(t1, S1, N_prime, missile_saturation, T1);
        t2         = t + 1/9 * h_temp;
        S2         = S + h_temp * (1/36 * k0 + 1/12 * k1);
        [k2, ~]    = FirstOrderODE(t2, S2, N_prime, missile_saturation, T1);
        t3         = t + 1/6 * h_temp;
        S3         = S + h_temp * (1/24 * k0 + 0.0 * k1 + 1/8 * k2);
        [k3, ~]    = FirstOrderODE(t3, S3, N_prime, missile_saturation, T1);
        t4         = t + 5/12 * h_temp;
        S4         = S + h_temp * (5/12 * k0 + 0.0 * k1 + (-25/16) * k2 + 25/16 * k3);
        [k4, ~]    = FirstOrderODE(t4, S4, N_prime, missile_saturation, T1);
        t5         = t + 1/2 * h_temp;
        S5         = S + h_temp * (1/20 * k0 + 0.0 * k1 + 0.0 * k2 + 1/4 * k3 + 1/5 * k4);
        [k5, a5]   = FirstOrderODE(t5, S5, N_prime, missile_saturation, T1);
        t6         = t + 5/6 * h_temp;
        S6         = S + h_temp * (-25/108 * k0 + 0.0 * k1 + 0.0 * k2 + 125/108 * k3 + (-65/27) * k4 + 125/54 * k5);
        [k6, a6]   = FirstOrderODE(t6, S6, N_prime, missile_saturation, T1);
        t7         = t + 1/6 * h_temp;
        S7         = S + h_temp * (31/300 * k0 + 0.0 * k1 + 0.0 * k2 + 0.0 * k3 + 61/225 * k4 + (-2/9) * k5 + 13/900 * k6);
        [k7, a7]   = FirstOrderODE(t7, S7, N_prime, missile_saturation, T1);
        t8         = t + 2/3 * h_temp;
        S8         = S + h_temp * (2 * k0 + 0.0 * k1 + 0.0 * k2 + (-53/6) * k3 + 704/45 * k4 + (-107/9) * k5 + 67/90 * k6 + 3 * k7);
        [k8, a8]   = FirstOrderODE(t8, S8, N_prime, missile_saturation, T1);
        t9         = t + 1/3 * h_temp;
        S9         = S + h_temp * ((-91/108) * k0 + 0.0 * k1 + 0.0 * k2 + 23/108 * k3 + (-976/135) * k4 + 311/54 * k5 + (-19/60) * k6 + 17/6 * k7 + (-1/12) * k8);
        [k9, a9]   = FirstOrderODE(t9, S9, N_prime, missile_saturation, T1);
        t10        = t + 1 * h_temp;
        S10        = S + h_temp * (2383/4100 * k0 + 0.0 * k1 + 0.0 * k2 + (-341/164) * k3 + 4496/1023 * k4 + (-301/82) * k5 + 2133/4100 * k6 + 45/82 * k7 + 45/164 * k8 + 18/41 * k9);
        [k10, a10] = FirstOrderODE(t10, S10, N_prime, missile_saturation, T1);
        t11        = t;
        S11        = S + h_temp * (3/205 * k0 + 0.0 * k1 + 0.0 * k2 + 0.0 * k3 + 0.0 * k4 + (-6/41) * k5 + (-3/205) * k6 + (-3/41) * k7 + 3/41 * k8 + 6/41 * k9 + 0.0 * k10);
        [k11, ~]   = FirstOrderODE(t11, S11, N_prime, missile_saturation, T1);
        t12        = t + 1 * h_temp;
        S12        = S + h_temp * ((-1777/4100) * k0 + 0.0 * k1 + 0.0 * k2 + (-341/164) * k3 + (4496/1025) * k4 + (-289/82) * k5 + 2193/4100 * k6 + 51/82 * k7 + 33/164 * k8 + 12/41 * k9 + 0.0 * k10 + 1 * k11);
        [k12, ~]   = FirstOrderODE(t12, S12, N_prime, missile_saturation, T1);
        ODEF_temp    = 41/840 * k0 + 34/105 * k5 + 9/35 * k6 + 9/35 * k7 + 9/280 * k8 + 9/280 * k9 + 41/840 * k10;
        ODEFhat_temp = 34/105 * k5 + 9/35 * k6 + 9/35 * k7 + 9/280 * k8 + 9/280 * k9 + 41/840 * k11 + 41/840 * k12;
        TE = h_temp * (ODEF_temp - ODEFhat_temp);
        S_temp = S + h_temp * (ODEF_temp);
        a_temp = norm(41/840 * a0 + 34/105 * a5 + 9/35 * a6 + 9/35 * a7 + 9/280 * a8 + 9/280 * a9 + 41/840 * a10);
        if all(S_temp(:) == 0.0) == true
            threshold = epsilon_absolute / t_guess;
        else
            mag = norm(S_temp);
            threshold = min(((epsilon_relative * mag) / t_guess), epsilon_absolute); 
        end
        TE = min_error + norm(TE);
        if TE == 0.0
            Fscale = Fmax;
            break
        end
        Fscale = beta * (threshold / TE)^(1 / 8);
        Fscale = min(max(Fscale, Fmin), Fmax);
        if TE <= threshold
            break
        end
        j = j + 1;
        if j >= attemptspertimestep
            break
        end
        h_temp = Fscale * h_temp;
    end
    if h_temp < h_min
        h = h_min;
        S = S + h * (ODEF_temp);
    else
        h = h_temp;
        S = S_temp;
    end
    t_next = t + h;
    t = t_next;
    a = a_temp;
    i = i + 1;
    h_next = h * Fscale;
    if h_next < h_min
        h = h_min;
    else
        h = h_next;
    end
end