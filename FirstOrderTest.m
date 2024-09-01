N_prime = 3;
missile_saturation = 200;
T1 = 0.1;


t0 = 0;
tf = 5;
r_target = [50; -100; 20];
v_target = [5; -1; 3];
r_missile = [50; 3; -12];
v_missile = [-3; -80; -5];
x_missile = [0; 0; 0];
S0 = [r_target; v_target; r_missile; v_missile; x_missile];
r_fuze = 10^-12;

[target_rmat, missile_rmat, time_mat, acceleration_mat, dv] = FirstOrderODESolver(t0, tf, S0, N_prime, missile_saturation, T1);
dist_mat = target_rmat - missile_rmat;
dist_mat = vecnorm(dist_mat, 1);

plot3(target_rmat(1, :), target_rmat(2, :), target_rmat(3, :))
hold on
plot3(missile_rmat(1, :), missile_rmat(2, :), missile_rmat(3, :), 'red')
figure;

plot(time_mat, dist_mat)
figure;

plot(time_mat, acceleration_mat)


fprintf('dv expenditure: %d m/s', dv)