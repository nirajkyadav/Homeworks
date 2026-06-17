% Importing the data
filename = 'Northridge Record.txt';
data = readtable(filename);

% Time and acceleration
t = data{:, 1};
acc = data{:, 2};

dt_original = 0.02; % Original time step = 0.02 sec
dt_target1 = 0.01;  % Target time step   = 0.01 sec
dt_target2 = 0.005; % Target time step   = 0.005 sec

% Ratio = f_target / f_original = (1/dt_target) / (1/dt_original) = dt_original / dt_target

% For dt = 0.01
[p1, q1] = rat(dt_original / dt_target1);  % Resampling ratio
acc1 = resample(acc, p1, q1);

t1 = (0:length(acc1)-1)' * dt_target1;     % Time column vector for 0.01
output1 = table(t1, acc1, 'VariableNames', {'time_sec', 'Acc_m_s2'});
% writetable(output1, 'Northridge_0_01.txt', 'Delimiter', '\t');

% For dt = 0.005
[p2, q2] = rat(dt_original / dt_target2);  % Resampling ratio
acc2 = resample(acc, p2, q2);

t2 = (0:length(acc2)-1)' * dt_target2;     % Time column vector for 0.005
output2 = table(t2, acc2, 'VariableNames', {'time_sec', 'Acc_m_s2'});
% writetable(output2, 'Northridge_0_005.txt', 'Delimiter', '\t');


figure;
plot(t, acc, 'k-', 'LineWidth', 1.5);
hold on;

plot(t1, acc1, 'r--', 'LineWidth', 1);
plot(t2, acc2, 'b:', 'LineWidth', 1);

grid on;
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');

legend('Original \Deltat = 0.02 s', ...
    'Resampled \Deltat = 0.01 s', ...
    'Resampled \Deltat = 0.005 s');

title('Northridge Record Original and Resampled Comparison');