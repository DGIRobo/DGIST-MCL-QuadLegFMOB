clear all
clc
close all
system('./run_linux')

filename = './data/data.csv';

T = readtable(filename);
VariableNames = T.Properties.VariableNames;

Arr = table2array(T);
[m,n] = size(Arr);

figure(1);
subplot(2,1,1);
plot(Arr(:,1), Arr(:,2), 'Color', 'b', 'LineWidth', 3);
hold on;
plot(Arr(:,1), Arr(:,4), 'Color', 'r', 'LineWidth', 3);
hold on;
plot(Arr(:,1), Arr(:,3), 'Color', 'm', 'LineWidth', 3);
title("GRF estimation under K_o=2*pi*100, f_{cutoff}=100Hz");
xlim([0 0.7]);xlabel('time (sec)');
ylabel('r-direction Ground Reaction Force (N)');
legend('real', 'MOB', 'FOB')

subplot(2,1,2);
plot(Arr(:,1), Arr(:,4)-Arr(:,2), 'Color', 'g', 'LineWidth', 3);
hold on;
plot(Arr(:,1), Arr(:,3)-Arr(:,2), 'Color', 'c', 'LineWidth', 3);
xlim([0 0.7]);
title("GRF estimation under K_o=2*pi*100, f_{cutoff}=100Hz");
xlabel('time (sec)');
ylabel('r-direction Ground Reaction Force (N)');
legend('MOB error', 'FOB error')

%{
subplot(2,2,3);
plot(Arr(:,1), Arr(:,5), 'Color', 'b', 'LineWidth', 3);
hold on;
plot(Arr(:,1), Arr(:,6), 'Color', 'r', 'LineWidth', 3);
xlim([0 1.2]);
xlabel('time (sec)');
ylabel('r-direction (m)');
legend('reference', 'real')

subplot(2,2,4)
plot(Arr(:,1), Arr(:,7), 'Color', 'b', 'LineWidth', 3);
hold on;
plot(Arr(:,1), Arr(:,8), 'Color', 'r', 'LineWidth', 3);
xlim([0 1.2]);
xlabel('time (sec)');
ylabel('theta-direction (rad)');
legend('reference', 'real')

figure(2);
subplot(2,2,1)
plot(Arr(:,1), Arr(:,9), 'Color', 'r', 'LineWidth', 3);
xlabel('time (sec)');
ylabel('tau_bi_m (Nm)');

subplot(2,2,2)
plot(Arr(:,1), Arr(:,10), 'Color', 'r', 'LineWidth', 3);
xlabel('time (sec)');
ylabel('tau_bi_b (Nm)');

subplot(2,2,3)
plot(Arr(:,1), Arr(:,11), 'Color', 'r', 'LineWidth', 3);
xlabel('time (sec)');
ylabel('qddot_m (Nm)');

subplot(2,2,4)
plot(Arr(:,1), Arr(:,12), 'Color', 'r', 'LineWidth', 3);
xlabel('time (sec)');
ylabel('qddot_b (Nm)');
%}
