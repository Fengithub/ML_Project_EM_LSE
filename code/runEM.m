% load data
data = importdata('../Data/wpbc.newdata.txt');
D = data.data;

% preprocessing
[ M, N, T, R ] = splitData( D );
[m, k] = size(M);
[n, k] = size(N);

% train Gamma GLM
num_iter = 5;
[ W, E, R1, replace_rate ] = EM( M, N, T, R, num_iter );

close all
% plot error
figure()
plot((0:num_iter)', E);
xlabel('# iterations')
ylabel('Error')
legend('All', 'Recur', 'Non-recur')

% plot replacement rate
figure()
plot((0:num_iter)', replace_rate);
xlabel('# iterations')
ylabel('% replacement')