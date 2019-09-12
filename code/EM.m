function [ W, E, R1, replace_rate ] = EM( M, N, T, R, num_iter, train_on_all )
%EM Summary of this function goes here
%   T is the time to recur (TTR) and R is the disease free time (DFS).

if nargin < 5
    num_iter = 5;
end

if nargin < 6
    train_on_all = 0;
end

R1 = R;
E = zeros(num_iter+1, 3);
replace_rate = zeros(num_iter+1, 1);
% train Gamma GLM
if train_on_all
    W = glmfit( [M(:,2:end); N(:,2:end)], [log(T); log(R)], 'normal' );
else
    W = glmfit( M(:,2:end), log(T), 'normal' );
end
%W = glmfit( M(:,2:end), T, 'Gamma', 'link', 'log' );
% initial error
[err, err1, err2] = calcError2( M, N, T, R, W, 1, 'log' );
E(1, :) = [err err1 err2];
replace_rate(1) = 0;

% predict on non-recur
for i = 1:num_iter
    y = exp(N * W);
    replace_rate(i+1) = sum(y >= R)/length(y);
    R1 = y .* (y >= R) + R1 .* (y < R);
    W = glmfit( [M(:,2:end); N(:,2:end)], [log(T); log(R1)], 'normal' );
    %W = glmfit( [M(:,2:end); N(:,2:end)], [T; R1], 'Gamma', 'link', 'log' );
    % training error
    [err, err1, err2] = calcError2( M, N, T, R, W, 1, 'log' );
    E(i+1, :) = [err err1 err2];
end

end

