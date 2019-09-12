
% load prognostic data
data = importdata('../Data/wpbc.newdata.txt');
D = data.data;
X = D(:, 4:end);
outcome = logical(D(:, 2));
time = D(:, 3);
[l, k] = size(D);
G = [(1:l)' (1:l)' D(:,1) ];

% cellfun(@(x)strcmp(x,'2'), A)
ERR = zeros(l, 3);
t  = zeros(size(D,1), 1);
for i = 1:l
    if outcome(i)
        fprintf('Sample: #%d:\n', i)
        % preprocessing
        [ M, N, T, R ] = splitData( D, G, [1:i-1 i+1:l]);

        % train EM
        W = glmfit( M(:,2:end), log(T), 'normal' );
        %W = glmfit( [M(:,2:end); N(:,2:end)], [log(T); log(R)], 'normal' );

        % validation error
        [ M2, N2, T2, R2 ] = splitData( D, G, i);
        [err, err1, err2] = calcError2( M2, N2, T2, R2, W, 1, 'log' );

        x = [1 D(i, 4:end)];
        t(i) = exp(x * W);
        ERR(i,:) = [err, err1, err2];
        fprintf('\n')
    end
end

outTable = [G(outcome==1, [3,2]), D(outcome==1, 2:3), t(outcome==1)];
ERR = calcError(t(outcome==1), D(outcome==1, 3), D(outcome==1, 2), 1);
save('leave1out.lm0.mat', 'outTable', 'W', 'ERR')