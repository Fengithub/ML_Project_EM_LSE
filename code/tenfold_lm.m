pool = 0;
train = 1;

% load group info
data = importdata('../Data/TenFoldCV_SampleIdx.txt');
G = data.data;

% load prognostic data
data = importdata('../Data/wpbc.newdata.txt');
D = data.data;

% cellfun(@(x)strcmp(x,'2'), A)
if train
    num_Group = max(G(:,2));
    ERR = zeros(10, 3);
    W = zeros(num_Group, size(D, 2) - 2);
    for i = 1:num_Group
        fprintf('Group: #%d:\n', i)
        % preprocessing
        [ M, N, T, R ] = splitData( D, G, [1:i-1 i+1:10]);
        % train SVM
        w = glmfit( M(:,2:end), log(T), 'normal' );
        %w = glmfit( [M(:,2:end); N(:,2:end)], [log(T); log(R)], 'normal' );
        W(i, :) = w';

        % validation error
        [ M2, N2, T2, R2 ] = splitData( D, G, i);
        [err, err1, err2] = calcError2( M2, N2, T2, R2, w, 1, 'log');
        ERR(i,:) = [err, err1, err2];
        fprintf('\n')
    end
    
    t  = zeros(size(D,1), 1);
    % predict
    for i = 1:size(D, 1)
        g = G(i, 2);
        w = W(g, :);
        x = [1 D(i, 4:end)];
        t(i) = min(exp(x * w'), 130);
    end
else
    % read trained params
    data = load('tenfold.lm.mat');
    t = data.outTable(:, 5);
end


outTable = [G(:, [3,2]), D(:, 2:3), t];
ERR = calcError(t, D(:, 3), D(:, 2), 1);
save('tenfold.lm.mat', 'outTable', 'W', 'ERR')