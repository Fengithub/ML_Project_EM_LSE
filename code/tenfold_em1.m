pool = 0;
train = 1;

% load group info
data = importdata('../Data/TenFoldCV_SampleIdx.txt');
G = data.data;

% load prognostic data
data = importdata('../Data/wpbc.newdata.txt');
D = data.data;
X = D(:, 4:end);
outcome = logical(D(:, 2));
time = D(:, 3);

% cellfun(@(x)strcmp(x,'2'), A)
if train
    num_Group = max(G(:,2));
    ERR = zeros(10, 3);
    W = zeros(num_Group, size(D, 2) - 2);
    for i = 1:num_Group
        fprintf('Group: #%d:\n', i)
        % preprocessing
        [ M, N, T, R ] = splitData( D, G, [1:i-1 i+1:10]);
        
        % load weightedLSE
        data = importdata(...
            sprintf('../Result/WeightedLSE/WeightedLSE_iter1_trainCV%d.txt', i));
        Y = data.data;
        outcome = logical(Y(:, 3));
        DFS = Y(outcome == 0, 4);
        pred = Y(outcome == 0, 5);
        R = max(DFS, pred);
        
        % train model
        [ w, E, R1, replace_rate ] = EM( M, N, T, R, 5, 1 );
        W(i, :) = w';

        % validation error
        [ M2, N2, T2, R2 ] = splitData( D, G, i);
        [err, err1, err2] = calcError2( M2, N2, T2, R2, w, 1, 'log' );
        %[err, err1, err2] = calcError2( M2, N2, T2, R2, w, 1, 'log');
        ERR(i,:) = [err, err1, err2];
        fprintf('\n')
    end
    
    t  = zeros(size(D,1), 1);
    % predict
    for i = 1:size(D, 1)
        g = G(i, 2);
        w = W(g, :);
        x = [1 D(i, 4:end)];
        t(i) = exp(x * w');
        %t(i) = x * w';
    end
else
    % read trained params
    data = load('tenfold.em1.mat');
    t = data.outTable(:, 5);
end


outTable = [G(:, [3,2]), D(:, 2:3), t];
ERR = calcError(t, D(:, 3), D(:, 2), 1);
save('tenfold.em1.mat', 'outTable', 'W', 'ERR')