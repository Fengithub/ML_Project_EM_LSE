data = importdata('../Result/error.txt');
D = data.data;
p = 1-bsxfun(@rdivide, D, D(1, :));