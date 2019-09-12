function [ err, err1, err2 ] = calcError( prediction, time, outcome, l )
%CALCERROR Summary of this function goes here
%   Detailed explanation goes here

% recur - only if late cases will be counted
err1 = sqrt(sum((prediction - time).^2.*(prediction > time | prediction < time-l).*outcome)) / sum(outcome) ;
% nonrecur - only if early cases will be counted
err2 = sqrt(sum((time - prediction).^2.*(prediction < time).*(1-outcome))) / sum(1-outcome);
% all
err = (err1 * sum(outcome) + err2 * sum(1-outcome)) / length(outcome);

fprintf('All: %4.3f, nonrecur: %4.3f, recur: %4.3f\n', err, err2, err1)
end

