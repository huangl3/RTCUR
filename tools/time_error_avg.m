% time_mat and error_mat both have size (k_max, i_max), this function is
% designed to get rid of (-1)'s and generate a proper vector pair 
% (time_avg, error_avg) for plot.
function [time_avg, error_avg] = time_error_avg(time_mat,error_mat,start_iter)
k_max = size(time_mat,1);
time_avg = cumsum(sum(time_mat, 1)./sum(time_mat>0, 1));
time_avg = time_avg - time_avg(start_iter);
error_avg = geo_mean(error_mat);
end_iter = floor(find(time_mat<0, 1)/k_max);
time_avg = time_avg(start_iter:end_iter);
error_avg = error_avg(start_iter:end_iter);
end