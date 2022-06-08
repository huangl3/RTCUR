function time = time_avg(times)
time = sum(times.*(times>0),2);
time = mean(time,3);
end