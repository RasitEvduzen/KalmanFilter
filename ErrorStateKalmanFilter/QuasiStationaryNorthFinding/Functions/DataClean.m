function [data] = DataClean(data,thmin,thmax)
data(abs(data) > thmax) = mean(data);
data(find((data) < (mean(data)-thmin))) = mean(data);
data(find((data) > (mean(data)+thmin))) = mean(data);
end