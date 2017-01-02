function v = avgradius( x )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
summ = 0;

for i = 1:length(x)
    summ = summ + sup(x(i)) - inf(x(i));
end

v = summ/length(x);

end

