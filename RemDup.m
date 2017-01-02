function [ newx, newy ] = RemDup( x,y )
%RemDup Removes duplicated tuples (x,y) from x and y.

duplicate = zeros(length(x),1);
for i = 1:(length(x)-1)
    actx = x(i); acty = y(i);
    for j = i+1:length(x)
        if (x(j) == actx) && (y(j) == acty)
            duplicate(j) = 1;
        end
    end
end

newindex = 1;
for i = 1:length(x)
    if ~duplicate(i)
        newx(newindex) = x(i);
        newy(newindex) = y(i);
        newindex = newindex + 1;
    end
end