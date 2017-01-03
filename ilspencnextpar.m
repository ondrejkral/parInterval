function inewpar = ilspencnextpar( i, ip, iter, option)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Compute specific parameters' sub-space
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  i ... index of specific sub-space
%  ip ... interval vector - parameters
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  inewpar ... interval vector - new sub-space of parametric vector
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

inewpar = intval(zeros(length(ip),1));

% computing where to split intervals - indices
v = zeros(iter,1);
switch(option)
    % choose interval with maximum width
    case 'MAX'

        % simulating division progress
        radius = rad(ip);   
        for i = 1:count
            [~,index] = max(radius);
            v(i) = index;
            radius(index) = radius(index)/2;
        end

    % random choice
    case 'RND'
         rnd = randperm(length(par));
         v = rnd(1:count);

    otherwise
        disp('ilspencnexpar: Invalid option value.')
        return;
end
% divide parameter at index position, upper or lower part determined
% by binary representation of i

state = decimalToBinaryVector(i,iter);
for j = 1:iter
    index = v(j);
    % 1 for upper part
    if state(j) == 1
        inewpar(index) = infsup(mid(ip(index)), ip(index).sup);
    else
    % 0 for lower part
        inewpar(index) = infsup(ip(index).inf, mid(ip(index)));
    end
end
        
end
