function iv = ilspencmonogetbound( A, b, ip, derivatives, i, option)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Compouting bounds of i-th component of the solution with respect
%  to parameters monotonicity.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A ... represenation of matrix A
%  b ... representation of vector b
%  ip ... interval vector - parameters
%  derivatives - interval vector - derivatives of the solution functions
%  with respect to each parameter
%  option - string - a way to compute enclosure of the i-th component
%                   'FAST' - faster but less sharp
%                   'SHARP' - slower but sharper
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iv ... interval vector - an enclosure of the i-th component of the
%  solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%  Implemented as described in monotonicity section in the thesis.
%
%ENDDOC====================================================================

iv = intval(NaN);

% Based on signs of derivates we can reduce some parameters to points.
lowerp = intval(zeros(length(derivatives),1));
upperp = intval(zeros(length(derivatives),1));
for j=1:length(derivatives)
    if derivatives(j) >= 0
        lowerp(j) = ip(j).inf;
        upperp(j) = ip(j).sup;
    elseif derivatives(j) <= 0
        lowerp(j) = ip(j).sup;
        upperp(j) = ip(j).inf;
    elseif derivatives(j) == 0
        lowerp(j) = mid(ip(j));
        upperp(j) = mid(ip(j));
    else
        lowerp(j) = ip(j);
        upperp(j) = ip(j);
    end
end

% We shrink some parameters. Now we deploy some method to solve our new 
% systems. Place for improvements here. Like if we shrink all parameters,
% we can compute bounds directly.

switch(option)
    case 'FAST'
        [Ar,br] = ilspencrelax(A,b,lowerp);
        lowerx = verifylss(Ar,br);

        [Ar,br] = ilspencrelax(A,b,upperp);
        upperx = verifylss(Ar,br);
    case 'SHARP'
        lowerx = ilspencresidual(A,b,lowerp,'SKALNA');
        upperx = ilspencresidual(A,b,upperp,'SKALNA');
    otherwise
        disp('ILSPENCMONOGETBOUND: Invalid option argument.')
end

iv = infsup(lowerx(i).inf,upperx(i).sup);
end