function iv = ilspencrump( iA, ib, iterations)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Bounds to interval system by Rump(2010).
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  iA ... interval matrix - matrix A of interval system
%  ib ... interval vector - vector b of interval system
%  iterations ... integer - how many epsilon-inflations
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iv ... interval vector - an enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

iv = intval(NaN);

% Approximate inverse.
R = inv(mid(iA));

%Approximate center of solution set.
x = R*mid(ib);

% Iteration matrix.
C = eye(dim(iA))-R*intval(iA);
Z = R*(ib-iA*intval(x));
X = Z; iter = 0;
while iter < iterations
    iter = iter+1;
    Y = X*infsup(0.9,1.1) + 1e-20*infsup(-1,1); % epsilon-inflation
    X = Z+C*Y; % interval iteration

    if all(in0(X,Y)), iv = x + X; 
        return;
    end
end
end

