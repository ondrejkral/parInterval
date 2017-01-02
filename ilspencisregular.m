function v = ilspencisregular( A, ip)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Checking regularity of parametric matrix.
%  Returns 1 if matrix is regular. 0 when cannot tell.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A ... represenation of matrix A
%  ip ... interval vector - parameters
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  v ... boolean - 0/1
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================v

v = 0;

% Radius vector of parametericvector.
radiusVector = rad(ip);

% Midpoint matrix inverse.
Acenter = ilspencmatrixcenter(A,ip);

% Checking invertibility
dim = size(Acenter);
if ~(rank(Acenter) == dim(1))
    return;
end

Acenterinv = inv(Acenter);

% Init.
[m, n, numparA] = ilspencmatrixdim(A);
M = intval(zeros(m,n));

% Meta-cell
A1 = A{1};

% "Relaxing". See that matrix is non-negative.
parfor k = 1:length(ip)
    if k <= numparA
    M = M + radiusVector(k)*abs(Acenterinv*intval(ilspencgetak(A1,A{k+1})));
    end
end

% Application of theorems in section 1.4 of the thesis.
if verspectrad(sup(M)) < 1 
    v = 1;
else 
    v = 0;
end
end

