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
%ENDDOC====================================================================

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
dimensions = ilspencmatrixdim(A);
M = intval(zeros(dimensions));

% "Relaxing". See that matrix is non-negative.
for k = 1:length(ip)
   M = M + radiusVector(k)*abs(Acenterinv*intval(ilspencgetak(A{1},A{k+1})));
end

% Application of theorems in section 1.4 of the thesis.
if verspectrad(sup(M)) < 1 
    v = 1;
else 
    v = 0;
end
end

