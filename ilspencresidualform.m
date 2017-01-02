function [ iAres, ibres] = ilspencresidualform( A, b, ip )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Converts interval linear system with dependencies
%  to residual form.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A ... represenation of matrix A
%  b ... representation of vector b
%  ip ... interval vector - parameters
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iAres ... interval matrix - matrix A in the residual form
%  ibres ... interval vector - vector b in the residual form
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

% Initial allocation.
dimensions = ilspencmatrixdim(A);
iAres = intval(zeros(dimensions));
ibres = intval(zeros(dimensions(1),1));

% Precondition matrix.
Cinv = inv(ilspencmatrixcenter(A,ip));

% "Center" of solution set.
x = Cinv*ilspencbcenter(b ,ip);

% Residual form.
for k = 1:length(ip)
    iAres = iAres + ip(k)*(Cinv*intval(ilspencgetak(A,k)));
end

for k = 1:length(ip)
    ibres = ibres + ip(k)*(Cinv*(ilspencgetbk(b,k) - ilspencgetak(A,k)*intval(x)));
end

end

