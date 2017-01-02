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
[m, n, numparA] = ilspencmatrixdim(A);
[~, numparb] = ilspencbdim(b);
iAres = intval(zeros(m,n));
ibres = intval(zeros(m,1));

% Precondition matrix.
Cinv = inv(ilspencmatrixcenter(A,ip));

% "Center" of solution set.
x = Cinv*ilspencbcenter(b ,ip);

% Meta-data cells
A1 = A{1};
b1 = b{1};

% Residual form.
parfor k = 1:length(ip)
    if k <= numparA
    iAres = iAres + ip(k)*(Cinv*intval(ilspencgetak(A1,A{k+1})));
    end
end

parfor k = 1:length(ip)
    if k <= numparb
    ibres = ibres + ip(k)*(Cinv*(ilspencgetbk(b1,b{k+1}) - ilspencgetak(A1,A{k+1})*intval(x)));
    end
end

end

