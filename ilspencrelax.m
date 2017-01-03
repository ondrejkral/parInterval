function [ Ar, br ] = ilspencrelax( A,b,ip )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Relaxing dependencies to normal interval system.
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
%  Ar ... interval matrix - matrix A without dependencies
%  br ... interval vector - vector b without dependencies
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

% init
[m,n, numparA] = ilspencmatrixdim(A);
[~, numparb] = ilspencbdim(b);
Ar = intval(zeros(m,n));
br = intval(zeros(m,1));

% Meta-cells
A1 = A{1};
b1 = b{1};

for k=1:length(ip)
    if k <= numparA
    Ar = Ar + ilspencgetak(A1,A{k+1})*ip(k);
    end
    
    if k <= numparb
    br = br + ilspencgetbk(b1,b{k+1})*ip(k);
    end
end
end
