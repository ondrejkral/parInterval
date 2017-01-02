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

[m,n] = ilspencmatrixdim(A);
Ar = intval(zeros(m,n));
br = intval(zeros(m,1));

for k=1:length(ip)
    Ar = Ar + ilspencgetak(A,k)*ip(k);
    br = br + ilspencgetbk(b,k)*ip(k);
end
end
