function V = ilspencmatrixcenter( A, ip)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Compute matrix at the center of the parameteric vector.
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
%  V ... double matrix - matrix A at the center of the parametric vector
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

% Computing center of parameteric vector.
parameterCenter = mid(ip);

% Allocation of the zero matrix with correct size. 
[m,n] = ilspencmatrixdim(A);
V = zeros(m,n);

% number of parameters in A
numparA = A{1}(4);

% Computing matrix at given center of parameteric vector.
% Not verified, approximate value.
for i = 1:length(ip)
    
    if i <= numparA
        Ak = ilspencgetak(A{1}, A{i+1});
    else
        Ak = 0;
    end
    V = V + Ak*parameterCenter(i);
end
end
