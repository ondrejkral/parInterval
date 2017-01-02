function [ m, n, pA] = ilspencmatrixdim( A )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Returns matrix dimensions.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A ... representation of matrix A
%  b ... representation of vector b
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  m ... integer - number of rows
%  n ... integer - number of columns
%  pA ... number of parameters in A
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

        m = A{1}(1);
        n = A{1}(2);
        
        pA = A{1}(4);
end

