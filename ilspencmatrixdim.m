function [ m, n ] = ilspencmatrixdim( A )
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
%  A ... represenation of matrix A
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  m ... integer - number of rows
%  n ... integer - number of columns
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

        m = A{1}(1);
        n = A{1}(2);
end

