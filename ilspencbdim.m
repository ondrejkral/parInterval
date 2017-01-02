function [l, pb] = ilspencbdim( b )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Returns length of vector b.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  b ... representation of vector b
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  l ... integer - length of b
%  pb ... number of parameters in b
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC==================================================================== 

        l = b{1}(1);
        pb = b{1}(3);
end

