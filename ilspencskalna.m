function iv = ilspencskalna( iA, ib )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Sklana's direct method for solving interval systems.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  iA ... interval matrix - matrix A of interval system
%  ib ... interval vector - vector b of interval system
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iv ... interval vector - an enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

I = eye(dim(iA));
R = mag(I - iA);

iv = infsup(-1,1)*(verifylss((I - R),mag(ib)));

end

