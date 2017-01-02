function iv = ilspencresidual( A, b, ip, option )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Solving parametric equations in residual form.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A ... represenation of matrix A
%  b ... representation of vector b
%  ip ... interval vector - parameters
%  option ... string - a method to solve a system in the resiudal form
%                       'RUMP' - Rump's iterative alghoritm
%                       'SKALNA' - Skalna's direct method
%                       'SIMPLE' - using relaxation and verifylss()
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iv ... interval vector - an enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
% "Center" x must be computed in same way as in ilspencresidualform() 
% to maintain rigorous envelope of solution.
%
%ENDDOC====================================================================

iv = intval(NaN);

% "Center" of solution set.
x = inv(ilspencmatrixcenter(A,ip))*ilspencbcenter(b ,ip);

% Computing residual form.
[Ares, bres] = ilspencresidualform(A, b, ip);

% Choosing a method to solve the residual form.
switch option
    case 'RUMP'
        y = ilspencrump(Ares, bres, 15);
    case 'SKALNA'
        y = ilspencskalna(Ares, bres);
    case 'SIMPLE'
        y = verifylss(Ares, bres);
    otherwise
        disp('ILSPENCRESIDUAL: Invalid option argument.')
        return;
end

% Combining results as proposed in the text.
iv = x + y;
    
end

