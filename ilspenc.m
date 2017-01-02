function iv = ilspenc( A, b, ip, option )

%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Main method for solving interval systems with linear dependencies.
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
%  iv ... output interval vector - an enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%  ILSPENCGETAK, ILSPENCGETBK, ILSPENCMATRIXDIM, ILSPENCVECTORDIM
%  Reimplement these in case of testing new representation 
%  of dependencies.
%  
%  FAST and DEFAULT same methods.
%
%ENDDOC====================================================================

iv = intval(NaN);

% Test of dataModel inicialization.
global dataModel;
if isempty(dataModel)
    dataModel = '3D';
end
if ~(strcmp(dataModel,'3D') || strcmp(dataModel,'cell'))
    disp('Invalid value of global variable dataModel.')
    return;
end

% Testing number of arguments.
if (nargin < 3 || nargin > 4)
    disp('Invalid number of arguments.') 
    return;
end
if (nargin == 3)
    option = 'DEFAULT';
end

% Testing regularity.
if ~ilspencisregular(A,ip)
    disp('Matrix is not regular.')
    return;
end

% Choosing method.
switch option
    case 'TIGHTEST'
        iv = ilspencmono(A,b,ip, 'NOIMPROVE');
    case 'TIGHT'
        x = ilspencbauerskeel(A,b,ip);
        iv = ilspencbsref(A,b,ip,x);
    case 'FAST'
        iv = ilspencresidual(A,b,ip, 'SKALNA');
    case 'DEFAULT'
        iv = ilspencresidual(A,b,ip, 'SKALNA');

    otherwise
        disp('ILSPENC: Invalid option argument.')
        return;
end
end

