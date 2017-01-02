function v = ilspencbcenter( b, ip )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Computing right side vector b at the center of the parametric vector.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  b ... representation of vector b
%  ip ... interval vector - parameters
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iv ... double vector - b at the center of the parametric vector
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

% Computing center of parameter vector.
parameterCenter = mid(ip);

% Allocation of zero vector with correct size. 
l = ilspencbdim(b);
v = zeros(l,1);

% number of parametrs in b
numparb = b{1}(3);
% Computing vector at given center of parameter vector.
for i = 1:length(ip);
    
    if i <= numparb
        bk = ilspencgetbk(b{1}, b{i+1});
    else
        bk = 0;
    end
    
    v = v + bk*parameterCenter(i);
end
end

