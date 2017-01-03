function iresult = ilspenciterate( A, b, ip, iterations, option )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Divide parameters' space into several sub-spaces and
%  solve them independently.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A ... represenation of matrix A
%  b ... representation of vector b
%  ip ... interval vector - parameters
%  iterations - integer - number of divisions of parameters' space
%  option - string - passing to the class divisor
%                       'MAX' - largest interval first
%                       'RND' - at random
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iresult ... interval vector - an enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================
iresult = intval(NaN);

if iterations == 0   
    iresult = ilspencresidual(A,b,ip, 'SKALNA');  
    
elseif iterations > 0
    
    newp = ilspencnextpar(1,ip,iterations, option);
    iresult = ilspencresidual(A,b,newp, 'SKALNA'); 
  
    % Iteratively union results.
    parfor i = 2:(2^iterations)
        
        % Getting next subspace.
        newp = ilspencnextpar(i,ip, iterations, option);

        v1 = ilspencresidual(A,b,newp, 'SKALNA'); 
        iresult = hull(iresult,v1);

    end
end
end
