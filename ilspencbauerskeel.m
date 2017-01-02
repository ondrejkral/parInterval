function iv = ilspencbauerskeel( A, b, ip )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Bauer-Skeel bounds to parametric interval systems.
%  Parametric generalization by Hladik (2012).
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
%  iv ... interval vector - an enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

% Inicialization of general variables.
dimensions = ilspencmatrixdim(A);
M = intval(zeros(dimensions));
C = intval(zeros(dimensions(1),1));
I = eye(dimensions);
radiusVector = rad(ip);

% precondition matrix
Acenter = ilspencmatrixcenter(A,ip);

% aproximate solution
Acenterinv = inv(Acenter);
x = intval(Acenterinv*ilspencbcenter(b ,ip));

% Meta-data cells
    A1 = A{1};
    b1 = b{1};
    numparA = A1(4);
    numparb = b1(3);


% Matrix M from Theorem 4.
parfor k = 1:length(ip)
    
    if k <= numparA
        Ak = ilspencgetak(A1, A{k+1});
        M = M + radiusVector(k)*abs(Acenterinv*intval(Ak));
    end

end

% Summation in interval enclosure from Theorem 4.
parfor k = 1:length(ip)
    
    if k <= numparA
        Ak = ilspencgetak(A1, A{k+1});
    else
        Ak = 0;
    end
    
    if k <= numparb
        bk = ilspencgetbk(b1, b{k+1})
    else
        bk = 0;
    end
    
    C = C + radiusVector(k)*abs(Acenterinv*(Ak*x - bk));
end

s = verifylss(I - M,C);

iv = hull(x - s, x + s);
end

