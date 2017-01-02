function iv = ilspenchbr( A, b, ip )
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Hansen-Bliek-Rohn bounds to parametric interval systems.
%  Hladik (2012).
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
%ENDDOC====================================================================

% Inicialization of general variables.
dimensions = ilspencmatrixdim(A);
M = intval(zeros(dimensions));
C = intval(zeros(dimensions(1),1));
I = eye(dimensions);
iv = intval(zeros(dimensions(1),1));

radiusVector = ilspencradius(ip);

% Precondition matrix.
Acenter = ilspencmatrixcenter(A,ip);
Acenterinv = inv(Acenter);

% approximate center of the solutioin
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

% Summation in 'x upper-index zero' definition from Theorem 5.
parfor k = 1:length(ip)
    
    if k <= numparb
        bk = ilspencgetbk(b1, b{k+1})
        C = C + radiusVector(k)*abs(Acenterinv*intval(bk));
    end

end

% x upper-index zero from Theorem 5.
x0 = verifylss(M0,abs(x) + C);

%M-asteriks from Theorem 5.
M0 = I - M;
M2 = inv(M0);

M2diag = diag(M2);

u = x0 + (x - abs(x)).*M2diag;
d = -x0 + (x + abs(x)).*M2diag;
coef = (1/(2*M2diag - 1))';

iv = hull(min(inf(d),inf(coef.*d)),max(sup(u),sup(coef.*u)));


end

