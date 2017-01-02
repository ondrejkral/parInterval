function iv = ilspenchbrref( A,b,ip,ix)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Refinement of Hans-Bliek-Rohn generalization by Hladik (2012).
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A ... represenation of matrix A
%  b ... representation of vector b
%  ip ... interval vector - parameters
%  ix ... interval vector - enclosure of the solution to refine
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iv ... interval vector - refined enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

% Inicialization of general variables.
[m, n, numparA] = ilspencmatrixdim(A);
[~, numparb] = ilspencbdim(b);
Y = intval(zeros(m,n)); y = intval(zeros(m,1));
Z = intval(zeros(m,n)); z = intval(zeros(m,1));
I = eye(m,n);

% Radius vecotr of parametric vecotr.
radiusvector = ilspencradius(ip);

% Precondition matrix.
Acenterinv = inv(ilspencmatrixcenter(A,ip));

% Approximate solution. 
x1 = Acenterinv*ilspencbcenter(b ,ip);

for k = 1:length(ip)
    if k <= numparA
        Ak = ilspencgetak(A1, A{k+1});
    else
        Ak = 0;
    end
    
    if k <= numparb
        bk = ilspencgetbk(b1, b{k+1});
    else
        bk = 0;
    end
    
    a = Acenterinv*(Ak*ix - bk);
    radiusK = radiusvector(k);
    
    parfor j = 1:dimensions(1)
        if inf(a(j)) >= 0
            Y(j,:) = Y(j,:) + radiusK*Acenterinv(j,:)*Ak;
            y(j) = y(j) + radiusK*Acenterinv(j,:)*bk;
        elseif sup(a(j)) <= 0
            Y(j,:) = Y(j,:) - radiusK*Acenterinv(j,:)*Ak;
            y(j) = y(j) - radiusK*Acenterinv(j,:)*bk;
        else
            Z(j,:) = Z(j,:) + radiusK*abs(Acenterinv(j,:)*Ak);
            z(j) = z(j) + radiusK*abs(Acenterinv(j,:)*bk);
        end
    end
end

% M-asteriks from refinement algorithm 2.
M0 = I - abs(Y) - Z;
M = inv(M0);

% 'x upper-index zero' from refinement algorithm 2.
x0 = verifylss(M0,abs(x1) - y + z);

M2diag = diag(M);

u = x0 + (x - abs(x)).*M2diag;
d = -x0 + (x + abs(x)).*M2diag;
coef = (1/(2*M2diag - 1))';

iv = hull(min(inf(d),inf(coef.*d)),max(sup(u),sup(coef.*u)));
end

