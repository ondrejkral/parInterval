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

%ILSPENCHBRREF Refinement of Hans-Bliek-Rohn generalization.
% Refinement by Hladik (2012).

% p = vector of interval values in which parameters lies
% x = enclosure obtained from Hans-Bliek-Rohn method

% Inicialization of general variables.
dimensions = ilspencmatrixdim(A);
Y = intval(zeros(dimensions)); y = intval(zeros(dimensions(1),1));
Z = intval(zeros(dimensions)); z = intval(zeros(dimensions(1),1));
I = eye(dimensions);

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

