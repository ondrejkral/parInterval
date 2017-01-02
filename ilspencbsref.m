function iv = ilspencbsref( A,b,ip,ix)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Refinement of Bauer-Skeel generalization by Hladik (2012)
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
%  iv ... interval vector - refined encolsure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC==================================================================== 

% Inicialization of general variables.
[m, n, numparA] = ilspencmatrixdim(A);
Y = intval(zeros(m,n)); y = intval(zeros(m,1));
Z = intval(zeros(m,n)); z = intval(zeros(m,1));
I = eye(m,n);

% Radius vector of parametric vector.
radiusvector = ilspencradius(ip);

% Precondition matrix.
Acenterinv = inv(ilspencmatrixcenter(A,ip));

% x-asterisk from Theorem 4.
x1 = intval(Acenterinv*ilspencbcenter(b ,ip));

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
    parfor j = 1:m
        if inf(a(j)) >= 0
            Y(j,:) = Y(j,:) + radiusK*Acenterinv(j,:)*Ak;
            y(j) = y(j) + radiusK*Acenterinv(j,:)*(Ak*x1 - bk);
        elseif sup(a(j)) <= 0
            Y(j,:) = Y(j,:) - radiusK*Acenterinv(j,:)*Ak;
            y(j) = y(j) - radiusK*Acenterinv(j,:)*(Ak*x1 - bk);
        else
            Z(j,:) = Z(j,:) + radiusK*abs(Acenterinv(j,:)*Ak);
            z(j) = z(j) + radiusK*abs(Acenterinv(j,:)*(Ak*x1 - bk));
        end
    end
end

refinement = verifylss(I - abs(Y) - Z, y + z);

iv = hull(x1 - refinement, x1 + refinement);
end

