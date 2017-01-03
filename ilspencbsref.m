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
[~, numparb] = ilspencbdim(b);
Y = intval(zeros(m,n)); y = intval(zeros(m,1));
Z = intval(zeros(m,n)); z = intval(zeros(m,1));
I = eye(m,n);

% Radius vector of parametric vector.
radiusvector = ilspencradius(ip);

% Precondition matrix.
Acenterinv = inv(ilspencmatrixcenter(A,ip));

% x-asterisk from Theorem 4.
x1 = intval(Acenterinv*ilspencbcenter(b ,ip));

% Meta-cells
A1 = A{1};
b1 = b{1};

parfor k = 1:length(ip) 
    Yadd = intval(zeros(m,n)); yadd = intval(zeros(m,1));
    Zadd = intval(zeros(m,n)); zadd = intval(zeros(m,1));
    
   if k <= numparA
        Ak = ilspencgetak(A1, A{k+1});
    else
        Ak = zeros(m,n);
    end
    
    if k <= numparb
        bk = ilspencgetbk(b1, b{k+1});
    else
        bk = zeros(m,1);
    end
    
    a = Acenterinv*(Ak*ix - bk);
    radiusK = radiusvector(k);
    for j = 1:m
        if inf(a(j)) >= 0
            Yadd(j,:) = radiusK*Acenterinv(j,:)*Ak;
            yadd(j) = radiusK*Acenterinv(j,:)*(Ak*x1 - bk);
        elseif sup(a(j)) <= 0
            Yadd(j,:) = - radiusK*Acenterinv(j,:)*Ak;
            yadd(j) = - radiusK*Acenterinv(j,:)*(Ak*x1 - bk);
        else
            Zadd(j,:) = radiusK*abs(Acenterinv(j,:)*Ak);
            zadd(j) = radiusK*abs(Acenterinv(j,:)*(Ak*x1 - bk));
        end
    end
    Y = Y + Yadd;
    Z = Z + Zadd;
    y = y + yadd;
    z = z + zadd;
end

refinement = verifylss(I - abs(Y) - Z, y + z);

iv = hull(x1 - refinement, x1 + refinement);
end

