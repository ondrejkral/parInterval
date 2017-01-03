function iV = ilspencgetak(A1, Ak,k)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Returning matrix of coeficients for k-th parameter.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A1 ... first cell block with meta-data
%  Ak ... k-th cell block in representation of matrix A
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iV ... double/interval matrix - coeficients of parameter - A^k matrix
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%  Additional layer for testing repesentation of parameters' linear 
%  dependencies. Because of MATLAB's lazy copying this shouldn't be
%  significant performance overhead.
%
%ENDDOC====================================================================

% Obtaining matrix dimension and other meta-data.
m = A1(1);
n = A1(2);

par = A1(3);
% par == 0: double values
% par == 1: intval values

% number of defined parameters;
%maxpar = A1(4);

switch(par)
    case 0
        iV = zeros(m,n);
        %if k <= maxpar
            for i = 1:length(Ak(1,:))
                column = Ak(:,i); 
                iV(column(1),column(2)) = column(3);
            end
        %end

    case 1
        iV = intval(zeros(m,n));
        %if k <= maxpar
            for i = 1:length(Ak(1,:))
                column = Ak(:,i); 
                % Must indexing with point value.
                iV(column(1).inf,column(2).inf) = column(3);
            end
        %end
    otherwise
        disp('Invalid parameter in data representation.')
        iV = intval(NaN);
end
end

