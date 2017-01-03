function iv = ilspencgetbk( b1, bk,k)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Returning vector of coeficients for k-th parameter.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  b1 ... first cell block with meta-data
%  bk ... k-th cell block in representation of vector b
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iv ... double/interval vector - an enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%   Additional layer for testing repesentation of parameters' linear 
%   dependencies. Because of MATLAB's lazy copying this shouldn't be
%   significant performance overhead.
%   Be careful, when k is out of index, scalar 0 is returned instead
%   allocating whole vector of zeros due to better performance.
%
%ENDDOC====================================================================

%Obtaining b-vector dimension and other meta-data.
m = b1(1);
par = b1(2);
% par == 0: double values
% par == 1: intval values

% number of defined parameters
%maxpar = b1(3);

switch(par)
    case 0
    iv = zeros(m,1);

    %if k <= maxpar
        for i = 1:length(bk(1,:))
            column = bk(:,i); 
            iv(column(1)) = column(2);
        end
    %end

    case 1
    iv = intval(zeros(m,1));

    %if k <= maxpar
        for i = 1:length(bk(1,:))
            column = bk(:,i);
            % Must indexing with infimum.
            iv(column(1).inf) = column(2);
        end 
    %end

otherwise
    disp('Invalid parameter in data representation.')
    iv = intval(NaN);
end                
end

