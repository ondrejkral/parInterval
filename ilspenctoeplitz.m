function [ A, newb, inewp ] = ilspenctoeplitz( ip, b)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Creates toeplitz system based on actual data model. Right side - vector
%  b - is without dependencies. Just normal vector.
%--------------------------------------------------------------------------
% .Input parameters.
%
%  b ... double/interval vector - vector b of interval system
%  ip ... interval vector - parameters
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  A ... represenation of matrix A
%  newb ... representation of vector b
%  inewp ... interval vector - parameters
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================

global dataModel;

A = intval(NaN);
newb = intval(NaN);

% Adding p1 interval for future affine compatibility.
inewp = intval(zeros(length(ip)+1,1));
inewp(1) = intval('1');
inewp(2:length(ip)+1) = ip;

ip = inewp;

% Check if parameters' count fit the dimension.
if length(ip) ~= 2*length(b)
    disp('Number of parameters does not fit the dimension!');
    return;
end

parametersCount = length(ip);
matrixmagnitude = length(b);

switch(dataModel)
    case '3D'
        
        % Matrix A representation.
        A = zeros(matrixmagnitude,matrixmagnitude,parametersCount);
        
        c = zeros(matrixmagnitude,1); c(1) = 1;
        r = zeros(matrixmagnitude,1);

        for i = (matrixmagnitude+2):(parametersCount)
            c = circshift(c,1);
            A(:,:,i) = toeplitz(c,r);
        end

        c = zeros(matrixmagnitude,1); c(1) = 1;
        r = zeros(matrixmagnitude,1); r(1) = 1;
        for i = 2:matrixmagnitude+1
            A(:,:,i) = toeplitz(c,r);
            c(1) = 0;
            r = circshift(r,1);
        end
        
        % Vector b representation.
        newb = b;
               
    case 'cell'
        
        % Matrix A representation.
        A{1} = [ matrixmagnitude; matrixmagnitude; 0;];
        A{2} = [1; 1; 0;];
        
        for i = 1:matrixmagnitude
            A{i+2} = zeros(3,matrixmagnitude - i + 1);
            for j = 1:(matrixmagnitude - i + 1)
                A{i+2}(:,j)=[j; (j+i-1); 1;];
            end
        end
        for i = 2:matrixmagnitude
            A{matrixmagnitude + i + 1} = zeros(3,matrixmagnitude - i + 1);
            for j = 1:(matrixmagnitude - i + 1)
                A{matrixmagnitude + i + 1}(:,j)=[(j+i-1); j; 1;];
            end
        end
        
        % Vector b representation.
        newb{1} = [length(b); isintval(b);];
        for j = 1:length(b) 
            newb{2}(:,j) =[ j; b(j);] ;
        end
        
    otherwise
        disp('Invalid value in dataModel.')
            
end                
end