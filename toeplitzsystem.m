function [ A, b, p ] = toeplitzsystem( ~, radius, matrixmagnitude)
% Number of parameters needed.
parametersCount = 2*matrixmagnitude;

% Random center of parameters. Second parameter is diagonal one, so we 
% made it bigger to avoid irregularity.
random = -10 + 20*rand(parametersCount,1);
random(2) = random(2) + 10*matrixmagnitude;

% parameters vector
p = midrad(random,radius);
% First parameter is affine one.
p(1) = intval('1');

% random b vector
randomBCenter = -10 + 20*rand(1,matrixmagnitude);
        
% matrix A representation
A{1} = [ matrixmagnitude; matrixmagnitude; 0; parametersCount;];
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

% vector b representation
b{1} = [matrixmagnitude; 0; 1;];
for j = 1:matrixmagnitude 
    b{2}(:,j) =[ j; randomBCenter(j);] ;
end
for i = 3:(parametersCount+1)
    b{i} = [1; 0;];
end
                         
end