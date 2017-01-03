function [A,b,p] = SymetricSystem(~, dimension, R)
% Number of parameters needed.
parametersCount = (dimension*dimension - dimension)/2 + dimension;
% idle parameter
parametersCount = parametersCount + 1;
% Random parameters center.
random = -10 + 20*rand(parametersCount,1);
% randomACenter = 10*coefMagnitudeMult + randomACenter;
% parameters vector
p = midrad(random,R);
p(1) = intval('1');
% random vector
randomBCenter = -10 + 20*rand(1,dimension);

A{1} = [dimension; dimension; 0; parametersCount];
b{1} = [dimension; 0; 1];

for i = 1:dimension
    b{2}(:,i)= [i; randomBCenter(i)];
end

A{2} = [1; 1; 0];
parameterIndex = 2;
for j = 1:dimension
    for i = 1:j
        if (i == j)
            A{parameterIndex+1} = [i; i; 1;]; 
            p(parameterIndex) = p(parameterIndex) + 10*dimension;
        else
            A{parameterIndex+1} = [i, j; j, i; 1, 1;];
        end
        parameterIndex = parameterIndex + 1;
    end
end
end
