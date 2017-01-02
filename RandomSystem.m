function [A,b,p] = RandomSystem(parCount, dimension, R)
% Number of parameters needed.
parametersCount = parCount + 1;
% Random parameters center.
random = -10 + 20*rand(parametersCount,1);
% parameters vector
p = midrad(random,R);
p(1) = intval('1');
% random vector
randomBCenter = -10 + 20*rand(1,dimension);

% number of occurencies of each parameter
r = randi([1 10],parCount,1);
 
global dataModel;

switch(dataModel)
    case 'cell'
        A{1} = [dimension; dimension; 0;];
        A{2} = zeros(3,dimension);
        for i = 1:dimension
            A{2}(:,i) = [i; i; 100*dimension];
        end

        
        for i = 1:parCount
            randomCoefs = -10 + 20*rand(r(i),1);
            [x,y] = RemDup(randi([1 dimension],r(i),1),randi([1 dimension],r(i),1));
            A{i + 2} = zeros(3,r(i));
            for j = 1:r(i)
                A{i + 2}(:,j) = [x(j); y(j); randomCoefs(j);];
            end
        end
       
        b{1} = [dimension; 0];        
        for i = 1:dimension
            b{2}(:,i)= [i; randomBCenter(i)];
        end
        
    case '3D'        
        A = zeros(dimension, dimension, parametersCount);
        
        % Making matrix diagonally dominant for better inversing.
        A(:,:,1) = 100*dimension*eye(dimension) ;
        
        for i = 1:parCount
            randomCoefs = -10 + 20*rand(r(i),1);
             [x,y] = RemDup(randi([1 dimension],r(i),1),randi([1 dimension],r(i),1));
            for j = 1:length(x)
                A(x(j),y(j),i+1) = randomCoefs(j);
            end
        end
        
        b = zeros(dimension,parametersCount);
        b(:,1) = randomBCenter;
    otherwise
        %invalid model
end
end
