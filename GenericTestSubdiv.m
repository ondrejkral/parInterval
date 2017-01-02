function GenericTestSubdiv(model, matrixType)
clc;
format long infsup;
rng(23101993, 'twister');
global dataModel;
dataModel = model;
iterations = 10;
radiusSample = [0.05, 0.1, 0.5, 1];
matrixDimSample = [5, 10, 15, 20, 25, 50, 100];
methodsCount = 7;
% Active methods set to none.
am = ones(methodsCount,1);

%coefMagMultSample = [5, 10, 15 20, 25];
coefmagnitude = 1; 

fileName = strcat('test-',datestr(datetime(),'dd-mm-yy-HH-MM'),model,matrixType,'subdiv','.txt');
fileLoc = 'C:\Users\ondre\Documents\MATLAB\';
fileAddr = strcat(fileLoc,fileName);

cFileName = strcat('comp-',datestr(datetime(),'dd-mm-yy-HH-MM'),model,matrixType,'.txt');
cFileAddr =  strcat(fileLoc,cFileName);
for l = 1:length(matrixDimSample)
     matrixdim = matrixDimSample(l);

     disp(strcat('DIM : ',num2str(matrixdim)));
        for c = 1:length(radiusSample)
            radius = radiusSample(c);
            fileID = fopen(fileAddr,'a');
            %fileID2 = fopen(cFileAddr,'a');
            
            t = zeros(1,methodsCount); 
            r = zeros(1,methodsCount); 
            %comp = zeros(1,methodsCount);
            
            disp(strcat('RAD : ',num2str(radius)));
            i = 0;
           while 1
                warning off;
                switch(matrixType)
                    case 'Sym'
                        [A, b, p] = SymetricSystem(coefmagnitude, matrixdim, radius);
                    case 'Toeplitz'                       
                        [A, b, p] = toeplitzsystem(coefmagnitude,radius,matrixdim);
                    otherwise
                        % invalid matrixType
                end
                warning on;
                if (verifyinvert(A, b, p))
                 
                   
                    
                    if am(1) == 1
                        %Subdivision
                        tic;
                        x1 = ilspenciterate(A,b,p,1, 'MAX');
                        t(1) = t(1) + toc;
                        r(1) = r(1) + avgradius(x1);
                    end
                    
                    if am(2) == 1
                        %Subdivision
                        tic;
                        x2 = ilspenciterate(A,b,p,2, 'MAX');
                        t(2) = t(2) + toc;
                        r(2) = r(2) + avgradius(x2);
                    end
                    
                    if am(3) == 1
                        %Subdivision
                        tic;
                        x3 = ilspenciterate(A,b,p,3, 'MAX');
                        t(3) = t(3) + toc;
                        r(3) = r(3) + avgradius(x3);
                    end
                    
                    if am(4) == 1
                        %Subdivision
                        tic;
                        x4 = ilspenciterate(A,b,p,4, 'MAX');
                        t(4) = t(4) + toc;
                        r(4) = r(4) + avgradius(x4);
                    end
                    
                    if am(5) == 1
                        %Subdivision
                        tic;
                        x5 = ilspenciterate(A,b,p,5, 'MAX');
                        t(5) = t(5) + toc;
                        r(5) = r(5) + avgradius(x5);
                    end
                    
                    if am(6) == 1
                        %Subdivision
                        tic;
                        x6 = ilspenciterate(A,b,p,6, 'MAX');
                        t(6) = t(6) + toc;
                        r(6) = r(6) + avgradius(x6);
                    end
                    
                    if am(7) == 1
                        %Subdivision
                        tic;
                        x7 = ilspenciterate(A,b,p,7, 'MAX');
                        t(7) = t(7) + toc;
                        r(7) = r(7) + avgradius(x7);
                    end      
                   
                % display number of current iteration
                    i = i +1;
                    disp(i); 
                    if i == iterations 
                        break
                    end
                end
               
            end
            
            % compute avarage time, radius and completeness over iterations
            resulttime = zeros(1,methodsCount);
            for z = 1:methodsCount
                resulttime(z) = t(z)/iterations;
            end

            resultradius = zeros(1,methodsCount);
            for z = 1:methodsCount
                resultradius(z) =r(z)/iterations;
            end
                
            
            % write data to file
            for a = 1:methodsCount
                if am(a) == 1
                    fprintf(fileID,'%.8f %.8f %d %d %d %d \n',resultradius(a),resulttime(a),matrixdim, radius, a);
                end
            end
%             % write completness to file
%              for a = 1:methodsCount
%                     fprintf(fileID2,'%.8f %d \n',resultcomp(a), a);
%             end
            %fclose(fileID2);
            fclose(fileID);
        end
end
end