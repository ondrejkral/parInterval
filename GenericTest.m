function GenericTest(matrixType)
clc;
format long infsup;
rng(23101993, 'twister');

% properties of testing
iterations = 10;
radiusSample = [0.05, 0.1, 0.5, 1];
matrixDimSample = [5, 10, 15, 20, 25, 50, 100];
methodsCount = 10;

% Active methods set to all.
am = ones(methodsCount,1);

am(8) = 0;
am(9) = 0;
am(10) = 0;

% coefMagMultSample = [5, 10, 15 20, 25];
coefmagnitude = 1; 


% logging results
fileName = strcat('test-',datestr(datetime(),'dd-mm-yy-HH-MM'),matrixType,'.txt');
%fileLoc = 'C:\Users\ondre\Documents\MATLAB\';
 fileLoc = '';
fileAddr = strcat(fileLoc,fileName);

% setting parallel options
% cores = feature('numcores');
cores = 2;

for cN = 1:cores
% maxNumCompThreads(cN);
poolobj = gcp('nocreate');
delete(poolobj);
parpool(cores - cN + 1);

    for l = 1:length(matrixDimSample)
         matrixdim = matrixDimSample(l);

         disp(strcat('DIM : ',num2str(matrixdim)));
            for c = 1:length(radiusSample)
                radius = radiusSample(c);
                fileID = fopen(fileAddr,'a');

                t = zeros(1,methodsCount); 
                r = zeros(1,methodsCount); 

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
                    if (ilspencisregular(A, p))

                        if am(5) == 1
                            %disp('Residual method with RUMP');
                            tic
                            x3 = ilspencresidual(A, b, p, 'RUMP');
                            elapsed = toc;
                            if( ~isnan(x3))
                                t(5) = t(5) + elapsed;
                                r(5) = r(5) + avgradius(x3);
                            else
                                continue;
                            end
                        end

                        if am(6) == 1
                            %disp('Residual method with SKALNA');
                            tic;
                            x4 = ilspencresidual(A, b, p, 'SKALNA');
                            t(6) = t(6) + toc;
                            r(6) = r(6) + avgradius(x4);
                        end

                        if am(1) == 1
                            %disp('Bauer-Skeel enclousure');
                            tic;
                            x1 = ilspencbauerskeel(A, b, p);
                            t(1) = t(1) + toc;
                            r(1) = r(1) + avgradius(x1);
                        end

                        if am(2) == 1
                            %disp('Bauer-Skeel refinement');
                            tic;
                            r1 = ilspencbsref(A, b, p, x1);
                            t(2) = t(2) + toc;
                            r(2) = r(2) + avgradius(r1); 
                        end

                        if am(3) == 1
                            %disp('Hans-Bliek-Rohn enclosure');
                            tic;
                            x2 = ilspenchbr(A, b, p);
                            t(3) = t(3) + toc;
                            r(3) = r(3) + avgradius(x2);
                        end

                        if am(4) == 1
                            %disp('Hans-Bliek-Rohn refinement');
                            tic;
                            r2 = ilspenchbrref(A, b, p, x2);
                            t(4) = t(4) + toc;
                            r(4) = r(4) + avgradius(r2);
                        end

                        if am(7) == 1
                            %disp('Monotonicity');
                            tic;
                            x5 = ilspencmono(A, b, p, 'NOIMPROVE');
                            t(7) = t(7) + toc;
                            r(7) = r(7) + avgradius(x5);
                        end

                        if am(8) == 1
                            %Subdivision
                            tic;
                            x6 = ilspenciterate(A,b,p,5, 'SKALNA');
                            t(8) = t(8) + toc;
                            r(8) = r(8) + avgradius(x6);
                        end
                        if am(9) == 1
                            %Subdivision
                            tic;
                            x7 = ilspenciterate(A,b,p,6, 'SKALNA');
                            t(9) = t(9) + toc;
                            r(9) = r(9) + avgradius(x7);
                        end
                        if am(10) == 1
                            %Subdivision
                            tic;
                            x8 = ilspenciterate(A,b,p,7, 'SKALNA');
                            t(10) = t(10) + toc;
                            r(10) = r(10) + avgradius(x8);
                        end                 

                    % display number of current iteration
                        i = i +1;
                        disp(i); 
                        if i == iterations 
                            break
                        end
                    end

                end

                % compute avarage time, radius over iterations
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
                        fprintf(fileID,'%.8f %.8f %d %d %d %d \n',resultradius(a),resulttime(a),matrixdim, radius, cN, a);
                    end
                end
                fclose(fileID);
            end
    end
end
end