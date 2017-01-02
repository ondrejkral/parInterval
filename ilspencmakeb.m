function V = ilspencmakeb(b)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Returns a vector in current representation and with par = 1.
%  Depends on which representation is active (3D, cell, ...), returns an
%  interval vector b in that format.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  b ... interval vector - simple right side vector b without dependencies
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  V ... right side vector b in active representation
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%ENDDOC====================================================================
switch(dataModel)
    case 'cell'
        V{1} = [length(b); 1;];
        V{2} = intval(zeros(2,length(b)));
        for i = 1:length(b)
            V{2}(:,i) = [intval(i); b(i);];
        end
                
    case '3D'
        V = b;
    otherwise
        disp('Invalid value in dataModel.')
        V = interval(NaN);
end
end

