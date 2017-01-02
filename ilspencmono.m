function iv = ilspencmono( A, b, ip, option)
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Method for solving interval linear system with dependencies
%  based on monotonicity approach.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  A ... represenation of matrix A
%  b ... representation of vector b
%  ip ... interval vector - parameters
%  option ... string - how to handle non-monotonous parameters
%                       'NOIMPROVE' - leave them alone
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  iv ... interval vector - an enclosure of the solution
%
%--------------------------------------------------------------------------
% .Implementation details.
% 
%  A template for handling this problem is presented. Different attitudes
%  in each step are possible.
%
%ENDDOC====================================================================

% First, we need to compute initial encolsure x*. 
% We can use some of our other methods.
x = ilspencresidual(A,b,ip,'SKALNA');

% Getting relaxed system for fast method.
[Ar, ~] = ilspencrelax(A,b,ip);

% Checking derivative signs of solution with respect to parameters.
D = intval(zeros(length(x),length(ip)));
for k = 1:length(ip)
    % Setting right side of system to bk - Ak*x, which is constant.
    % Our methods behave differently, depending on data model used.
    % Thus we need to create right side vector with respect to it.
    % db = ilspencmakeb(ilspencgetbk(b,k) - ilspencgetak(A,k)*x);
    % But not in this case, verifylss() can operate with b as interval
    % vector.
    db = ilspencgetbk(b,k) - ilspencgetak(A,k)*x;
    
    % Solving derivate signs system by some of our other methods
    % as proposed in MCM method (Skalna 2008).
    D(:,k) = verifylss(Ar,db);
end

% Now we can comopute tight enclosures for dimensions of x in which 
% is x with respect to p monotonous.
% Several approaches to handling non-monotonous parameters can be applied.
iv = intval(zeros(length(x),1));
switch(option)
    case 'NOIMPROVE'
        % No improvement for non-monotonous components.
        for i= 1:length(x)
            newxi = ilspencmonogetbound(A, b, ip, D(i,:), i, 'SHARP');
            iv(i) = newxi;
        end
    otherwise
        disp('ILSPENCMONO: Invalid option argument.')
        iv = intval(NaN);
end
end


