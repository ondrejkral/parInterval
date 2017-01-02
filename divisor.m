classdef divisor < handle
%BEGINDOC==================================================================
% .Author.
%
%  Ondrej Kral
%
%--------------------------------------------------------------------------
% .Description.
%
%  Object for managing division of parameters' space.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  ipar ... interval vector - vector of parameters
%  count ... integer - number of divisions
%  option ... string - which parameter split first
%                       'MAX' - largest first
%                       'RND' - at random
%
%--------------------------------------------------------------------------
% .Implementation details.
%
% Calling nextPar() on class divisor retruns interval vector of parameters
% for specific sub-space.
%
%ENDDOC====================================================================
    
    properties (SetAccess = private)
        paramvector % parameters
        state % actual subspace, vector of ones and zeros
        divindex % indicies of parameters to divide
    end
    
    methods
        % constructor
        function obj = divisor(ipar, count, option)
            obj.paramvector = ipar;
            obj.divindex = findIndex(ipar, count, option);
            obj.state = zeros(count,1);
        end
        
        % next parameters' vector in subspace
        function v = nextPar(obj)
            v = getPar(obj);
            updateState(obj);
        end
        
        % updating actual subspace, "addition with carry"
        function updateState(obj)
            for i = 1:length(obj.state)
                if obj.state(i) == 1
                   obj.state(i) = 0;
                else
                    obj.state(i) = 1;
                    break;
                end
            end                
        end
        % getting subspace from state zero-one values
        function p = getPar(obj)
            p = obj.paramvector;
            % divide parameter at chosen positions
            for i = 1:length(obj.divindex)
                index = obj.divindex(i);
                % 1 for upper part
                if obj.state(i) == 1
                    p(index) = infsup(mid(p(index)), p(index).sup);
                else
                % 0 for lower part
                    p(index) = infsup(p(index).inf, mid(p(index)));
                end
            end
        end
    end
end

% division rule here
function v = findIndex(par, count, option)
            v = zeros(count,1);
            switch(option)
                % choose interval with maximum width
                case 'MAX'
                    
                    % simulating division progress
                    radius = rad(par);   
                    for i = 1:count
                        [~,index] = max(radius);
                        v(i) = index;
                        radius(index) = radius(index)/2;
                    end
                
                % random choice
                case 'RND'
                     rnd = randperm(length(par));
                     v = rnd(1:count);

                otherwise
                    disp('divisor.findIndex: Invalid option value.')
                    v = NaN;
                    return;
            end
end