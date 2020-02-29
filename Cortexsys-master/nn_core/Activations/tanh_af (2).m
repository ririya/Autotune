classdef tanh_af < handle
    properties(GetAccess = 'public', SetAccess = 'private')
       costType;
       classID;
       alpha;
       rho;
       eps;
    end
    
    properties(GetAccess = 'public', SetAccess = 'public')
       defs;
    end
     
    properties(GetAccess = 'public', Constant = true)
       type = 1;
       name = 'tanh';
    end

    methods
        function obj = tanh_af(varargin) 
            obj.classID = cputime;
            
            if (numel(varargin) == 3)
                obj.costType = varargin{2};
                obj.defs = varargin{1};
                obj.alpha = precision(varargin{3}, obj.defs);
            elseif (numel(varargin) == 2)
                obj.costType = varargin{2};
                obj.defs = varargin{1};
            elseif  (numel(varargin) == 1) 
                obj.defs = varargin{1};
            end  
            % Setup training hyperparameters
            obj.alpha = precision(0.0125, obj.defs); % gradient descent learning rate
            obj.rho = precision(0.95, obj.defs); % AdaDelta hyperparameter (don't generally need to modify)
            obj.eps = precision(1e-7, obj.defs); % AdaDelta hyperparameter (don't generally need to modify)
        end
        
        function [W, b] = initializeWeightsBiases(obj, L_out, L_in, varargin)
            % Bengio, X. Glorot, Understanding the difficulty of training deep feedforward neuralnetworks, AISTATS 2010
            % http://deeplearning.net/tutorial/mlp.html
            if numel(varargin) == 1
                if strcmp(varargin{1}, 'lstm')
                    half_interval = sqrt(6/(L_out+4*L_in));
                end
            else
                half_interval = sqrt(6/(L_out+L_in));
            end
            W = half_interval*(2*rand(L_out, L_in, obj.defs.PRECISION) - 1);
            b = precision(zeros(L_out, 1), obj.defs);
        end
        
        function g = activ(~,z,~)
            g = tanh(z);
        end
        
        function gp = grad(~,a)
           gp = 1-a.^2;
        end
        
        % Gradient at the output layer %
        function gp = ograd(obj,z)
            if (obj.costType == obj.defs.COSTS.SQUARED_ERROR)    
                gp = obj.grad(z); % squared error
            end
        end
        
        function J = cost(~,Y, A, m, t)
            J = squaredErrorCostFun(Y, A, m, t);
        end
    end
end
