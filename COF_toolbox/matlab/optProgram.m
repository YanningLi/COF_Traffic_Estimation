% Yanning Li, Sep 01, 2015
% This class construct a the optimization program. It sets the linear
% constraints of the optimization program, has multiple functions to set
% the objective function, and solve the program using cplex.


classdef optProgram < handle
    
    properties
        % configuration
        start_time;        
        end_time;
        
        % Decision variable locator, s atruct
        % dv_index.(linkStr) = [start of upstream, start of downstream, start of initial;
        %                  end of upstream, end of downstream, end of initial
        dv_index;
        dv_index_link_max; % the maximum index for links
        
        % there are additionary decision variables 
        % ( the L1 error e=|q1-Rq1|) at merges or diverges 
        % that we would like to add in the decision variable.
        % dv_index.(juncStr) = [start index; end index]
        dv_index_max;    % the maximal index
                
        % a matrix that contains all inequality constraints from all links
        Aineq;
        bineq;
        size_Aineq; % save the size of Aineq
        
        % a matrix that contains all equality constraints for all links
        Aeq;
        beq;
        size_Aeq;   % save the size of Aeq
        
        % upper, lower, and type of each decision variable
        lb;
        ub;
        ctype;
        
        H;  % for defining quadratic objective
        f;  % the linear objective
        
    end
    
    
    methods
        %===============================================================
        % initialize object
        function self = optProgram(~)
            self.start_time = 0;
            
            self.end_time = 0;
            
            self.dv_index = struct;
            self.dv_index_max = 0;
            self.dv_index_link_max = 0;
            
            
            % preallocate memory to speed up
            self.Aineq = zeros(100000,2000);
            self.bineq = zeros(100000,1);
            self.size_Aineq = [0 0];
            
            self.Aeq = zeros(0,0);
            self.beq = zeros(0,1);
            self.size_Aeq = [0 0];
            
            self.lb = zeros(0,1);
            self.ub = zeros(0,1);
            self.ctype = '';
            
            self.f = [];
            self.H = [];
            
        end
        
        
        %===============================================================
        % set configuration for simulation
        function setConfig(self, start_time, end_time)
            
            self.start_time = start_time;
            self.end_time = end_time;
            
        end
        
        
        %===============================================================
        % set Linear constraints
        function setConstraints(self, net, e_max)

            global INDEX_UP INDEX_DOWN INDEX_INI
            
            % set inequality for each link
            for link = net.link_labels'
                
                linkStr = sprintf('link_%d',link);
                
                % set inequality constraints
                constraints_ineq = setIneqConstraints(...
                    net.network_hwy.(linkStr).para_vf,...
                    net.network_hwy.(linkStr).para_w,...
                    net.network_hwy.(linkStr).para_kc,...
                    net.network_hwy.(linkStr).para_km,...
                    net.network_hwy.(linkStr).para_postkm*1000,...
                    self.start_time, self.end_time,...
                    net.network_hwy.(linkStr).BC_us,...
                    net.network_hwy.(linkStr).BC_ds,...
                    net.network_hwy.(linkStr).T_us,...
                    net.network_hwy.(linkStr).T_ds,...
                    net.network_hwy.(linkStr).X_grid_cum,...
                    net.network_hwy.(linkStr).IC, e_max);
                
                % Get the model and data constraints
                link_model_constraints = constraints_ineq.setModelMatrix;
                link_data_constraints = constraints_ineq.setDataMatrix;
                link_constraints = [link_model_constraints; link_data_constraints];
                
                % number of rows to be added in the constraints
                num_row_constraints = size(link_constraints, 1);
                
                % number of col should be same length as decision vairable
                % the last column is the right hand side value Ax<=b
                num_col_constraints = size(link_constraints,2) - 1;
                
                % add the constraints of this link to the full constraints
                % The rest of the element will be filled with 0
                self.Aineq(self.size_Aineq(1)+1:self.size_Aineq(1) + num_row_constraints,...
                    self.size_Aineq(2)+1:self.size_Aineq(2)+ num_col_constraints) ...
                    = link_constraints(:,1:num_col_constraints);
                
                self.bineq(self.size_Aineq(1)+1:self.size_Aineq(1) + num_row_constraints,1)...
                    = link_constraints(:,num_col_constraints+1);
                
                % update Aineq matrix size
                self.size_Aineq = self.size_Aineq + [num_row_constraints, num_col_constraints];
                
                % update index for decision variable
                self.dv_index.(linkStr) = zeros(2,3);
                % add index for upstream boundary flows
                self.dv_index.(linkStr)(1, INDEX_UP) = self.dv_index_max + 1;
                self.dv_index.(linkStr)(2, INDEX_UP) = self.dv_index_max + ...
                   length(net.network_hwy.(linkStr).BC_us);
                self.dv_index_max = self.dv_index_max + length(net.network_hwy.(linkStr).BC_us);
                % add index for downstream boundary flows
                self.dv_index.(linkStr)(1, INDEX_DOWN) = self.dv_index_max + 1;
                self.dv_index.(linkStr)(2, INDEX_DOWN) = self.dv_index_max +...
                   length(net.network_hwy.(linkStr).BC_ds);
                self.dv_index_max = self.dv_index_max + length(net.network_hwy.(linkStr).BC_ds);
                % add index for initial densities
                self.dv_index.(linkStr)(1, INDEX_INI) = self.dv_index_max + 1;
                self.dv_index.(linkStr)(2, INDEX_INI) = self.dv_index_max + ...
                    length(net.network_hwy.(linkStr).IC);
                self.dv_index_max = self.dv_index_max + length(net.network_hwy.(linkStr).IC);
                    
            end
            
            % save the maximal index for the decision variable associated
            % with links (up, down, and initial conditions)
            self.dv_index_link_max = self.dv_index_max;
            
            % set equality for each junction using conservation
            % We did not preallocate memory for equality matrix
            constraints_eq = setEqConstraints(net, self.dv_index, self.dv_index_max);
            self.Aeq = constraints_eq.EqMatrix;
            self.beq = 0*self.Aeq(:,1);
            self.size_Aeq = size(self.Aeq);
            
            % if we have junction, then add auxiliary decision variable
            % which is the L1 norm e = |q1 - Rq2|
            for junc = net.junc_labels'

                % set up auxilary variables for e = q2-Rq1 at each merge
                % or diverge
                juncStr = sprintf('junc_%d',junc);
                num_steps = length(net.network_junc.(juncStr).T);
                
                if strcmp(net.network_junc.(juncStr).type_junc,'merge') ||...
                   strcmp(net.network_junc.(juncStr).type_junc,'diverge')
               
                    self.dv_index.(juncStr) = zeros(2,1);
                    self.dv_index.(juncStr)(1) = self.dv_index_max + 1;
                    self.dv_index.(juncStr)(2) = self.dv_index_max + ...
                        length(net.network_junc.(juncStr).T);
                    self.dv_index_max = self.dv_index_max + ...
                        length(net.network_junc.(juncStr).T);
                end
                
                % Add additional constraints e = |q1-Rq2|
                tmpMatrix = zeros(0, self.dv_index_max);
                    
                if strcmp(net.network_junc.(juncStr).type_junc,'merge')
                        
                    % parameter conditions
                    R_priority = net.network_junc.(juncStr).ratio(2)...
                        /net.network_junc.(juncStr).ratio(1);
                    
                    for step = 1:num_steps
                            
                        num_row = size(tmpMatrix,1);
                        linkStr = sprintf('link_%d',net.network_junc.(juncStr).inlabel(2));
                        tmpMatrix(num_row+1, self.dv_index.(linkStr)(1, INDEX_DOWN) + step - 1)= 1;
                        linkStr = sprintf('link_%d',net.network_junc.(juncStr).inlabel(1));
                        tmpMatrix(num_row+1, self.dv_index.(linkStr)(1, INDEX_DOWN) + step - 1)= -R_priority;
                        tmpMatrix(num_row+1, self.dv_index.(juncStr)(1) - 1 + step)= -1;
                        
                        linkStr = sprintf('link_%d',net.network_junc.(juncStr).inlabel(2));
                        tmpMatrix(num_row+1, self.dv_index.(linkStr)(1, INDEX_DOWN) + step - 1)= -1;
                        linkStr = sprintf('link_%d',net.network_junc.(juncStr).inlabel(1));
                        tmpMatrix(num_row+1, self.dv_index.(linkStr)(1, INDEX_DOWN) + step - 1)= R_priority;
                        tmpMatrix(num_row+1, self.dv_index.(juncStr)(1) - 1 + step)= -1;
                             
                    end
                    
                elseif strcmp(net.network_junc.(juncStr).type_junc,'diverge')
                    
                    % parameter conditions
                    R_priority = net.network_junc.(juncStr).ratio(2)...
                        /net.network_junc.(juncStr).ratio(1);
                    
                    for step = 1:num_steps
                            
                        num_row = size(tmpMatrix,1);
                        linkStr = sprintf('link_%d',net.network_junc.(juncStr).outlabel(2));
                        tmpMatrix(num_row+1, self.dv_index.(linkStr)(1, INDEX_UP) + step-1)= 1;
                        linkStr = sprintf('link_%d',net.network_junc.(juncStr).outlabel(1));
                        tmpMatrix(num_row+1, self.dv_index.(linkStr)(1, INDEX_UP) + step-1)= -R_priority;
                        tmpMatrix(num_row+1, self.dv_index.(juncStr)(1) - 1 + step)= -1;
                        
                        linkStr = sprintf('link_%d',net.network_junc.(juncStr).outlabel(2));
                        tmpMatrix(num_row+1, self.dv_index.(linkStr)(1, INDEX_UP) + step-1)= -1;
                        linkStr = sprintf('link_%d',net.network_junc.(juncStr).outlabel(1));
                        tmpMatrix(num_row+1, self.dv_index.(linkStr)(1, INDEX_UP) + step-1)= R_priority;
                        tmpMatrix(num_row+1, self.dv_index.(juncStr)(1) - 1 + step)= -1;

                    end
                        
                end
                    
                % Add to Aineq
                self.Aineq( self.size_Aineq(1) + 1: self.size_Aineq(1) + size(tmpMatrix,1),...
                           1 : size(tmpMatrix,2)) = tmpMatrix;
                self.bineq( self.size_Aineq(1) + 1: self.size_Aineq(1) + size(tmpMatrix,1) ) = 0;
                 
                self.size_Aineq(1) = self.size_Aineq(1) + size(tmpMatrix, 1);
                self.size_Aineq(2) = size(tmpMatrix,2);
                
                % augment equality matrix to same number of columns
                self.Aeq(:, self.size_Aeq(2)+1:self.size_Aineq(2)) = 0;
                self.size_Aeq = size(self.Aeq);
                
            end
            
            % We preallocated the memory for the matrix with all entries initialized as 0
            % When Aineq is constructed, truncate the last 0 rows
            self.Aineq = self.Aineq(1:self.size_Aineq(1),1:self.size_Aineq(2));
            self.bineq = self.bineq(1:self.size_Aineq(1),1);
            
            % to match cplex Ax <= b format
            self.Aineq = -self.Aineq;
            self.bineq = -self.bineq;
            
            % set upper and lower bounds for initial and boundary condition
            % variables.
            for link = net.link_labels'
                
                linkStr = sprintf('link_%d',link);
                
                % lower bound for q and rho
                self.lb(self.dv_index.(linkStr)(1, INDEX_UP):self.dv_index.(linkStr)(2, INDEX_UP)) = 0;
                self.lb(self.dv_index.(linkStr)(1, INDEX_DOWN):self.dv_index.(linkStr)(2, INDEX_DOWN)) = 0;
                self.lb(self.dv_index.(linkStr)(1, INDEX_INI):self.dv_index.(linkStr)(2, INDEX_INI)) = 0;
                
                % upper bound for q and rho
                self.ub(self.dv_index.(linkStr)(1, INDEX_UP):self.dv_index.(linkStr)(2, INDEX_UP)) = ...
                    1.0*net.network_hwy.(linkStr).para_qmax;
                self.ub(self.dv_index.(linkStr)(1, INDEX_DOWN):self.dv_index.(linkStr)(2, INDEX_DOWN)) = ...
                    1.0*net.network_hwy.(linkStr).para_qmax;
                self.ub(self.dv_index.(linkStr)(1, INDEX_INI):self.dv_index.(linkStr)(2, INDEX_INI)) = ...
                    1.0*net.network_hwy.(linkStr).para_km;
            end
            
            % set the upper and lower bound for auxiliary vairlabels e at
            % junctions
            for junc = net.junc_labels'
                
                if strcmp(net.network_junc.(juncStr).type_junc,'merge') ||...
                        strcmp(net.network_junc.(juncStr).type_junc,'diverge')
                    juncStr = sprintf('junc_%d',junc);
                    
                    self.lb( self.dv_index.(juncStr)(1):self.dv_index.(juncStr)(2)) = 0;
                    self.ub( self.dv_index.(juncStr)(1):self.dv_index.(juncStr)(2)) = 100;
                end
                
            end
            
           
            % set variable type, here assume all continuous
            self.ctype(1:self.dv_index_max) = 'C';
        end
        
        
        %===============================================================
        % solve using cplex
        % H and f are objective functions
        % min x'*H*x + f*x
        % s.t. Ax < b
        function [x, fval, exitflag, output] = solveProgram(self)
            
            try
                % CPLEX
                % min f
                % s.t. Aineq*x <= bineq
                %      Aeq*x = beq
                %                 options = cplexoptimset('Diagnostics', 'on', 'TolXInteger', 10e-10);
                %                 [x, fval, exitflag, output] = ...
                %                     cplexmiqp(obj.H, obj.f, obj.Aineq, obj.bineq, obj.Aeq, obj.beq,...
                %                     [ ], [ ], [ ], obj.lb, obj.ub, obj.ctype, [ ], options);
                %                 [x, fval, exitflag, output] = ...
                %                     cplexqp(obj.H, obj.f, obj.Aineq, obj.bineq, obj.Aeq, obj.beq,...
                %                      obj.lb, obj.ub);
                [x, fval, exitflag, output] = ...
                    cplexlp(self.f, self.Aineq, self.bineq, self.Aeq, self.beq,...
                    self.lb, self.ub);
                
                fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
                fprintf ('Solution value = %f \n', fval);
            catch m
                throw (m);
            end
            
            if (~isempty(x))
                %fprintf('\n Objective function:\n Min:\n');
            end
            
        end
        
        %===============================================================
        % Add entropy condition
        % Intuition: 
        % min sum  -w(i)T(i)*(alpha*(q1(i)+q2(i)) - beta*|q2(i)-R*q1(i)|))...
        % need w(i) > w(i+1) while alpha and beta satisfying conditions
        % input:
        %       network: the network object
        %       junction: nx1 vector, the junctions that we want to add
        %           entropic condition; NOTE: for now, the theory only
        %           support one junction
        function addEntropy(self, net, junction)
            global INDEX_UP INDEX_DOWN
            
            for junc = junction'
                                
                juncStr = sprintf('junc_%d',junc);
                
                % parameters for setting entropic condition at junctions
                num_steps = length(net.network_junc.(juncStr).T);
                                
                % initialize the objective function
                if isempty(self.f) 
                    % initialize with 0
                    self.f = zeros(self.dv_index_max,1);
                end
                
                % if connection
                if strcmp(net.network_junc.(juncStr).type_junc, 'connection')
                                        
                    linkStr = sprintf('link_%d', net.network_junc.(juncStr).inlabel); 
                    
                    % for a connection, simply put a linear weight
                    f_junc = zeros(self.dv_index_max, 1);
                    f_junc(self.dv_index.(linkStr)(1, INDEX_DOWN):...
                            self.dv_index.(linkStr)(2, INDEX_DOWN),1) =...
                            -(num_steps:-1:1)'.*net.network_junc.(juncStr).T;
                        
                    self.f = self.f + f_junc;
                    
                % if merge/diverge For each Merge or diverge junction, 
                % we add varaibel e = |q2-Rq1| for each step
                elseif strcmp(net.network_junc.(juncStr).type_junc, 'merge') ||...
                       strcmp(net.network_junc.(juncStr).type_junc, 'diverge') 
                    
                    % H is the matrix for defining quadratic functions
                    if isempty(self.H)
                        self.H = zeros(self.dv_index_max);
                    end
                    
                    if isempty(self.f)
                        self.f = zeros(self.dv_index_max,1);
                    end
                    
                    % set entropy for this junction
                    f_junc = zeros(self.dv_index_max ,1);
                    
                    % parameter conditions
                    R_priority = net.network_junc.(juncStr).ratio(2)...
                        /net.network_junc.(juncStr).ratio(1);
                    
                    beta_pri = 1;
                    
                    % compute the exponential factor
                    T_junc = net.network_junc.(juncStr).T;
                    if R_priority <= 1
                        E = ( (R_priority+2)*max(T_junc)/min(T_junc) + 1 )...
                                /( R_priority + 1 ) + 0.1;
                    else
                        E = ( R_priority + (1+2*R_priority)*max(T_junc)/min(T_junc) )...
                            /( R_priority + 1 ) + 0.1;
                    end
                    
                    % compute the alpha parameter
                    if R_priority <=1
                        alpha_pri = 0.5*( (E+1)*beta_pri/(E-1) +...
                            (E+E*R_priority-1)*beta_pri );
                    else
                        alpha_pri = 0.5*( R_priority*(E+1)*beta_pri/(E-1) + ...
                            (E+E*R_priority-R_priority)*beta_pri );
                    end
                    
                    % generate the weight vector
                    weight = 0.01*ones(num_steps, 1);
                    for j = 1:num_steps-1
                        weight(j+1) = E*weight(j);
                    end
                    weight(1:num_steps) = weight(num_steps:-1:1);
                    
                    % add entropic objective component
                    % sum -w(i)T(i)alpha{q1(i)+q2(i)}
                    if strcmp(net.network_junc.(juncStr).type_junc,'merge')
                        % entropy condition:
                        for link = net.network_junc.(juncStr).inlabel'
                            linkStr = sprintf('link_%d', link);
                            f_junc(self.dv_index.(linkStr)(1, INDEX_DOWN):...
                                self.dv_index.(linkStr)(2, INDEX_DOWN),1) =...
                                weight.*T_junc;
                        end
                    elseif strcmp(net.network_junc.(juncStr).type_junc,'diverge')
                        % entropy condition:
                        for link = net.network_junc.(juncStr).outlabel'
                            f_junc(self.dv_index.(linkStr)(1, INDEX_UP):...
                                self.dv_index.(linkStr)(2, INDEX_UP),1) =...
                                weight.*T_junc;
                        end
                    end
                    
                    self.f = self.f - alpha_pri*f_junc;
                    
                    % add the component: sum w(i) x T(i) x beta x e(i)
                    % Use e = |q1-Rq2|
                    self.f = self.f + beta_pri*[ zeros(self.dv_index.(juncStr)(1)-1,1);...
                        weight.*T_junc;...
                        zeros( self.dv_index_max - self.dv_index.(juncStr)(2) ,1)];
                    
                end     
            
            end
        end
        
        
        %===============================================================
        % maximum upstream flow for one link
        function maxUpflow(self, links)
            global INDEX_UP
            
            if isempty(self.f)
                % if not defiend by any function yet
                self.f = zeros(self.dv_index_max,1);
            end
            
            f_upflow = zeros(self.dv_index_max,1);
            
            for link = links'
                
                linkStr = sprintf('link_%d', link);
                
                T_grid = net.network_hwy.(linkStr).T_us;
            
                num_steps = length(T_grid);
                f_upflow(self.dv_index.(linkStr)(1,INDEX_UP):...
                         self.dv_index.(linkStr)(2,INDEX_UP)) = -(num_steps:-1:1)'.*T_grid;
            end
            
            self.f = self.f + f_upflow;

        end
        
        
        %===============================================================
        % maximum downstream flow for one link
        function maxDownflow(self, net, links)
            global INDEX_DOWN
            
            if isempty(self.f)
                % if not defiend by any function yet
                self.f = zeros(self.dv_index_max,1);
            end
            
            f_downflow = zeros(self.dv_index_max,1);
            
            for link = links'
                
                linkStr = sprintf('link_%d', link);
                
                T_grid = net.network_hwy.(linkStr).T_ds;
            
                num_steps = length(T_grid);
                f_downflow(self.dv_index.(linkStr)(1,INDEX_DOWN):...
                         self.dv_index.(linkStr)(2,INDEX_DOWN)) = -(num_steps:-1:1)'.*T_grid;
            end
            
            self.f = self.f + f_downflow;
 
        end
        
        

        
    end
    
end













