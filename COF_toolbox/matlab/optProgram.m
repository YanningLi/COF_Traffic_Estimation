% Yanning Li, Sep 01, 2015
% This class construct a the optimization program. It sets the linear
% constraints of the optimization program, has multiple functions to set
% the objective function, and solve the program using cplex.

% TODO:
% The entropy condition works for only one junction. 



classdef optProgram < handle
    
    properties
        % configuration
        start_time;        
        end_time;
        
        % Decision variable locator, s atruct
        % Dec.(linkStr) = [start of upstream, start of downstream, start of initial;
        %                  end of upstream, end of downstream, end of initial]
        Dec;
        Dec_max;    % the maximal index
        
        % auxiliary variables ( the L1 error e=|q1-Rq1|) at merges or diverges 
        % that we would like to add in the decision variable.
        % L1_error_index.(juncStr) = [start index; end index]
        L1_error_index;
        L1_error_index_max;
        
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
        function obj = optProgram(~)
            obj.start_time = 0;
            
            obj.end_time = 0;
            
            obj.Dec = struct;
            obj.Dec_max = 0;
            obj.L1_error_index = struct;
            obj.L1_error_index_max = 0;
            
            
            % preallocate memory to speed up
            obj.Aineq = zeros(100000,2000);
            obj.bineq = zeros(100000,1);
            obj.size_Aineq = [0 0];
            
            obj.Aeq = zeros(0,0);
            obj.beq = zeros(0,1);
            obj.size_Aeq = [0 0];
            
            obj.lb = zeros(0,1);
            obj.ub = zeros(0,1);
            obj.ctype = '';
            
            obj.f = [];
            obj.H = [];
            
        end
        
        
        %===============================================================
        % set configuration for simulation
        function setConfig(obj, start_time, end_time)
            
            obj.start_time = start_time;
            obj.end_time = end_time;
            
        end
        
        
        %===============================================================
        % set Linear constraints
        function setConstraints(obj, network, e_max)
            global LINKLABEL TIME V_PARA W_PARA K_C_PARA K_M_PARA Q_MAX_PARA
            global POSTKM P_MIN P_MAX DENS 
            global INDEX_UP INDEX_DOWN INDEX_INI
            
            % set inequality for each link
            for link = network.link_labels'
                
                linkStr = sprintf('link_%d',link);
                
                % boundary conditions
                q_us = network.q_upstream(network.q_upstream(:,LINKLABEL)==link,:);
                q_ds = network.q_downstream(network.q_downstream(:,LINKLABEL)==link,:);
                
                % initial conditions
                Con_ini = network.Con_initial(network.Con_initial(:,TIME)==obj.start_time &...
                    network.Con_initial(:,LINKLABEL)==link, :);
                init_seg = unique([Con_ini(:,P_MIN),Con_ini(:,P_MAX)]);
                num_ini_seg = size(Con_ini,1);
                
                % set inequality constraints
                constraints_ineq = setIneqConstraints(...
                    network.network_hwy.(linkStr).paras(V_PARA),...
                    network.network_hwy.(linkStr).paras(W_PARA),...
                    network.network_hwy.(linkStr).paras(K_C_PARA),...
                    network.network_hwy.(linkStr).paras(K_M_PARA),...
                    network.network_hwy.(linkStr).profile(1,POSTKM)*1000,...
                    obj.start_time, obj.end_time,...
                    q_us, q_ds,...
                    init_seg, Con_ini(:,DENS), e_max);
                
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
                obj.Aineq(obj.size_Aineq(1)+1:obj.size_Aineq(1) + num_row_constraints,...
                    obj.size_Aineq(2)+1:obj.size_Aineq(2)+ num_col_constraints) ...
                    = link_constraints(:,1:num_col_constraints);
                
                obj.bineq(obj.size_Aineq(1)+1:obj.size_Aineq(1)+size_constraints(1),1)...
                    = link_constraints(:,num_col_constraints+1);
                
                % update Aineq matrix size
                obj.size_Aineq = obj.size_Aineq + [num_row_constraints, num_col_constraints];
                
                % update index for decision variable
                obj.Dec.(linkStr) = zeros(2,3);
                % add index for upstream boundary flows
                obj.Dec.(linkStr)(1, INDEX_UP) = obj.Dec_max + 1;
                obj.Dec.(linkStr)(2, INDEX_UP) = obj.Dec_max + size(q_us, 1);
                obj.Dec_max = obj.Dec_max + size(q_us, 1);
                % add index for downstream boundary flows
                obj.Dec.(linkStr)(1, INDEX_DOWN) = obj.Dec_max + 1;
                obj.Dec.(linkStr)(2, INDEX_DOWN) = obj.Dec_max + size(q_ds, 1);
                obj.Dec_max = obj.Dec_max + size(q_ds, 1);
                % add index for initial densities
                obj.Dec.(linkStr)(1, INDEX_INI) = obj.Dec_max + 1;
                obj.Dec.(linkStr)(2, INDEX_INI) = obj.Dec_max + num_ini_seg;
                obj.Dec_max = obj.Dec_max + num_ini_seg;
                    
            end
            
            
            % set equality for each junction using conservation
            constraints_eq = setEqConstraints(network,obj.Dec);
            obj.Aeq = constraints_eq.EqMatrix;
            obj.beq = 0*obj.Aeq(:,1);
            
            % if we have junction, then add auxiliary decision variable
            % which is the L1 norm e = |q1 - Rq2|
            for junc = network.junc_labels'

                % number of auxilary variables for e = q2-Rq1 at each merge
                % or diverge
                juncStr = sprintf('junc_%d',junc);
                if strcmp(network.network_junc.(juncStr).type_junc,'merge') ||...
                   strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                    obj.L1_error_index.(juncStr) = zeros(2,1);
                    obj.L1_error_index.(juncStr)(1,1) = obj.L1_error_index_max + 1;
                    obj.L1_error_index.(juncStr)(1,1) = obj.L1_error_index_max + ...
                        length(network.network_junc.(juncStr).T);
                    obj.L1_error_index_max = obj.L1_error_index_max + ...
                        length(network.network_junc.(juncStr).T);
                end
                
                
                % Add additional constraints e = |q1-Rq2|
                tmpMatrix = zeros(0, obj.L1_error_index_max);
                    
                if strcmp(network.network_junc.(juncStr).type_junc,'merge')
                        
                    links = [network.network_junc.(juncStr).inlabel...
                            network.network_junc.(juncStr).outlabel];
                    
                    % parameter conditions
                    R_priority = network.network_junc.(juncStr).ratio(2)...
                        /network.network_junc.(juncStr).ratio(1);
                    
                    for step = 1:num_steps
                            
                        num_row = size(tmpMatrix,1);
                        linkStr = sprintf('link_%d',links(2));
                        tmpMatrix(num_row+1, obj.Dec.(linkStr)(1, INDEX_DOWN) + step-1)= 1;
                        linkStr = sprintf('link_%d',links(1));
                        tmpMatrix(num_row+1, obj.Dec.(linkStr)(1, INDEX_DOWN) + step-1)= -R_priority;
                        tmpMatrix(num_row+1, obj.L1_error_index.(juncStr)(1,1) - 1 + step)= -1;
                        
                        linkStr = sprintf('link_%d',links(2));
                        tmpMatrix(num_row+1, obj.Dec.(linkStr)(1, INDEX_DOWN) + step-1)= -1;
                        linkStr = sprintf('link_%d',links(1));
                        tmpMatrix(num_row+1, obj.Dec.(linkStr)(1, INDEX_DOWN) + step-1)= R_priority;
                        tmpMatrix(num_row+1, obj.L1_error_index.(juncStr)(1,1) - 1 + step)= -1;
                             
                    end
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                    
                    links = [network.network_junc.(juncStr).outlabel...
                        network.network_junc.(juncStr).inlabel];
                    
                    % parameter conditions
                    R_priority = network.network_junc.(juncStr).ratio(2)...
                        /network.network_junc.(juncStr).ratio(1);
                    
                    for step = 1:num_steps
                            
                        num_row = size(tmpMatrix,1);
                        linkStr = sprintf('link_%d',links(2));
                        tmpMatrix(num_row+1, obj.Dec.(linkStr)(1, INDEX_UP) + step-1)= 1;
                        linkStr = sprintf('link_%d',links(1));
                        tmpMatrix(num_row+1, obj.Dec.(linkStr)(1, INDEX_UP) + step-1)= -R_priority;
                        tmpMatrix(num_row+1, obj.L1_error_index.(juncStr)(1,1) - 1 + step)= -1;
                        
                        linkStr = sprintf('link_%d',links(2));
                        tmpMatrix(num_row+1, obj.Dec.(linkStr)(1, INDEX_UP) + step-1)= -1;
                        linkStr = sprintf('link_%d',links(1));
                        tmpMatrix(num_row+1, obj.Dec.(linkStr)(1, INDEX_UP) + step-1)= R_priority;
                        tmpMatrix(num_row+1, obj.L1_error_index.(juncStr)(1,1) - 1 + step)= -1;

                    end
                        
                end
                    
                % Add to Aineq
                obj.Aineq( obj.size_Aineq(1) + 1: obj.size_Aineq(1) + size(tmpMatrix,1),...
                           1 : size(tmpMatrix,2)) = tmpMatrix;
                obj.size_Aineq = size(obj.Aineq);
                
                obj.bineq = zeros(size(tmpMatrix,1),1);
                
                % augment equality matrix to same size
                obj.Aeq(:, obj.size_Aeq(2)+1:obj.size_Aineq(2)) = 0;
                obj.size_Aeq = size(obj.Aeq);
                
            end
            
            % if no junctions, still update the number of decision
            % variables
            obj.L1_error_index_max = obj.Dec_max;
            
            
            % We preallocated the memory for the matrix with all entries initialized as 0
            % When Aineq is constructed, truncate the last 0 rows
            obj.Aineq = obj.Aineq(1:obj.size_Aineq(1),1:obj.size_Aineq(2));
            obj.bineq = obj.bineq(1:obj.size_Aineq(1),1);
            
            % to match cplex Ax <= b format
            obj.Aineq = -obj.Aineq;
            obj.bineq = -obj.bineq;
            
            
            % set upper and lower bounds for initial and boundary condition
            % variables.
            for link = network.link_labels'
                
                linkStr = sprintf('link_%d',link);
                
                % lower bound
                obj.lb(obj.Dec.(linkStr)(1, INDEX_UP):obj.Dec.(linkStr)(2, INDEX_UP)) = 0;
                obj.lb(obj.Dec.(linkStr)(1, INDEX_DOWN):obj.Dec.(linkStr)(2, INDEX_DOWN)) = 0;
                obj.lb(obj.Dec.(linkStr)(1, INDEX_INI):obj.Dec.(linkStr)(2, INDEX_INI)) = 0;
                
                % upper bound for q and rho
                obj.lb(obj.Dec.(linkStr)(1, INDEX_UP):obj.Dec.(linkStr)(2, INDEX_UP)) = ...
                    1.0*network.network_hwy.(linkStr).paras(Q_MAX_PARA);
                obj.lb(obj.Dec.(linkStr)(1, INDEX_DOWN):obj.Dec.(linkStr)(2, INDEX_DOWN)) = ...
                    1.0*network.network_hwy.(linkStr).paras(Q_MAX_PARA);
                obj.lb(obj.Dec.(linkStr)(1, INDEX_INI):obj.Dec.(linkStr)(2, INDEX_INI)) = ...
                    1.0*network.network_hwy.(linkStr).paras(K_M_PARA);
            end
            
            % lower and upper bound and ctype for e
            obj.lb( obj.Dec_max+1: obj.L1_error_index_max) = 0; 
            obj.ub( obj.Dec_max+1: obj.L1_error_index_max) = 100; 
                    
            % set variable type, here assume all continuous
            obj.ctype(1:obj.Dec_max) = 'C';
            obj.ctype( obj.Dec_max+1: obj.L1_error_index_max) = 'C';     
        end
        
        
        %===============================================================
        % solve using cplex
        % H and f are objective functions
        % min x'*H*x + f*x
        % s.t. Ax < b
        function [x, fval, exitflag, output] = solveProgram(obj)
            
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
                    cplexlp(obj.f, obj.Aineq, obj.bineq, obj.Aeq, obj.beq,...
                    obj.lb, obj.ub);
                
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
        %       E: float, the exponential weight derived from conditions, read
        %           paper for details
        %       
        function addEntropy(obj, network, junction, E)
            global INDEX_UP INDEX_DOWN
            
            for junc = junction'
                                
                juncStr = sprintf('junc_%d',junc);
                
                % parameters for setting entropic condition at junctions
                num_steps = length(network.network_junc.(juncStr).T);
                                
                % generate the weight
                weight = 0.01*ones(num_steps, 1);
                for j = 1:num_steps-1
                    weight(j+1) = E*weight(j);
                end
                weight(1:num_steps) = weight(num_steps:-1:1);
                
                % update the objective function
                if length(obj.f) ~= obj.L1_error_index_max
                    % pad with 0
                    obj.f( length(obj.f)+1: obj.L1_error_index_max ) = 0;
                end
                
                % if connection
                if strcmp(network.network_junc.(juncStr).type_junc, 'connection')
                    
                    inLink = network.network_junc.(juncStr).inlabel;
                    
                    linkStr = sprintf('link_%d', inLink); 
                    % for a connection, simply put a linear weight
                    f_junc = zeros( obj.L1_error_index_max, 1);
                    f_junc(obj.Dec.(linkStr)(1, INDEX_DOWN):...
                            obj.Dec.(linkStr)(2, INDEX_DOWN),1) =...
                            -(num_steps:-1:1)'.*network.network_junc.(juncStr).T;
                        
                    obj.f = obj.f + f_junc;
                    
                % if merge/diverge For each Merge or diverge junction, 
                % we add varaibel e = |q2-Rq1| for each step
                elseif strcmp(network.network_junc.(juncStr).type_junc, 'merge') ||...
                       strcmp(network.network_junc.(juncStr).type_junc, 'diverge') 
                    
                    
                    if isempty(obj.H)
                        obj.H = zeros(obj.L1_error_index_max);
                    end
                    
                    % set entropy for this junction
                    f_junc = zeros(obj.L1_error_index_max ,1);
                    
                    if strcmp(network.network_junc.(juncStr).type_junc, 'merge')
                        links = [network.network_junc.(juncStr).inlabel...
                            network.network_junc.(juncStr).outlabel];
                    elseif strcmp(network.network_junc.(juncStr).type_junc, 'diverge')
                        links = [network.network_junc.(juncStr).outlabel...
                        network.network_junc.(juncStr).inlabel];
                    end
                    
                    % parameter conditions
                    R_priority = network.network_junc.(juncStr).ratio(2)...
                        /network.network_junc.(juncStr).ratio(1);
                    
                    beta_pri = 1;
                    
                    % compute the alpha parameter
                    if R_priority <=1
                        alpha_pri = 0.5*( (E+1)*beta_pri/(E-1) +...
                            (E+E*R_priority-1)*beta_pri );
                    else
                        alpha_pri = 0.5*( R_priority*(E+1)*beta_pri/(E-1) + ...
                            (E+E*R_priority-R_priority)*beta_pri );
                    end
      
                    % add entropic objective component
                    if strcmp(network.network_junc.(juncStr).type_junc,'merge')
                        % entropy condition:
                        for link = links
                            linkStr = sprintf('link_%d', link);
                            f_junc(obj.Dec.(linkStr)(1, INDEX_DOWN):...
                                obj.Dec.(linkStr)(2, INDEX_DOWN),1) =...
                                weight.*network.network_junc.(juncStr).T;
                        end
                    elseif strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                        % entropy condition:
                        for link = links
                            f_junc(obj.Dec.(linkStr)(1, INDEX_UP):...
                                obj.Dec.(linkStr)(2, INDEX_UP),1) =...
                                weight.*network.network_junc.(juncStr).T;
                        end
                   
                    
                    obj.f = obj.f - alpha_pri*f_junc;
                    % Use e = |q1-Rq2|
                    obj.f = obj.f + beta_pri*[ zeros(obj.L1_error_index.(juncStr)(1)-1,1);...
                        weight.*network.network_junc.(juncStr).T;...
                        zeros( obj.L1_error_index_max - obj.L1_error_index.(juncStr)(2) ,1)];
                    end
                    
                end     %end for each junction
            
            end
        end
        
        
        %===============================================================
        % maximum upstream flow for one link
        function maxUpflow(obj,links, T_in)
            global NUM_UP
            
            num_steps = length(T_in);
            leng_linkVar = obj.Dec(size(obj.Dec,1),size(obj.Dec,2),2);
            
            if isempty(obj.f)
                % if not defiend by any function yet
                obj.f = zeros(leng_linkVar,1);
                f_upflow = zeros(leng_linkVar,1);
            elseif length(obj.f) > leng_linkVar
                % defined by entropy function; save the total length
                tt_len = length(obj.f);
                f_upflow = zeros(tt_len,1);
            elseif length(obj.f) == leng_linkVar
                % defined by maxDown or connection
                f_upflow = zeros(leng_linkVar,1);
            end
            
            
            for i=1:length(links)
                f_upflow(obj.Dec(links(i),NUM_UP,1):...
                    obj.Dec(links(i),NUM_UP,2)) = -(num_steps:-1:1)'.*T_in;
            end
            
            obj.f = obj.f + f_upflow;
            
        end
        
        
        %===============================================================
        % maximum downstream flow for one link
        function maxDownflow(obj, links, T_out)
            global NUM_DOWN
            
            num_steps = length(T_out);
            leng_linkVar = obj.Dec(size(obj.Dec,1),size(obj.Dec,2),2);
            
            if isempty(obj.f)
                % if not defiend by any function yet
                obj.f = zeros(leng_linkVar,1);
                f_downflow = zeros(leng_linkVar,1);
            elseif length(obj.f) > leng_linkVar
                % defined by entropy function; save the total length
                tt_len = length(obj.f);
                f_downflow = zeros(tt_len,1);
            elseif length(obj.f) == leng_linkVar
                f_downflow = zeros(leng_linkVar,1);
            end
            
            for i=1:length(links)
                f_downflow(obj.Dec(links(i),NUM_DOWN,1):...
                    obj.Dec(links(i),NUM_DOWN,2)) = -(num_steps:-1:1)'.*T_out;
            end
            
            obj.f = obj.f + f_downflow;
            
        end
        
        
        %===============================================================
        % modified for new discretization class _grid
        
        % onramp control at merges (only merges);
        % Assume only the onramp flow is controllable and the two mainline
        % flows hit their maximum using entropy condition;
        function addOnRampControl(obj,network,junction)
            global LINKTYPE NUM_DOWN ONRAMP MAINLINE
            
            for i = 1:length(junction)
                
                junc = junction(i); % the array of junction to add entropy
                juncStr = sprintf('junc_%d',junc);
               
                % parameters for merges and diverges
                num_steps = length(network.network_junc.(juncStr).T);
                
                leng_linkVar = obj.Dec(size(obj.Dec,1),size(obj.Dec,2),2);
                
                % set the length of the objective function vector
                if isempty(obj.f)
                    % if not defiend by any function yet
                    obj.f = zeros(leng_linkVar,1);
                    f_onRampCtrl = zeros(leng_linkVar,1);
                elseif length(obj.f) > leng_linkVar
                    % defined by entropy function; save the total length
                    tt_len = length(obj.f);
                    f_onRampCtrl = zeros(tt_len,1);
                elseif length(obj.f) == leng_linkVar
                    f_onRampCtrl = zeros(leng_linkVar,1);
                end
                
                links = network.network_junc.(sprintf('junc_%d', junc)).inlabel;
                
                tmpLinkStr1 =  sprintf('link_%d',links(1));
                tmpLinkStr2 =  sprintf('link_%d',links(2));
                if (network.network_hwy.(tmpLinkStr1).profile(LINKTYPE) == ONRAMP &&...
                        network.network_hwy.(tmpLinkStr2).profile(LINKTYPE) == MAINLINE)
                    onRampLink = links(1);
                    inLink = links(2);
                elseif  (network.network_hwy.(tmpLinkStr1).profile(LINKTYPE) == MAINLINE &&...
                        network.network_hwy.(tmpLinkStr2).profile(LINKTYPE) == ONRAMP)
                    onRampLink = links(2);
                    inLink = links(1);
                else
                    error('Junc_%d has no onramp',junc);
                end
                
                % entropy condition:
                % Min F = -sum(weight*qmain.*T*2) - sum(weight*qon*T);
                f_onRampCtrl(obj.Dec(inLink,NUM_DOWN,1):...
                    obj.Dec(inLink,NUM_DOWN,2)) = ...
                    -(num_steps:-1:1)'.*network.network_junc.(juncStr).T;
                
                f_onRampCtrl(obj.Dec(onRampLink,NUM_DOWN,1):...
                    obj.Dec(onRampLink,NUM_DOWN,2)) = ...
                    - network.network_junc.(juncStr).T*2*num_steps;
                    %-(num_steps:-1:1)'.*network.network_junc.(juncStr).T*2;
                
                obj.f = obj.f + f_onRampCtrl;
                
                
            end
            
            
        end
        
        
    end
    
end













