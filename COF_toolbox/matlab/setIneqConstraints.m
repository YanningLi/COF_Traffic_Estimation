% Yanning Li, Sep 01, 2015

% This class sets the inequality constraints for the optimization program.
% It is based on the setIneqConstraints_grid.m which speeds up the program by
% only adding the minimum set of constarints. See the code description for
% setIneqConstraints_grid.m for details.

% Sep 01, 2015
% Remark: this version does not support internal and density condition with
% uneven discretization of the boundary flows.


classdef setIneqConstraints
    
    properties
        
        n_max_us;  % num of upstream steps 
        n_max_ds;  % num of downstream steps
        m_max;  % num of internal conditions
        b_max;  % num of initial conditions
        u_max;  % num of density conditions
        
        v;
        w;
        k_c;
        k_m;
        
        start_time;
        end_time;   % end of simulation time window
        
        xi;     % post meter of the upstream position
        chi;    % post meter of the downstream position
        
        % up and downstream time grid
        T_us;       %upstream_steps x 1
        T_us_cum;   % cumulative; upstream_steps+1 x 1
        T_ds;
        T_ds_cum;
        
        e_max;  % error for historical data, assume to be larger +-20%
        list;
        
        ini_segments;   % uneven space descritization of link b_max+1 x 1
        rho_ini;    % b_max x 1
        indiced_rho_ini;    
        X;
        
        % for boundary conditions
        qin_meas;
        qout_meas;
        
        % for internal conditions
        x_min;
        x_max;
        t_min;
        t_max;
        v_meas;
        r_meas;
        
        % for density conditions
        x_min_u;
        x_max_u;
        t_u;
        dens_meas;
        
        % auxiliary binary variables
        nb_min;
        nb_max;
        nb_min_u;
        nb_max_u;
        
        % auxiliary variables
        num_aux;
        size_row;
    end
    
    methods
        
        function si = setIneqConstraints(...
                v, w, k_c, k_m, postm,...
                start_time, end_time,...
                qin_meas, qout_meas,...
                ini_segments, rho_ini, e_max)
            global T_MIN T_MAX
            
            
            si.n_max_us = size(qin_meas,1);   % boundary conditions
            si.n_max_ds = size(qout_meas,1);
            
            si.b_max = size(rho_ini,1);   % initial conditions
            si.m_max = 0;                 % disable internal condition
            si.u_max = 0;                 % disable density condition
            
            si.v = v;
            si.w = w;
            si.k_c = k_c;
            si.k_m = k_m;
            
            si.start_time = start_time;
            si.end_time = end_time;
            
            si.xi = 0;  % relative postm for each link
            si.chi = postm;
            
            % if still using uniform discretization
            si.T_us = qin_meas(:,T_MAX)-qin_meas(:,T_MIN);
            si.T_us_cum = [0; qin_meas(:,T_MAX)];
            si.T_ds = qout_meas(:,T_MAX)-qout_meas(:,T_MIN);
            si.T_ds_cum = [0; qout_meas(:,T_MAX)];
            
            
            % set boundary condition
            si.qin_meas = qin_meas;
            si.qout_meas = qout_meas;
            
            % disable internal condition
            si.x_min = [];
            si.x_max = [];
            si.t_min = [];
            si.t_max = [];
            si.v_meas = [];
            si.r_meas = [];
            
            %disable density condition
            si.x_min_u = [];
            si.x_max_u = [];
            si.t_u = [];
            si.dens_meas = [];
            
            si.e_max = e_max;   % assume 0 error from the data
            
            % set initial condition
            if ~iscolumn(ini_segments)
                ini_segments = ini_segments';
            end
            si.ini_segments = ini_segments;
            
            if (length(ini_segments)~=length(rho_ini)+1 && ~isempty(rho_ini))
                error('Check ini_segments and rho_ini dimensions');
            end
            
            si.X = si.ini_segments(2:length(si.ini_segments)) - si.ini_segments(1:length(si.ini_segments)-1);
            
            if (isempty(rho_ini) || ~all(~isnan(rho_ini)))
                sprintf('Warning: Not defining a complete set of initial value conditions could significantly increase the computational time!\n');
            end
            
            if ~iscolumn(rho_ini)
                rho_ini = rho_ini';
            end
            si.rho_ini = rho_ini;
            
            if (~isempty(rho_ini))
                block_ini = double(rho_ini > k_c); 
                block_ini(isnan(rho_ini)) = NaN;
                si.indiced_rho_ini = groupSameElement(block_ini,inf);
            else
                si.indiced_rho_ini = [];
            end
            
            si.num_aux = 0;
            
            si.size_row = si.n_max_us + si.n_max_ds +...
                          si.b_max + 2*si.m_max + 2*si.u_max + 1; 
            %This row size if defined as a temporary one to compute the number of needed binary variables
            [si.nb_min,si.nb_max, si.nb_min_u, si.nb_max_u] = si.getBinaryvar;
            
            si.size_row = si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + 2*si.u_max +...
                sum(si.nb_min) + sum(si.nb_max) + sum(si.nb_min_u) + sum(si.nb_max_u) + sum(si.num_aux) + 1;  
            si.list = [];
            
        end
        
        %==========================================================================
        %Method to substract to arrays with a "null" condition
        
        function [output] = substractArray(~,value1,value2)
            if ( (~isempty(value1))&& (~isempty(value2)) )
                output = value1-value2;
            else
                output = []; %return "null"
            end
        end
        
        %==========================================================================
        %Equations of function Mgamma, explicit solution of upstream condition
        
        function[array] = mgamma(si,n,t,x)
            
            array = zeros(1,si.size_row);  %Initialize the array
            
            if((si.T_us_cum(n) + (x - si.xi)/si.v <= t) &&...
                    (si.T_us_cum(n+1) + (x-si.xi)/si.v >= t ) && (n<= si.n_max_us))
                % in characteristic domain
                array(1,1:n-1) = si.T_us(1:n-1);
                array(1,n) = t - (x-si.xi) / si.v - si.T_us_cum(n);
                return
                
            elseif((si.T_us_cum(n+1) + (x-si.xi)/si.v < t) && (n <= si.n_max_us))
                % in Fan domain
                array(1,1:n) = si.T_us(1:n);
                array(1,si.size_row) = -si.k_c*si.v * (t - si.T_us_cum(n+1) - (x-si.xi)/si.v);
                return
                
            else
                array = [];     %return a "null"
                return
            end
            
        end
        
        %==========================================================================
        %Equations of function Mbeta, explicit solution of downstream condition
        
        function [array] = mbeta(si,n,t,x)
            array = zeros(1,si.size_row);
            
            if( (si.T_ds_cum(n) + (x - si.chi)/si.w <= t) &&...
                    (si.T_ds_cum(n+1) + (x-si.chi)/si.w >= t )&& (n <= si.n_max_ds))
                % in characteristic domain
                array(1,si.n_max_us + si.n_max_ds +1:...
                    si.n_max_us + si.n_max_ds + si.b_max) = -si.X; %initial number of vehicles
                array(1,si.n_max_us +1:si.n_max_us + n-1) = si.T_ds(1:n-1);
                array(1,si.n_max_us + n) = t - (x-si.chi)/si.w - si.T_ds_cum(n);
                array(1,si.size_row) = si.k_m*(x-si.chi);
                return
                
            elseif( (si.T_ds_cum(n+1) + (x-si.chi)/si.w < t) && (n <= si.n_max_ds))
                % in Fan domain
                array(1,si.n_max_us + si.n_max_ds +1:...
                    si.n_max_us + si.n_max_ds + si.b_max) = -si.X; %initial number of vehicles
                array(1,si.n_max_us +1:si.n_max_us + n) = si.T_ds(1:n);
                array(1,si.size_row) = -si.k_c*si.v * (t - si.T_ds_cum(n+1) - (x-si.chi)/(si.v));
                return
            else
                array = []; %return a "null"
                return
            end
        end
        
        %==========================================================================
        %Equations of function Mmu, explicit solutions of internal conditions
        
        function [array] = mmu(si, m, t, x)
            array = zeros(1,si.size_row);
            
            if( (x >= si.x_min(m) + si.v_meas(m)*(t-si.t_min(m))) &&...
                ( x >= si.x_max(m) + si.v*(t-si.t_max(m))) && (x <= si.x_min(m) + si.v*(t-si.t_min(m))) )
                array(1,si.n_max_us + si.n_max_ds + si.b_max + m) = 1;
                array(1,si.n_max_us + si.n_max_ds + si.b_max + si.m_max + m) =...
                    t - (x - si.x_min(m) - si.v_meas(m)*(t - si.t_min(m))) / (si.v - si.v_meas(m)) - si.t_min(m) ;
                return
            end
            
            if ( (x<= si.x_min(m) + si.v_meas(m)*(t-si.t_min(m))) &&...
                 ( x <= si.x_max(m) + si.w*(t-si.t_max(m))) && (x >= si.x_min(m) + si.w*(t-si.t_min(m))) )
                array(1,si.n_max_us + si.n_max_ds + si.b_max + m) = 1;
                array(1,si.n_max_us + si.n_max_ds + si.b_max + si.m_max + m) =...
                    t - (x - si.x_min(m) - si.v_meas(m)*(t-si.t_min(m)))/(si.w - si.v_meas(m)) - si.t_min(m);
                array(1,si.size_row) = -si.k_c*(si.v-si.w)*(x - si.x_min(m) - si.v_meas(m)*(t-si.t_min(m)))/(si.w - si.v_meas(m)) ;
                return
            end
            
            if( (x < si.x_max(m) + si.v*(t - si.t_max(m))) && ( x > si.x_max(m) + si.w*(t - si.t_max(m))) )
                array(1,si.n_max_us + si.n_max_ds + si.b_max + m) = 1;
                array(1, si.n_max_us + si.n_max_ds + si.b_max + si.m_max + m) = si.t_max(m) - si.t_min(m);
                array(1, si.size_row) = -(t-si.t_max(m))*si.k_c*(si.v - (x-si.x_max(m))/(t-si.t_max(m)+0.000001));
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        %==========================================================================
        %Equations of function MTau1, explicit solutions of initial solutions applying when p<pc
        
        function [array] = mtau1 (si,b,t,x)
            array = zeros(1,si.size_row);
            
            if ( (si.xi + sum(si.X(1:b-1)) + t*si.v <= x) && (x <=si.xi+ sum(si.X(1:b)) + t*si.v) && ( b <= si.b_max))
                array(1,si.n_max_us + si.n_max_ds + 1 : si.n_max_us + si.n_max_ds + b-1) = -si.X(1:b-1);
                
                array(1,si.n_max_us + si.n_max_ds + b) = (t*si.v - x + sum(si.X(1:b-1)) + si.xi);
                return
            end
            
            if ( (si.xi+ sum(si.X(1:b-1)) + t*si.w <= x) && (x < si.xi+ sum(si.X(1:b-1)) + t*si.v) && (b <= si.b_max))
                array(1,si.n_max_us + si.n_max_ds + 1: si.n_max_us + si.n_max_ds + b-1) = -si.X(1:b-1);
                
                array(1,si.size_row) = -si.k_c*(t*si.v - x + sum(si.X(1:b-1)) + si.xi);
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        %==========================================================================
        %Equations of function MTau2, explicit solutions of initial solutions applying when p>pc
        
        function [array] = mtau2 (si,b,t,x)
            array = zeros(1,si.size_row);
            if ( (si.xi + sum(si.X(1:b-1)) + t*si.w <= x) &&...
                    (x <= si.xi + sum(si.X(1:b)) + t*si.w) &&...
                    ( b <= si.b_max))
                
                array(1,si.n_max_us + si.n_max_ds + 1 : si.n_max_us + si.n_max_ds + b - 1 ) = -si.X(1:b-1);
                array(1,si.n_max_us + si.n_max_ds + b) = (t*si.w - x + sum(si.X(1:b-1)) + si.xi);
                array(1,si.size_row) = si.k_m*t*si.w;
                return
            end
            
            if ( (si.xi + sum(si.X(1:b)) + t*si.w < x) &&...
                    (x <= si.xi + sum(si.X(1:b)) + t*si.v) &&...
                    (b <= si.b_max))
                
                array(1,si.n_max_us + si.n_max_ds + 1:si.n_max_us + si.n_max_ds+b-1) = -si.X(1:b-1);
                array(1,si.n_max_us + si.n_max_ds + b) = -si.X(b);
                array(1,si.size_row) = -(si.k_c*(t*si.w - x + sum(si.X(1:b)) + si.xi) - si.k_m*t*si.w);
                return
            else
                array = []; %Return a "NULL" value
                return
            end
        end
        
        
        %==========================================================================        
        %Equations of function Mups1, explicit solutions of density constraints applying when p<pc        
        
        function [array] = mups1 (si,u,t,x)
            
            array = zeros(1,si.size_row);
            
            if ( (si.x_min_u(u) + (t-si.t_u(u))*si.v <= x) &&...
                 (x <=si.x_max_u(u) + (t-si.t_u(u))*si.v) &&...
                 (t>=si.t_u(u)) && ( u <= si.u_max))
                
                array(1,si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + u ) = 1;
                array(1,si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + si.u_max + u) =...
                    (t-si.t_u(u))*si.v - x + si.x_min_u(u);
                return
            end
            
            if ( (si.x_min_u(u) + (t-si.t_u(u))*si.w <= x) &&...
                 (x <= si.x_min_u(u) + (t-si.t_u(u))*si.v) &&...
                 (t>=si.t_u(u)) && (u <= si.u_max))
                array(1,si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + u ) = 1;
                array(1,si.size_row) = -si.k_c*((t-si.t_u(u))*si.v - x + si.x_min_u(u));
                return
            else
                
                array = []; %Return a "NULL"
                return
                
            end
            
        end
      
        %==========================================================================        
        %Equations of function Mups2, explicit solutions of density constraints applying when p>pc        
        
        function [array] = mups2 (si,u,t,x)
            
            array = zeros(1,si.size_row);
            
            if ( (si.x_min_u(u) + (t-si.t_u(u))*si.w <= x) &&...
                 (x <=si.x_max_u(u) + (t-si.t_u(u))*si.w) &&...
                 (t>=si.t_u(u)) && ( u <= si.u_max))
                
                array(1,si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + u ) = 1;
                array(1,si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + si.u_max + u) =...
                        (t-si.t_u(u))*si.w - x + si.x_min_u(u);
                array(1,si.size_row) = si.k_m*(t-si.t_u(u))*si.w;
                return
            end
            
            if ( (si.x_max_u(u) + (t-si.t_u(u))*si.w <= x) &&...
                 (x <= si.x_max_u(u) + (t-si.t_u(u))*si.v) &&...
                 (t>=si.t_u(u)) && (u <= si.u_max))
                
                array(1,si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + u ) = 1;
                array(1,si.size_row) = -(si.k_c*((t-si.t_u(u))*si.w - x +...
                                si.x_max_u(u)) - si.k_m*(t-si.t_u(u))*si.w);
                return
            else
                array = []; %Return a "NULL"
                return
            end
            
        end
        
        %==========================================================================
        %Equations that define the upstream condition
        
        function[array] = gamma(si,t,x)
            
            array = zeros(1,si.size_row);  %Initialize the array
            
            if t < si.start_time || t > si.end_time
                array = [];
                return
            end
            
            n = sum(si.T_us_cum <= t);  % in step n interval [ , )
            if t == si.end_time     % the final time point
                n = n-1;    % in the last step interval
            end
            
            
            if( (abs(x - si.xi)<0.00001) && (si.T_us_cum(n)<=t) &&...
                (t<=si.T_us_cum(n+1)) && (n <= si.n_max_us) )
                array(1,1:n-1) = si.T_us(1:n-1);
                array(1,n) = t-si.T_us_cum(n);
                return
            else
                array = []; %return a "NULL"
                return
            end
        end
        
        %==========================================================================
        %Equations that define the downstream condition
        
        function [array] = beta(si,t,x)
            array = zeros(1,si.size_row);
            
            if t < si.start_time || t > si.end_time
                array = [];
                return
            end
            
            n = sum(si.T_ds_cum <= t);  % in step n interval [ , )
            if t == si.end_time     % the final time point
                n = n-1;    % in the last step interval
            end
            
            if( (abs(x-si.chi)<0.01) && (si.T_ds_cum(n) <= t) &&...
                (t<= si.T_ds_cum(n+1)) && (n <= si.n_max_ds) )
                array(1,si.n_max_us + si.n_max_ds +1:...
                      si.n_max_us + si.n_max_ds + si.b_max) = -si.X;   %Initial number of vehicles
                array(1,si.n_max_us +1:si.n_max_us + n-1) = si.T_ds(1:n-1);
                array(1,si.n_max_us + n ) = t-si.T_ds_cum(n);
                return
            else
                array = [];
                return
            end
        end
        
        %==========================================================================
        %Equations that define the internal condition
        
        function [array] = mu(si, m ,t ,x)
            
            array = zeros(1,si.size_row);
            if ( (abs(x - (si.v_meas(m)*(t-si.t_min(m)) + si.x_min(m))) < 0.01) &&...
                    (t>= si.t_min(m)) && (t<= si.t_max(m)) )
                array(1, si.n_max_us + si.n_max_ds + si.b_max + m) = 1;
                array(1, si.n_max_us + si.n_max_ds + si.b_max + si.m_max + m ) = t - si.t_min(m);
                return
            else
                array = [];
                return
            end
        end
        
        %==========================================================================
        %Equations that define the density condition       

        %new variables x_min_u, x_max_u, t_u, u_max, m_max_vm0
      
        function [array] = ups(si, u ,t ,x)
            
            if ( (abs(t-si.t_u(u)) < 0.0001) && (x>= si.x_min_u(u)) && (x<= si.x_max_u(u)) )
                array = zeros(1,si.size_row);
                array(1, si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + u) = 1;
                array(1, si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + si.u_max + u ) = -(x - si.x_min_u(u+1));
                return
            else
                
                array = [];
                return
                
            end
            
        end
      
        %==========================================================================
        %Equations that define the initial condition
        
        function[array] = tau(si, b, t, x)
            if( (abs(t-si.t0)<0.01) &&  ( sum(si.X(1:b)) <= x) && ( x <=sum(si.X(1:b+1))) && (b <= si.b_max) )
                array = zeros(1,si.size_row);
                array(1, 2*(si.n_max+1)+1:2*(si.n_max+1)+b) = -si.X(1:b);
                array(1,2*(si.n_max+1) +  b +1) = -(x-sum(si.X(1:b)));
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        %==========================================================================
        %Function to create the model constraints
        % output: the model constraints matrix
        
        function  [list] = setModelMatrix(si)
            
            % first allocate memory for list
            
            list = zeros(100000,si.size_row);
            % initialize rows counts
            rows = 0; 
            
            %==============================================
            % Solution associated with n-th upstream boundary conditions
            for n=1:si.n_max_us
                
                % at downstream points including the begining and ending
                for p=1:si.n_max_ds+1
                    
                    if si.T_us_cum(n) + (si.chi-si.xi)/si.v <= si.T_ds_cum(p) &&...
                       si.T_us_cum(n+1) + (si.chi-si.xi)/si.v >= si.T_ds_cum(p)
                        % point (si.T_ds_cum(p), si.chi)
                        array = si.mgamma(n,si.T_ds_cum(p),si.chi);
                        array2 = si.substractArray(array, si.beta(si.T_ds_cum(p), si.chi));
                        if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                            rows = rows+1;
                            list(rows,:) = array2;
                        end
                    end
                    
                end
                
                % freeflow speed intersection at downstream point
                % point (si.T_us_cum(n) + (si.chi-si.xi)/si.v, si.chi) 
                % where solution >= value condition 
                array = si.mgamma(n, si.T_us_cum(n) + (si.chi-si.xi)/si.v, si.chi);
                array2 = si.substractArray(array, si.beta(si.T_us_cum(n) + (si.chi-si.xi)/si.v, si.chi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % There is no need to pose the inequality for upstream
                % points or initial time points
                
                % internal condition points
                for m=1:si.m_max
                    
                    array = si.mgamma(n,si.t_min(m), si.x_min(m));
                    array2 = si.substractArray(array, si.mu(m, si.t_min(m), si.x_min(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        
                        list(rows,:) = array2;
                    end
                    
                    array = si.mgamma(n,si.t_max(m), si.x_max(m));
                    array2 = si.substractArray(array, si.mu(m, si.t_max(m), si.x_max(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.T_us_cum(n) * si.v - si.v_meas(m) * si.t_min(m) +...
                              si.x_min(m) - si.xi) / (si.v - si.v_meas(m));
                    x_temp = si.x_min(m) + si.v_meas(m)*(t_temp - si.t_min(m));
                    array = si.mgamma(n,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(m, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
          
                end
                
                % density condition points
                for u=1:si.u_max
                    
                    array = si.mgamma(n,si.t_u(u), si.x_min_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    array = si.mgamma(n,si.t_u(u), si.x_max_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mgamma(n,si.t_u(u), si.xi + si.v*(si.t_u(u)-si.T_us_cum(n)));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi +...
                        si.v*(si.t_u(u)-si.T_us_cum(n))));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
     
                end
                
            end
            
            
            %==============================================
            % Solution associated with n-th downstram boundary conditions
            for n=1:si.n_max_ds
                
                % at upstream points
                for p=1:si.n_max_us+1
                    
                    if si.T_ds_cum(n) + (si.xi-si.chi)/si.w <= si.T_us_cum(p) &&...
                            si.T_ds_cum(n+1) + (si.xi-si.chi)/si.w >= si.T_us_cum(p)
                        % point (si.T_us_cum(p), si.xi)
                        array = si.mbeta(n,si.T_us_cum(p), si.xi);
                        array2 = si.substractArray(array, si.gamma(si.T_us_cum(p), si.xi));
                        if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                            rows = rows+1;
                            
                            list(rows,:) = array2;
                        end
                    end
                    
                end
                
                % intersection point at upstream (si.T_ds_cum(n) + (si.xi-si.chi)/si.w,si.xi)
                array = si.mbeta(n,si.T_ds_cum(n) + (si.xi-si.chi)/si.w, si.xi);
                array2 = si.substractArray(array, si.gamma(si.T_ds_cum(n) + (si.xi-si.chi)/si.w, si.xi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    
                    list(rows,:) = array2;
                end
                
                % There is no need to pose the inequality for upstream
                % points or initial time points
                
                % internal condition points
                for m=1:si.m_max
                    array = si.mbeta(n,si.t_min(m), si.x_min(m));
                    array2 = si.substractArray(array, si.mu(m, si.t_min(m), si.x_min(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mbeta(n,si.t_max(m), si.x_max(m));
                    array2 = si.substractArray(array, si.mu(m, si.t_max(m), si.x_max(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.T_ds_cum(n) * si.w - si.v_meas(m) * si.t_min(m) +...
                        si.x_min(m) - si.chi) / (si.w - si.v_meas(m));
                    x_temp = si.x_min(m) + si.v_meas(m)*(t_temp - si.t_min(m));
                    array = si.mbeta(n,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(m, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                end
                
                % density condition points
                for u=1:si.u_max
                    
                    array = si.mbeta(n,si.t_u(u), si.x_min_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mbeta(n,si.t_u(u), si.x_max_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mbeta(n,si.t_u(u), si.chi + si.w*(si.t_u(u)-si.T_ds_cum(n)));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.chi +...
                        si.w*(si.t_u(u)-si.T_ds_cum(n))));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                end
                
            end
            
            
            %==============================================
            % Solution associated with initial conditions
            for k=1:si.b_max
                
                % at upstream points
                for p=1:si.n_max_us+1
                    
                    %TAU1 rho < rho_c
                    
                    % upstream time grid points
                    array = si.mtau1(k,si.T_us_cum(p),si.xi);
                    array2 = si.substractArray(array,si.gamma(si.T_us_cum(p),si.xi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = si.mtau2(k,si.T_us_cum(p),si.xi);
                    array2 = si.substractArray(array,si.gamma(si.T_us_cum(p),si.xi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % at upstream intersection points
                % tau 1
                array = si.mtau1(k,si.start_time + (-sum(si.X(1:k)))/si.w,si.xi);
                array2 = si.substractArray(array,si.gamma(si.start_time +...
                    (-sum(si.X(1:k)))/si.w,si.xi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % tau 2
                array = si.mtau2(k,si.start_time + (-sum(si.X(1:k)))/si.w,si.xi);
                array2 = si.substractArray(array,si.gamma(si.start_time +...
                    (-sum(si.X(1:k)))/si.w,si.xi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % at downstream points
                for p=1:si.n_max_ds+1
                   
                    % Tau1
                    array = si.mtau1(k,si.T_ds_cum(p),si.chi);
                    array2 = si.substractArray(array,si.beta(si.T_ds_cum(p),si.chi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                   
                    %tau2
                    array = si.mtau2(k,si.T_ds_cum(p),si.chi);
                    array2 = si.substractArray(array,si.beta(si.T_ds_cum(p),si.chi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % intersection point at downstream
                array = si.mtau1(k,si.start_time + (si.chi-(sum(si.X(1:k-1))+si.xi))/si.v,si.chi);
                array2 = si.substractArray(array,si.beta(si.start_time +...
                    (si.chi-(sum(si.X(1:k-1))+si.xi))/si.v,si.chi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = si.mtau2(k,si.start_time + (si.chi-(sum(si.X(1:k-1))+si.xi))/si.v,si.chi);
                array2 = si.substractArray(array,si.beta(si.start_time +...
                    (si.chi-(sum(si.X(1:k-1))+si.xi))/si.v,si.chi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % for internal points
                for p=1:si.m_max
                    
                    %TAU 1
                    array = si.mtau1(k,si.t_min(p), si.x_min(p));
                    array2 = si.substractArray(array, si.mu(p, si.t_min(p), si.x_min(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau1(k,si.t_max(p), si.x_max(p));
                    array2 = si.substractArray(array, si.mu(p, si.t_max(p), si.x_max(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - (si.xi + sum(si.X(1:k))) - si.v_meas(p) * si.t_min(p)) ...
                        / (si.v - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mtau1(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - (si.xi + sum(si.X(1:k-1))) - si.v_meas(p) * si.t_min(p))...
                        / (si.w - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mtau1(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - (si.xi + sum(si.X(1:k-1))) - si.v_meas(p) * si.t_min(p))...
                        / (si.v - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mtau1(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - (si.xi + sum(si.X(1:k))) - si.v_meas(p) * si.t_min(p))...
                        / (si.w - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mtau1(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = si.mtau2(k,si.t_min(p), si.x_min(p));
                    array2 = si.substractArray(array, si.mu(p, si.t_min(p), si.x_min(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau2(k,si.t_max(p), si.x_max(p));
                    array2 = si.substractArray(array, si.mu(p, si.t_max(p), si.x_max(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - (si.xi + sum(si.X(1:k))) - si.v_meas(p) * si.t_min(p))...
                        / (si.v - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mtau2(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - (si.xi + sum(si.X(1:k-1))) - si.v_meas(p) * si.t_min(p))...
                        / (si.w - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mtau2(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - (si.xi + sum(si.X(1:k-1))) - si.v_meas(p) * si.t_min(p))...
                        / (si.v - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mtau2(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - (si.xi + sum(si.X(1:k))) - si.v_meas(p) * si.t_min(p))...
                        / (si.w - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mtau2(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % for density points
                for u=1:si.u_max
                    
                    %TAU 1
                    array = si.mtau1(k,si.t_u(u), si.x_min_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau1(k,si.t_u(u), si.x_max_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau1(k,si.t_u(u), si.xi + sum(si.X(1:k-1)) +...
                        (si.t_u(u)-si.start_time)*si.v);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi + sum(si.X(1:k-1)) +...
                        (si.t_u(u)-si.start_time)*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau1(k,si.t_u(u), si.xi + sum(si.X(1:k)) +...
                        (si.t_u(u)-si.start_time)*si.v);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi + sum(si.X(1:k)) +...
                        (si.t_u(u)-si.start_time)*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau1(k,si.t_u(u), si.xi + sum(si.X(1:k-1)) +...
                        (si.t_u(u)-si.start_time)*si.w);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi + sum(si.X(1:k-1)) +...
                        (si.t_u(u)-si.start_time)*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau1(k,si.t_u(u), si.xi + sum(si.X(1:k)) +...
                        (si.t_u(u)-si.start_time)*si.w);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi + sum(si.X(1:k)) +...
                        (si.t_u(u)-si.start_time)*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = si.mtau2(k,si.t_u(u), si.x_min_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau2(k,si.t_u(u), si.x_max_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau2(k,si.t_u(u), si.xi + sum(si.X(1:k-1)) +...
                        (si.t_u(u)-si.start_time)*si.v);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi + sum(si.X(1:k-1)) +...
                        (si.t_u(u)-si.start_time)*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau2(k,si.t_u(u), si.xi + sum(si.X(1:k)) +...
                        (si.t_u(u)-si.start_time)*si.v);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi + sum(si.X(1:k)) +...
                        (si.t_u(u)-si.start_time)*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau2(k,si.t_u(u), si.xi + sum(si.X(1:k-1)) +...
                        (si.t_u(u)-si.start_time)*si.w);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi + sum(si.X(1:k-1)) +...
                        (si.t_u(u)-si.start_time)*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mtau2(k,si.t_u(u), si.xi + sum(si.X(1:k)) +...
                        (si.t_u(u)-si.start_time)*si.w);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.xi + sum(si.X(1:k)) +...
                        (si.t_u(u)-si.start_time)*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
            end
            
            
            %==============================================
            % Solution associated with internal conditions
            for m=1:si.m_max
                
                % upstream points
                for p=1:si.n_max_us+1
                    
                    array = si.mmu(m,si.T_us_cum(p), si.xi);
                    array2 = si.substractArray(array, si.gamma(si.T_us_cum(p), si.xi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.xi - si.x_min(m) + si.w*si.t_min(m)) / (si.w);
                    x_temp = si.xi;
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.gamma(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.xi - si.x_max(m) + si.w*si.t_max(m)) / (si.w);
                    x_temp = si.xi;
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.gamma(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % downstream points
                for p=1:si.n_max_ds+1
                    
                    array = si.mmu(m,si.T_ds_cum(p), si.chi);
                    array2 = si.substractArray(array, si.beta(si.T_ds_cum(p), si.chi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.chi - si.x_min(m) + si.v*si.t_min(m)) / (si.v);
                    x_temp = si.chi;
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.beta(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.chi - si.x_max(m) + si.v*si.t_max(m)) / (si.v);
                    x_temp = si.chi;
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.beta(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % other internal points
                for p=1:si.m_max
                    
                    array = si.mmu(m,si.t_min(p), si.x_min(p));
                    array2 = si.substractArray(array, si.mu(p,si.t_min(p), si.x_min(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mmu(m,si.t_max(p), si.x_max(p));
                    array2 = si.substractArray(array, si.mu(p,si.t_max(p), si.x_max(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(m) - si.x_min(p) + si.v_meas(p) * si.t_min(p) -...
                        si.v_meas(m) * si.t_min(m)) / (si.v_meas(p) - si.v_meas(m));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (si.x_max(m) - si.x_min(p) + si.v_meas(p) * si.t_min(p) -...
                        si.v * si.t_max(m)) / (si.v_meas(p) - si.v);
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (si.x_min(m) - si.x_min(p) + si.v_meas(p) * si.t_min(p) -...
                        si.v * si.t_min(m)) / (si.v_meas(p) - si.v);
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (si.x_max(m) - si.x_min(p) + si.v_meas(p) * si.t_min(p) -...
                        si.v * si.t_max(m)) / (si.v_meas(p) - si.w);
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (si.x_min(m) - si.x_min(p) + si.v_meas(p) * si.t_min(p) -...
                        si.v * si.t_min(m)) / (si.v_meas(p) - si.w);
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mmu(m,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
      
                    
                end
                
                
                % density condition points
                for u=1:si.u_max
                    
                    array = si.mmu(m,si.t_u(u), si.x_min_u(u));
                    array2 = si.substractArray(array, si.ups(u,si.t_u(u), si.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mmu(m,si.t_u(u), si.x_max_u(u));
                    array2 = si.substractArray(array, si.ups(u,si.t_u(u), si.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mmu(m,si.t_u(u), si.x_min(m)+(si.t_u(u)-si.t_min(m))*si.w);
                    array2 = si.substractArray(array, si.ups(u,si.t_u(u),...
                                               si.x_min(m)+(si.t_u(u)-si.t_min(m))*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mmu(m,si.t_u(u), si.x_max(m)+(si.t_u(u)-si.t_max(m))*si.w);
                    array2 = si.substractArray(array, si.ups(u,si.t_u(u),...
                        si.x_max(m)+(si.t_u(u)-si.t_max(m))*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mmu(m,si.t_u(u), si.x_min(m)+(si.t_u(u)-si.t_min(m))*si.v);
                    array2 = si.substractArray(array, si.ups(u,si.t_u(u),...
                        si.x_min(m)+(si.t_u(u)-si.t_min(m))*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mmu(m,si.t_u(u), si.x_max(m)+(si.t_u(u)-si.t_max(m))*si.v);
                    array2 = si.substractArray(array, si.ups(u,si.t_u(u),...
                        si.x_max(m)+(si.t_u(u)-si.t_max(m))*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
            end
            
            
            % Solution associated with density conditions
            for k=1:si.u_max

                % upstream points
                for p=1:si.n_max_us+1
                    
                    %UPS1
                    array = si.mups1(k,si.T_us_cum(p),si.xi);
                    array2 = si.substractArray(array,si.gamma(si.T_us_cum(p),si.xi));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    %UPS2
                    array = si.mups2(k,si.T_us_cum(p),si.xi);
                    array2 = si.substractArray(array,si.gamma(si.T_us_cum(p),si.xi));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % upstream intersection points
                array = si.mups1(k,si.t_u(k) + (si.xi-si.x_max_u(k))/si.w,si.xi);
                array2 = si.substractArray(array,si.gamma(si.t_u(k) +...
                        (si.xi-si.x_max_u(k))/si.w,si.xi));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = si.mups2(k,si.t_u(k) + (si.xi-si.x_max_u(k))/si.w,si.xi);
                array2 = si.substractArray(array,si.gamma(si.t_u(k) +...
                    (si.xi-si.x_max_u(k))/si.w,si.xi));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % downstream points
                for p=1:si.n_max_ds+1
                    
                    %UPS1
                    array = si.mups1(k,si.T_ds_cum(p),si.chi);
                    array2 = si.substractArray(array,si.beta(si.T_ds_cum(p),si.chi));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    %USP2
                    array = si.mups2(k,si.T_ds_cum(p),si.chi);
                    array2 = si.substractArray(array,si.beta(si.T_ds_cum(p),si.chi));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % downstream intersection points
                array = si.mups1(k,si.t_u(k) + (si.chi-si.x_min_u(k))/si.v,si.chi);
                array2 = si.substractArray(array,si.beta(si.t_u(k) +...
                    (si.chi-si.x_min_u(k))/si.v,si.chi));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = si.mups2(k,si.t_u(k) + (si.chi-si.x_min_u(k))/si.v,si.chi);
                array2 = si.substractArray(array,si.beta(si.t_u(k) +...
                    (si.chi-si.x_min_u(k))/si.v,si.chi));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % internal condition points
                for p=1:si.m_max
                    
                    %UPS 1
                    array = si.mups1(k,si.t_min(p), si.x_min(p));
                    array2 = si.substractArray(array, si.mu(p, si.t_min(p), si.x_min(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups1(k,si.t_max(p), si.x_max(p));
                    array2 = si.substractArray(array, si.mu(p, si.t_max(p), si.x_max(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - si.x_max_u(k) - si.v_meas(p) * si.t_min(p) +...
                        si.v*si.t_u(k)) / (si.v - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mups1(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - si.x_min_u(k) - si.v_meas(p) * si.t_min(p) +...
                        si.w*si.t_u(k)) / (si.w - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mups1(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - si.x_min_u(k) - si.v_meas(p) * si.t_min(p) +...
                        si.v*si.t_u(k)) / (si.v - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mups1(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - si.x_max_u(k) - si.v_meas(p) * si.t_min(p) +...
                        si.w*si.t_u(k)) / (si.w - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mups1(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %UPS 2
                    array = si.mups2(k,si.t_min(p), si.x_min(p));
                    array2 = si.substractArray(array, si.mu(p, si.t_min(p), si.x_min(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups2(k,si.t_max(p), si.x_max(p));
                    array2 = si.substractArray(array, si.mu(p, si.t_max(p), si.x_max(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - si.x_max_u(k) - si.v_meas(p) * si.t_min(p) +...
                        si.v*si.t_u(k)) / (si.v - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mups2(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - si.x_min_u(k) - si.v_meas(p) * si.t_min(p) +...
                        si.w*si.t_u(k)) / (si.w - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mups2(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - si.x_min_u(k) - si.v_meas(p) * si.t_min(p) +...
                        si.v*si.t_u(k)) / (si.v - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mups2(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (si.x_min(p) - si.x_max_u(k) - si.v_meas(p) * si.t_min(p) +...
                        si.w*si.t_u(k)) / (si.w - si.v_meas(p));
                    x_temp = si.v_meas(p) * (t_temp - si.t_min(p)) + si.x_min(p);
                    array = si.mups2(k,t_temp, x_temp);
                    array2 = si.substractArray(array, si.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % density condition points
                for u=1:si.u_max
                    
                    %UPS 1
                    array = si.mups1(k,si.t_u(u), si.x_min_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups1(k,si.t_u(u), si.x_max_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups1(k,si.t_u(u), si.x_min_u(k) + (si.t_u(u)-si.t_u(k))*si.v);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(k) +...
                        (si.t_u(u)-si.t_u(k))*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups1(k,si.t_u(u), si.x_max_u(k) + (si.t_u(u)-si.t_u(k))*si.v);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(k) +...
                        (si.t_u(u)-si.t_u(k))*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups1(k,si.t_u(u), si.x_min_u(k) + (si.t_u(u)-si.t_u(k))*si.w);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(k) +...
                        (si.t_u(u)-si.t_u(k))*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups1(k,si.t_u(u), si.x_max_u(k) + (si.t_u(u)-si.t_u(k))*si.w);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(k) +...
                        (si.t_u(u)-si.t_u(k))*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    % %UPS 2
                    array = si.mups2(k,si.t_u(u), si.x_min_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups2(k,si.t_u(u), si.x_max_u(u));
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups2(k,si.t_u(u), si.x_min_u(k) + (si.t_u(u)-si.t_u(k))*si.v);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(k) +...
                        (si.t_u(u)-si.t_u(k))*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups2(k,si.t_u(u), si.x_max_u(k) + (si.t_u(u)-si.t_u(k))*si.v);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(k) +...
                        (si.t_u(u)-si.t_u(k))*si.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups2(k,si.t_u(u), si.x_min_u(k) + (si.t_u(u)-si.t_u(k))*si.w);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_min_u(k) +...
                        (si.t_u(u)-si.t_u(k))*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = si.mups2(k,si.t_u(u), si.x_max_u(k) + (si.t_u(u)-si.t_u(k))*si.w);
                    array2 = si.substractArray(array, si.ups(u, si.t_u(u), si.x_max_u(k) +...
                        (si.t_u(u)-si.t_u(k))*si.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
            end
            
            % truncate the list matrix to the first rows with
            % inequality information.
            list = list(1:rows,:);
            
            
        end
        
        %==========================================================================
        %Function to create the data constraints
        % output: the data constraints matrix
        
        function [list2] = setDataMatrix(si)
            
            global FLUX
            
            %Use sparse matrix to speed up computation
            list2 = zeros(10000,si.size_row);
            rows = 0;   %initialize row counts
            
            % upstream data constraints
            if (~isempty(si.qin_meas))
                
                for n=1:si.n_max_us
                    if ~isnan(si.qin_meas(n, FLUX))
                        array = zeros(1,si.size_row);
                        array(1,n) = 1;
                        array(1,si.size_row) = si.qin_meas(n,FLUX)*(1-si.e_max);
                        
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for n=1:si.n_max_us
                    if ~isnan(si.qin_meas(n, FLUX))
                    array = zeros(1,si.size_row);
                    array(1,n) = -1;
                    array(1,si.size_row) = -si.qin_meas(n,FLUX)*(1+si.e_max);
                    rows = rows+1;
                    list2(rows,:)=array;
                    end
                end
                
            end
            
            % downstream data constraints
            if (~isempty(si.qout_meas) )
                
                for n=1:si.n_max_ds
                    if ~isnan(si.qout_meas(n,FLUX))
                    array = zeros(1,si.size_row);
                    array(1,si.n_max_us + n) = 1;
                    array(1,si.size_row) = si.qout_meas(n,FLUX)*(1-si.e_max);
                    rows = rows+1;
                    list2(rows,:)=array;
                    end
                end
                
                for n=1:si.n_max_ds
                    if ~isnan(si.qout_meas(n,FLUX))
                    array = zeros(1,si.size_row);
                    array(1,si.n_max_ds + n) = -1;
                    array(1,si.size_row) = -si.qout_meas(n,FLUX)*(1+si.e_max);
                    rows = rows+1;
                    list2(rows,:)=array;
                    end
                end
                
            end
            
            % Initial data constraints
            if(~isempty(si.rho_ini) )
                
                for n=1:si.b_max
                    if ~isnan(si.rho_ini(n))
                        array = zeros(1,si.size_row);
                        array(1,si.n_max_us + si.n_max_ds + n) = 1;
                        array(1,si.size_row) = si.rho_ini(n)*(1-si.e_max);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for n=1:si.b_max
                    if ~isnan(si.rho_ini(n))
                        array = zeros(1,si.size_row);
                        array(1,si.n_max_us + si.n_max_ds + n) = -1;
                        array(1,si.size_row) = -si.rho_ini(n)*(1+si.e_max);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
            end
            
            % internal data constraints for rate of change r
            if(~isempty(si.r_meas))
                
                for m = 1:si.m_max
                    if ~isnan(si.r_meas(m))
                        array = zeros(1,si.size_row);
                        array(1,si.n_max_us + si.n_max_ds + si.b_max + si.m_max + m) = 1;
                        array(1,si.size_row) = si.r_meas(m);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for m = 1:si.m_max
                    if ~isnan(si.r_meas(m))
                        array = zeros(1,si.size_row);
                        array(1,si.n_max_us + si.n_max_ds + si.b_max + si.m_max + m) = -1;
                        array(1,si.size_row) = - si.r_meas(m);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
            end
            
            % density data constraints for density rho_u
            if(~isempty(si.dens_meas))
                
                for u = 1:si.u_max
                    if ~isnan(si.dens_meas(u))
                        array = zeros(1,si.size_row);
                        array(1,si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + ...
                                si.u_max + u) = 1;
                        array(1,si.size_row) = si.dens_meas(u)*(1-si.e_max);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for u = 1:si.u_max
                    if ~isnan(si.dens_meas(u))
                        array = zeros(1,si.size_row);
                        array(1,si.n_max_us + si.n_max_ds + si.b_max + 2*si.m_max + ...
                                si.u_max + u) = -1;
                        array(1,si.size_row) = -si.dens_meas(u)*(1+si.e_max);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
            end
            
            
            % truncate matrix
            list2 = list2(1:rows,:);
            
            
        end
        
        %==========================================================================
        %Function to create the MILP constraints
        % The following functions need to be modified to enable different
        % discretization in up/downstream boundaries.
        
%         function [list3] = setMatrix3(si)
%             
%             list3 = zeros(0,0);
%             
%             %Binary variables useful counter as a reference
%             countB=0;
%             Cmax = 500000; 
%             
%             % nb: number of binary variables per point
%             nb = ceil(log2(2*(si.b_max+1)+2*(si.n_max+1)+(si.m_max+1)+(si.u_max+1)));
%             
%             % comb_mat: combination matrix for possible binary combinations
%             comb_mat = zeros(2^nb,nb);
%             
%             for i = 1:2^nb
%                 comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
%             end
%             
%             for k = 0:si.m_max
%                 
%                 %x_min and t_min point------------------
%                 if (k==0 || (si.x_min(k+1)~=si.x_max(k) && si.t_min(k+1)~=si.t_max(k)))
%                     tempM = zeros(0,0);
%                     %Upper bound constraints
%                     %L1<=M1
%                     %L1<=M2
%                     temp_array = si.mu(k, si.t_min(k+1), si.x_min(k+1));
%                     
%                     for block = 1:size(si.indiced_rho_ini,1) % For initial
%                         
%                         if (si.indiced_rho_ini(block,3) == 0)
%                             
%                             for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
%                                 %-1 to be consistent with other functions si.mtau...
%                                 arraym = si.mtau1(b,si.t_min(k+1), si.x_min(k+1));
%                                 array2 = si.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                     break;  % Only need to consider the first one
%                                 end
%                             end
%                             
%                         elseif (si.indiced_rho_ini(block,3) == 1)
%                             
%                             for b = si.indiced_rho_ini(block,2)-1:-1:si.indiced_rho_ini(block,1)-1
%                                 % For k>k_c case, From the last one
%                                 arraym = si.mtau2(b,si.t_min(k+1), si.x_min(k+1));
%                                 array2 = si.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                     break;
%                                 end
%                             end
%                             
%                         elseif (isnan(si.indiced_rho_ini(block,3)))
%                             
%                             for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
%                                 % For undefined initial blocks, check all
%                                 arraym = si.mtau1(b,si.t_min(k+1), si.x_min(k+1));
%                                 array2 = si.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                 end
%                                 
%                                 arraym = si.mtau2(b,si.t_min(k+1), si.x_min(k+1));
%                                 array2 = si.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                 end
%                             end
%                             
%                         end
%                     end 
%                     
%                     for n=0:si.n_max    % For Downstream
%                         
%                         arraym = si.mbeta(n,si.t_min(k+1), si.x_min(k+1));
%                         array2 = si.substractArray(arraym,temp_array);
%                         if(~isempty(array2) && arraym(1,si.n_max +1 + n+1) ~= si.Tvec(n+1))
%                             % The second condition ensure we only consider
%                             % the solution at the characteristic domain
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             break;  %only one
%                         end
%                         
%                     end
%                     
%                     for n=0:si.n_max    % For upstream
%                         
%                         arraym = si.mgamma(n,si.t_min(k+1), si.x_min(k+1));
%                         array2 = si.substractArray(arraym, temp_array);
%                         if(~isempty(array2) && arraym(1,si.size_row)==0)
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             break;
%                         end
%                         
%                     end
%                     
%                     for n=0:si.m_max    % For internal
%                         
%                         if(n~=k)
%                             arraym = si.mmu(n,si.t_min(k+1), si.x_min(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                         end
%                         
%                     end
%                     
%                     index_u = find(si.t_u~=si.t0+si.Tcum(si.n_max+2));
%                     for i=1:length(index_u) % For density
%                 
%                         n = index_u(i)-1;
%                         
%                         arraym = si.mups1(n,si.t_min(k+1), si.x_min(k+1));
%                         array2 = si.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                         
%                         arraym = si.mups2(n,si.t_min(k+1), si.x_min(k+1));
%                         array2 = si.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                         
%                     end
%                     
%                     rowsM = size(tempM,1);
%                     
%                     %Define all the possible combinations according to the solutions that apply
%                     %at the specified point
%                     
%                     %Lower bound constraints
%                     %C*b1 +L1>=M1
%                     %C*b1 +L1>=M2
%                     
%                     for i = 1:rowsM
%                         
%                         array = temp_array;
%                         
%                         %Decode the binary combinations
%                         
%                         for counter=1:si.nb_min(k+1)
%                             if comb_mat(i,nb-si.nb_min(k+1)+counter) == 1
%                                 array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) +2*(si.u_max+1)+ countB+ counter) = -Cmax;
%                                 
%                             elseif comb_mat(i,nb-si.nb_min(k+1)+counter) == 0
%                                 array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) +2*(si.u_max+1)+ countB+ counter) = Cmax;
%                                 
%                             end
%                         end
%                         array(1,si.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
%                         
%                         % Add constraint to the MILP matrix
%                         
%                         array2 = si.substractArray(array,tempM(i,:));
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                     end
%                     
%                     %Define the last constraint to bound the binary combination
%                     
%                     array = zeros(1,si.size_row);
%                     
%                     for counter=1:si.nb_min(k+1)
%                         
%                         array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) + 2*(si.u_max+1) + countB + counter) = -2^(si.nb_min(k+1)-counter);
%                         
%                     end
%                     
%                     array(1,si.size_row) = -(rowsM-1); %RHS (maximum possible value of the binary comb)
%                     rows = size(list3,1);
%                     list3(rows+1,:) = array;
%                     
%                     countB = countB+si.nb_min(k+1);
%                     
%                 end
%                 
%                 %x_max and t_max point-------------------------
%                 
%                 tempM = zeros(0,0);
%                 
%                 %Upper bound constraints
%                 %L1<=M1
%                 %L1<=M2
%                 
%                 temp_array = si.mu(k, si.t_max(k+1), si.x_max(k+1));
%                 
%                 for block = 1:size(si.indiced_rho_ini,1) % For initial
%                     
%                     if (si.indiced_rho_ini(block,3) == 0)
%                         
%                         for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
%                             %-1 to be consistent with other functions si.mtau...
%                             arraym = si.mtau1(b,si.t_max(k+1), si.x_max(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                                 break;  % Only need to consider the first one
%                             end
%                         end
%                         
%                     elseif (si.indiced_rho_ini(block,3) == 1)
%                         
%                         for b = si.indiced_rho_ini(block,2)-1:-1:si.indiced_rho_ini(block,1)-1
%                             % For k>k_c case, From the last one
%                             arraym = si.mtau2(b,si.t_max(k+1), si.x_max(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                                 break;
%                             end
%                         end
%                         
%                     elseif isnan(si.indiced_rho_ini(block,3))
%                         
%                         for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
%                             % For undefined initial blocks, check all
%                             arraym = si.mtau1(b,si.t_max(k+1), si.x_max(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                             
%                             arraym = si.mtau2(b,si.t_max(k+1), si.x_max(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                         end
%                         
%                     end
%                     
%                 end
%                 
%                 for n=0:si.n_max
%                     arraym = si.mbeta(n,si.t_max(k+1), si.x_max(k+1));
%                     array2 = si.substractArray(arraym,temp_array);
%                     if(~isempty(array2) && array(1,si.n_max+1 + n+1) ~= si.Tvec(n+1))
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                         break;
%                     end
%                 end
%                 
%                 for n=0:si.n_max
%                     arraym = si.mgamma(n,si.t_max(k+1), si.x_max(k+1));
%                     array2 = si.substractArray(arraym, temp_array);
%                     if(~isempty(array2)  && arraym(1,si.size_row)==0)
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                         break;
%                     end
%                     
%                 end
%                 
%                 for n=0:si.m_max
%                     
%                     
%                     if(n~=k)
%                         arraym = si.mmu(n,si.t_max(k+1), si.x_max(k+1));
%                         array2 = si.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             
%                         end
%                         
%                     end
%                     
%                     
%                 end
%                 
%                 index_u = find(si.t_u~=si.t0+si.Tcum(si.n_max+2));
%                 for i=1:length(index_u) % For density
%                     
%                     n = index_u(i)-1;
%                     
%                     arraym = si.mups1(n,si.t_max(k+1), si.x_max(k+1));
%                     array2 = si.substractArray(arraym,temp_array);
%                     if(~isempty(array2))
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                     end
%                     
%                     arraym = si.mups2(n,si.t_max(k+1), si.x_max(k+1));
%                     array2 = si.substractArray(arraym,temp_array);
%                     if(~isempty(array2))
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                     end
%                     
%                 end
%                 
%                 rowsM = size(tempM,1);
%                 
%                 %Define all the possible combinations according to the solutions that apply
%                 %at the specified point
%                 
%                 %Lower bound constraints
%                 %C*b1 +L1>=M1
%                 %C*b1 +L1>=M2
%                 
%                 for i = 1:rowsM
%                     
%                     array = temp_array;
%                     
%                     %Decode the binary combinations
%                     for counter=1:si.nb_max(1,k+1)
%                         if comb_mat(i,nb-si.nb_max(1,k+1)+counter) == 1
%                             array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) +2*(si.u_max+1)+ countB + counter) = -Cmax;
%                             
%                         elseif comb_mat(i,nb-si.nb_max(1,k+1)+counter) == 0
%                             array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) +2*(si.u_max+1)+ countB + counter) = Cmax;
%                             
%                         end
%                     end
%                     array(1,si.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
%                     
%                     % Add constraint to the matrix constraints
%                     
%                     array2 = si.substractArray(array,tempM(i,:));
%                     rows = size(list3,1);
%                     list3(rows+1,:) = array2;
%                 end
%                 
%                 %Define the last constraint to bound the binary combination
%                 
%                 array = zeros(1,si.size_row);
%                 
%                 % Define constraint matrix elements for binary terms
%                 for counter=1:si.nb_max(1,k+1)
%                     
%                     array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) + 2*(si.u_max+1)+ countB + counter) = -2^(si.nb_max(1,k+1)-counter);
%                     
%                 end
%                 
%                 array(1,si.size_row) = -(rowsM-1); %RHS
%                 rows = size(list3,1);
%                 list3(rows+1,:) = array;
%                 
%                 countB = countB + si.nb_max(k+1);
%             end
% 
%         end
%         
%         
%         function [list3] = setMatrix3u(si)
%             
%             list3 = zeros(0,0);
%             
%             %Binary variables useful counter as a reference
%             countB=0;
%             Cmax = 500000; 
%             
%             % nb: number of binary variables per point
%             nb = ceil(log2(2*(si.b_max+1)+2*(si.n_max+1)+(si.m_max+1)+(si.u_max+1)));
%             
%             % comb_mat: combination matrix for possible binary combinations
%             comb_mat = zeros(2^nb,nb);
%             
%             for i = 1:2^nb
%                 comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
%             end
%             
%             % Binary inequalities induced by density conditions
%             for k = 0:si.u_max
%                  %x_min_u and t_u point------------------
%                 if (k==0 || (si.x_min_u(k+1)~=si.x_max_u(k) && si.t_u(k+1)~=si.t_u(k)))
%                     tempM = zeros(0,0);
%                     %Upper bound constraints
%                     %L1<=M1
%                     %L1<=M2
%                     temp_array = si.ups(k, si.t_u(k+1), si.x_min_u(k+1));
%                     
%                     for block = 1:size(si.indiced_rho_ini,1) % For initial
%                         
%                         if (si.indiced_rho_ini(block,3) == 0)
%                             
%                             for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
%                                 %-1 to be consistent with other functions si.mtau...
%                                 arraym = si.mtau1(b,si.t_u(k+1), si.x_min_u(k+1));
%                                 array2 = si.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                     break;  % Only need to consider the first one
%                                 end
%                             end
%                             
%                         elseif (si.indiced_rho_ini(block,3) == 1)
%                             
%                             for b = si.indiced_rho_ini(block,2)-1:-1:si.indiced_rho_ini(block,1)-1
%                                 % For k>k_c case, From the last one
%                                 arraym = si.mtau2(b,si.t_u(k+1), si.x_min_u(k+1));
%                                 array2 = si.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                     break;
%                                 end
%                             end
%                             
%                         elseif (isnan(si.indiced_rho_ini(block,3)))
%                             
%                             for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
%                                 % For undefined initial blocks, check all
%                                 arraym = si.mtau1(b,si.t_u(k+1), si.x_min_u(k+1));
%                                 array2 = si.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                 end
%                                 
%                                 arraym = si.mtau2(b,si.t_u(k+1), si.x_min_u(k+1));
%                                 array2 = si.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                 end
%                             end
%                             
%                         end
%                     end 
%                     
%                     for n=0:si.n_max    % For Downstream
%                         
%                         arraym = si.mbeta(n,si.t_u(k+1), si.x_min_u(k+1));
%                         array2 = si.substractArray(arraym,temp_array);
%                         if(~isempty(array2) && arraym(1,si.n_max + n+2) ~= si.Tvec(n+1))
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             break;  %only one
%                         end
%                         
%                     end
%                     
%                     for n=0:si.n_max    % For upstream
%                         
%                         arraym = si.mgamma(n,si.t_u(k+1), si.x_min_u(k+1));
%                         array2 = si.substractArray(arraym, temp_array);
%                         if(~isempty(array2) && arraym(1,si.size_row)==0)
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             break;
%                         end
%                         
%                     end
%                     
%                     for n=0:si.m_max    % For internal
%                         
%                         arraym = si.mmu(n,si.t_u(k+1), si.x_min_u(k+1));
%                         array2 = si.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                         
%                     end
%                     
%                     index_u = find(si.t_u <= si.t0+si.Tcum(si.n_max+2));
%                     for i=1:length(index_u) % For density
%                 
%                         n = index_u(i)-1;
%                         
%                         
%                             arraym = si.mups1(n,si.t_u(k+1), si.x_min_u(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                             
%                             arraym = si.mups2(n,si.t_u(k+1), si.x_min_u(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                        
%                         
%                     end
%                     
%                     rowsM = size(tempM,1);
%                     
%                     %Define all the possible combinations according to the solutions that apply
%                     %at the specified point
%                     
%                     %Lower bound constraints
%                     %C*b1 +L1>=M1
%                     %C*b1 +L1>=M2
%                     
%                     for i = 1:rowsM
%                         
%                         array = temp_array;
%                         
%                         %Decode the binary combinations
%                         
%                         for counter=1:si.nb_min_u(k+1)
%                             if comb_mat(i,nb-si.nb_min_u(k+1)+counter) == 1
%                                 array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) +2*(si.u_max+1) +sum(si.nb_min) +sum(si.nb_max) + countB+ counter) = -Cmax;
%                                 
%                             elseif comb_mat(i,nb-si.nb_min_u(k+1)+counter) == 0
%                                 array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) +2*(si.u_max+1)+sum(si.nb_min) +sum(si.nb_max)+ countB+ counter) = Cmax;
%                                 
%                             end
%                         end
%                         array(1,si.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
%                         
%                         % Add constraint to the MILP matrix
%                         
%                         array2 = si.substractArray(array,tempM(i,:));
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                     end
%                     
%                     %Define the last constraint to bound the binary combination
%                     array = zeros(1,si.size_row);
%                     
%                     for counter=1:si.nb_min_u(k+1)
%                         
%                         array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) +2*(si.u_max+1)+sum(si.nb_min)+sum(si.nb_max)+ countB + counter) = -2^(si.nb_min_u(k+1)-counter);
%                         
%                     end
%                     
%                     array(1,si.size_row) = -(rowsM-1); %RHS (maximum possible value of the binary comb)
%                     rows = size(list3,1);
%                     list3(rows+1,:) = array;
%                     
%                     countB = countB+si.nb_min_u(k+1);
%                     
%                 end
%                 
%                 %x_max_u and t_u point-------------------------
%                 
%                 tempM = zeros(0,0);
%                 
%                 %Upper bound constraints
%                 %L1<=M1
%                 %L1<=M2
%                 
%                 temp_array = si.ups(k, si.t_u(k+1), si.x_max_u(k+1));
%                 
%                 for block = 1:size(si.indiced_rho_ini,1) % For initial
%                     
%                     if (si.indiced_rho_ini(block,3) == 0)
%                         
%                         for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
%                             %-1 to be consistent with other functions si.mtau...
%                             arraym = si.mtau1(b,si.t_u(k+1), si.x_max_u(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                                 break;  % Only need to consider the first one
%                             end
%                         end
%                         
%                     elseif (si.indiced_rho_ini(block,3) == 1)
%                         
%                         for b = si.indiced_rho_ini(block,2)-1:-1:si.indiced_rho_ini(block,1)-1
%                             % For k>k_c case, From the last one
%                             arraym = si.mtau2(b,si.t_u(k+1), si.x_max_u(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                                 break;
%                             end
%                         end
%                         
%                     elseif isnan(si.indiced_rho_ini(block,3))
%                         
%                         for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
%                             % For undefined initial blocks, check all
%                             arraym = si.mtau1(b,si.t_u(k+1), si.x_max_u(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                             
%                             arraym = si.mtau2(b,si.t_u(k+1), si.x_max_u(k+1));
%                             array2 = si.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                         end
%                         
%                     end
%                     
%                 end
%                 
%                 for n=0:si.n_max
%                     arraym = si.mbeta(n,si.t_u(k+1), si.x_max_u(k+1));
%                     array2 = si.substractArray(arraym,temp_array);
%                     if(~isempty(array2) && array(1,si.n_max + n+2) ~= si.Tvec(n+1))
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                         break;
%                     end
%                 end
%                 
%                 for n=0:si.n_max
%                     arraym = si.mgamma(n,si.t_u(k+1), si.x_max_u(k+1));
%                     array2 = si.substractArray(arraym, temp_array);
%                     if(~isempty(array2)  && arraym(1,si.size_row)==0)
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                         break;
%                     end
%                     
%                 end
%                 
%                 for n=0:si.m_max
%                     
%                     arraym = si.mmu(n,si.t_u(k+1), si.x_max_u(k+1));
%                     array2 = si.substractArray(arraym,temp_array);
%                     if(~isempty(array2))
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                     end
%                     
%                 end
%                 
%                 index_u = find(si.t_u <= si.t0+si.Tcum(si.n_max+2));
%                 for i=1:length(index_u) % For density
%                     
%                     n = index_u(i)-1;
%                     
%                
%                         arraym = si.mups1(n,si.t_u(k+1), si.x_max_u(k+1));
%                         array2 = si.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                         
%                         arraym = si.mups2(n,si.t_u(k+1), si.x_max_u(k+1));
%                         array2 = si.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                  
%                     
%                 end
%                 
%                 rowsM = size(tempM,1);
%                 
%                 %Define all the possible combinations according to the solutions that apply
%                 %at the specified point
%                 
%                 %Lower bound constraints
%                 %C*b1 +L1>=M1
%                 %C*b1 +L1>=M2
%                 
%                 for i = 1:rowsM
%                     
%                     array = temp_array;
%                     
%                     %Decode the binary combinations
%                     for counter=1:si.nb_max_u(1,k+1)
%                         if comb_mat(i,nb-si.nb_max_u(1,k+1)+counter) == 1
%                             array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) +  2*(si.u_max+1) + sum(si.nb_min) + sum(si.nb_max) +countB + counter) = -Cmax;
%                             
%                         elseif comb_mat(i,nb-si.nb_max_u(1,k+1)+counter) == 0
%                             array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1) + 2*(si.u_max+1) + sum(si.nb_min) + sum(si.nb_max) + countB + counter) = Cmax;
%                             
%                         end
%                     end
%                     array(1,si.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
%                     
%                     % Add constraint to the matrix constraints
%                     
%                     array2 = si.substractArray(array,tempM(i,:));
%                     rows = size(list3,1);
%                     list3(rows+1,:) = array2;
%                 end
%                 
%                 %Define the last constraint to bound the binary combination
%                 
%                 array = zeros(1,si.size_row);
%                 
%                 % Define constraint matrix elements for binary terms
%                 for counter=1:si.nb_max_u(1,k+1)
%                     
%                     array(1,2*(si.n_max + 1) + si.b_max + 1 +2*(si.m_max+1)+ 2*(si.u_max+1) + sum(si.nb_min) + sum(si.nb_max) + countB + counter) = -2^(si.nb_max_u(1,k+1)-counter);
%                     
%                 end
%                 
%                 array(1,si.size_row) = -(rowsM-1); %RHS
%                 rows = size(list3,1);
%                 list3(rows+1,:) = array;
%                 
%                 countB = countB + si.nb_max_u(k+1);
%           
%             end
%             
%         end
%         
%         
%         function [list4] = setMatrixAux(si)
%              list4 = zeros(0,0);
%              
%              for u=0:si.u_max
%                  array = zeros(1,si.size_row);
%                  array(1,2*(si.n_max+1) + (si.b_max+1) + 2*(si.m_max + 1) + (si.u_max+1) + u + 1) = -1;
%                  array(1,2*(si.n_max+1) + (si.b_max+1) + 2*(si.m_max + 1) + 2*(si.u_max+1) + sum(si.nb_min) + sum(si.nb_max)...
%                      + sum(si.nb_min_u) + sum(si.nb_max_u) + u + 1) = si.k_m-si.k_c;
%                  array(1,si.size_row) = -si.k_c;
%                  rows = size(list4,1);
%                  list4(rows+1,:) = array;
%              end
%              
%              for u=0:si.u_max
%                  array = zeros(1,si.size_row);
%                  array(1,2*(si.n_max+1) + (si.b_max+1) + 2*(si.m_max + 1) + (si.u_max+1) + u + 1) = 1;
%                  array(1,2*(si.n_max+1) + (si.b_max+1) + 2*(si.m_max + 1) + 2*(si.u_max+1) + sum(si.nb_min) + sum(si.nb_max)...
%                      + sum(si.nb_min_u) + sum(si.nb_max_u) + u + 1) = -si.k_c;
%                  array(1,si.size_row) = 0;
%                  rows = size(list4,1);
%                  list4(rows+1,:) = array;
%              end
%             
%         end
%             
        
        
        %=========================================================================
        %Get Number Binary variables
        
        function [nb_min,nb_max, nb_min_u, nb_max_u] = getBinaryvar(si)
            
            nb_min = zeros(1,si.m_max);
            nb_max = zeros(1,si.m_max);
            
            nb_min_u = zeros(1,si.u_max);
            nb_max_u = zeros(1,si.u_max);
            
            % Binaries induced by internal conditions
            for k = 1:si.m_max
                % For x_min
                if (k==1 || (si.x_min(k)~=si.x_max(k-1) &&...
                        si.t_min(k)~=si.t_max(k-1)))
                    
                    countM=0;
                    
                    % For initial
                    for block = 1:size(si.indiced_rho_ini,1) 
                        
                        if (si.indiced_rho_ini(block,3) == 0)   % For k<k_c case, only the first one
                            
                            for b = si.indiced_rho_ini(block,1):si.indiced_rho_ini(block,2)
                                array = si.mtau1(b,si.t_min(k), si.x_min(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;  
                                end
                            end
                            
                        elseif (si.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                            
                            for b = si.indiced_rho_ini(block,2):-1:si.indiced_rho_ini(block,1)
                                array = si.mtau2(b,si.t_min(k), si.x_min(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;
                                end
                            end
                            
                        elseif isnan(si.indiced_rho_ini(block,3))
                            
                            for b = si.indiced_rho_ini(block,1):si.indiced_rho_ini(block,2)
                                array = si.mtau1(b,si.t_min(k), si.x_min(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                                
                                array = si.mtau2(b,si.t_min(k), si.x_min(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                            end
                            
                        end
                        
                    end
                    
                    % for downstream
                    % in downstream characteristic domain
                    for n=1:si.n_max_ds
                        
                        array = si.mbeta(n,si.t_min(k), si.x_min(k));
                        if(~isempty(array) && array(1,si.n_max_us + n) ~= si.T_ds(n))
                            % If in characteristic domain
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % for upstream
                    % in upstream characteristic domain
                    for n=1:si.n_max_us
                        
                        array = si.mgamma(n,si.t_min(k), si.x_min(k));
                        if(~isempty(array) && array(1,n) ~= si.T_us(n) )
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % for other internal conditions
                    for n=1:si.m_max
                        
                        if (n~=k)
                            array = si.mmu(n,si.t_min(k), si.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                    
                    % for density conditions
                    index_u = find(si.t_u~= si.end_time);
                    for i=1:length(index_u) % For density
                
                        n = index_u(i);
                        
                        array = si.mups1(n,si.t_min(k), si.x_min(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                        array = si.mups2(n,si.t_min(k), si.x_min(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                    end
                    
                    nb_min(1,k) = ceil(log2(countM));
                    
                end
                
                countM=0;
                
                %Upper bound constraints
                %L1<=M1
                %L1<=M2
                
                % For initial
                for block = 1:size(si.indiced_rho_ini,1)
                    
                    if (si.indiced_rho_ini(block,3) == 0)   % For k<k_c case, only the first one
                        
                        for b = si.indiced_rho_ini(block,1):si.indiced_rho_ini(block,2)
                            array = si.mtau1(b,si.t_min(k), si.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif (si.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                        
                        for b = si.indiced_rho_ini(block,2):-1:si.indiced_rho_ini(block,1)
                            array = si.mtau2(b,si.t_min(k), si.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif isnan(si.indiced_rho_ini(block,3))
                        
                        for b = si.indiced_rho_ini(block,1):si.indiced_rho_ini(block,2)
                            array = si.mtau1(b,si.t_min(k), si.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = si.mtau2(b,si.t_min(k), si.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                    
                end
                
                % for downstream
                % in downstream characteristic domain
                for n=1:si.n_max_ds
                    
                    array = si.mbeta(n,si.t_min(k), si.x_min(k));
                    if(~isempty(array) && array(1,si.n_max_us + n) ~= si.T_ds(n))
                        % If in characteristic domain
                        countM = countM+1;
                        break;
                    end
                    
                end
                
                % for upstream
                % in upstream characteristic domain
                for n=1:si.n_max_us
                    
                    array = si.mgamma(n,si.t_min(k), si.x_min(k));
                    if(~isempty(array) && array(1,n) ~= si.T_us(n) )
                        countM = countM+1;
                        break;
                    end
                    
                end
                
                % for other internal conditions
                for n=0:si.m_max
                    
                    if (n~=k)
                        array = si.mmu(n,si.t_min(k), si.x_min(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                    end
                    
                end
                
                % for density conditions
                index_u = find(si.t_u~= si.end_time);
                for i=1:length(index_u) % For density
                    
                    n = index_u(i);
                    
                    array = si.mups1(n,si.t_min(k), si.x_min(k));
                    if(~isempty(array))
                        countM = countM+1;
                    end
                    
                    array = si.mups2(n,si.t_min(k), si.x_min(k));
                    if(~isempty(array))
                        countM = countM+1;
                    end
                    
                end
                
                nb_max(1,k) = ceil(log2(countM));
            end
            
            
            %Binaries induced by density conditions
            for k = 1:si.u_max
                % For x_min
                if (k==1 || (si.x_min_u(k)~=si.x_max_u(k-1) && si.t_u(k)~=si.t_u(k-1)))
                    
                    countM=0;
                    
                    for block = 1:size(si.indiced_rho_ini,1) % For initial
                        
                        if (si.indiced_rho_ini(block,3) == 0)   % For k<k_c case, only the first one
                            
                            for b = si.indiced_rho_ini(block,1):si.indiced_rho_ini(block,2)
                                array = si.mtau1(b,si.t_u(k), si.x_min_u(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;  
                                end
                            end
                            
                        elseif (si.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                            
                            for b = si.indiced_rho_ini(block,2):-1:si.indiced_rho_ini(block,1)
                                array = si.mtau2(b,si.t_u(k), si.x_min_u(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;
                                end
                            end
                            
                        elseif isnan(si.indiced_rho_ini(block,3))
                            
                            for b = si.indiced_rho_ini(block,1):si.indiced_rho_ini(block,2)
                                array = si.mtau1(b,si.t_u(k), si.x_min_u(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                                
                                array = si.mtau2(b,si.t_u(k), si.x_min_u(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                            end
                            
                        end
                        
                    end
                    
                    % downstream
                    for n=1:si.n_max_ds
                        
                        array = si.mbeta(n,si.t_u(k), si.x_min_u(k));
                        if(~isempty(array) && array(1,si.n_max_us + n) ~= si.T_ds(n))
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % upstream
                    for n=1:si.n_max_us
                        
                        array = si.mgamma(n,si.t_u(k), si.x_min_u(k));
                        if(~isempty(array) && array(1,n) ~= si.T_us(n) )
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % internal conditions
                    for n=1:si.m_max
                        
                        array = si.mmu(n,si.t_u(k), si.x_min_u(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                    end
                    
                    index_u = find(si.t_u~= si.end_time);
                    for i=1:length(index_u) % For density
                
                        n = index_u(i);
                        
                  
                            array = si.mups1(n,si.t_u(k), si.x_min_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = si.mups2(n,si.t_u(k), si.x_min_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                  
                        
                    end
                    
                    nb_min_u(1,k) = ceil(log2(countM));
                    
                end
                
                countM=0;
                
                %Upper bound constraints
                %L1<=M1
                %L1<=M2
                
                for block = 1:size(si.indiced_rho_ini,1) % For initial
                    
                    if (si.indiced_rho_ini(block,3) == 0)
                        
                        for b = si.indiced_rho_ini(block,1):si.indiced_rho_ini(block,2)
                            array = si.mtau1(b,si.t_u(k), si.x_max_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;  % Only need to consider the first one
                            end
                        end
                        
                    elseif (si.indiced_rho_ini(block,3) == 1)
                        
                        for b = si.indiced_rho_ini(block,2):-1:si.indiced_rho_ini(block,1)
                            % For k>k_c case, From the last one
                            array = si.mtau2(b,si.t_u(k), si.x_max_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif isnan(si.indiced_rho_ini(block,3))
                        
                        for b = si.indiced_rho_ini(block,1)-1:si.indiced_rho_ini(block,2)-1
                            array = si.mtau1(b,si.t_u(k), si.x_max_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = si.mtau2(b,si.t_u(k), si.x_max_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                end
                
                % downstream
                for n=1:si.n_max_ds
                    array = si.mbeta(n,si.t_u(k), si.x_max_u(k));
                    if(~isempty(array) && array(1,si.n_max_us + n) ~= si.T_ds(n))
                        countM = countM+1;
                        break;
                    end
                end
                
                % upstream
                for n=1:si.n_max_us
                    array = si.mgamma(n,si.t_u(k), si.x_max_u(k));
                    if(~isempty(array) && array(1,n) ~= si.T_us(n) )
                        countM = countM+1;
                        break;
                    end
                end
                
                % internal
                for n=1:si.m_max
                    
                        array = si.mmu(n,si.t_u(k), si.x_max_u(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
              
                end
                
                index_u = find(si.t_u ~= si.end_time);
                for i=1:length(index_u) % For density
                    
                    n = index_u(i);
                    
                 
                        array = si.mups1(n,si.t_u(k), si.x_max_u(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                        array = si.mups2(n,si.t_u(k), si.x_max_u(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
               
                    
                end
                
                nb_max_u(1,k) = ceil(log2(countM));
            end
      
        end
        
        
    end
end





