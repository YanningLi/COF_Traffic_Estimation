% Yanning Li, Sep 02, 2015
% This class construct a the post solution object
% It uses the Berkeley toolbox to compute the Moskowitz function and
% visualize the result. It also has functions which can check if the
% solution is entropic and update the boundary discretization grid.


classdef postSolution < handle
    
    properties
        
        % column vector, solution obtained from the cplex solver
        x;  
        
        % save a copy of the index of decision variables
        dv_index;
        
        % save a copy of the network structure
        % NOTE: this is not good if we have a large IC BC data set
        net;
        
        % struct, .(linkStr), contains the density for each link
        k;  
        
        % struct, .(linkStr), contains the car ID for each link
        N;
        
        % struct, .(linkStr), discretization in space in meters
        x_mesh_m; 
        
        % the mesh grid is same for all links; different from the boundary
        % discretization.
        t_mesh_s; 
        
        % struct, .(linkStr), fundamental diagram
        fd;
        
        % the resolution we would like to plot in the visualization (mesh)
        dx_res;
        dt_res;
        
        
        % the reference point for sampling M. Note, this is used for
        % determing which boundary conditions we should use to compute the
        % sending and receiving function.
        t_ref;
        
    end
    
    methods
        
        %===============================================================
        % compute the Moskowitz function
        % input:
        %       x: float, dv_index_max x 1, solution from cplex
        %       net: the network object constructed by initNetwork
        %       dv_index: the decision varialbe struct
        %       end_time: end time of the simulation
        %       dx_res/dt_res: float, resolution for visualization
        function self=postSolution(x, net, dv_index, end_time, dx_res, dt_res)
            global INDEX_UP INDEX_DOWN INDEX_INI
            
            self.x = x;
            self.dv_index = dv_index;
            self.net = net;
            
            for link = net.link_labels'
                
                linkStr = sprintf('link_%d',link);
                
                % extract the computed result 
                q_in = x(dv_index.(linkStr)(1,INDEX_UP): dv_index.(linkStr)(2,INDEX_UP));
                q_out = x(dv_index.(linkStr)(1,INDEX_DOWN): dv_index.(linkStr)(2,INDEX_DOWN));
                p_ini = x(dv_index.(linkStr)(1,INDEX_INI): dv_index.(linkStr)(2,INDEX_INI));
                
                % creat a triangular fundamental diagram
                self.fd.(linkStr) = LH_Tfd(net.network_hwy.(linkStr).para_vf,...
                    net.network_hwy.(linkStr).para_w,...
                    net.network_hwy.(linkStr).para_km);
                
                us_position = 0;
                ds_position = net.network_hwy.(linkStr).para_postkm*1000;
                pbEnv = LH_general(self.fd.(linkStr),us_position,ds_position);
                
                
                %===========================================
                % extract initial segment vector
                ini_seg = net.network_hwy.(linkStr).X_grid_cum;
                
                % Berkeley toolbox need ini_seg and p_ini to be a row
                % vector
                if ~isrow(ini_seg)
                    ini_seg = ini_seg';
                end
                if ~isrow(p_ini)
                    p_ini = p_ini';
                end
                
                pbEnv.setIniDens(ini_seg,p_ini);
                
                % set upstream boundary condition
                time_grid = net.network_hwy.(linkStr).T_us_cum;
                if ~isrow(time_grid)
                    time_grid = time_grid';
                end
                if ~isrow(q_in)
                    q_in = q_in';
                end
                pbEnv.setUsFlows(time_grid,q_in);
                
                % set downstream boundary condition
                time_grid = net.network_hwy.(linkStr).T_ds_cum;
                if ~isrow(time_grid)
                    time_grid = time_grid';
                end
                if ~isrow(q_out)
                    q_out = q_out';
                end
                pbEnv.setDsFlows(time_grid,q_out);
                
                %===========================================
                % specify resolution
                nx = floor((ds_position)/dx_res);
                dx = (ds_position)/nx;
                self.dt_res = dt_res;
                self.dx_res = dx_res;
                
                self.x_mesh_m.(linkStr) = 0:dx:ds_position;
                self.t_mesh_s = 0:dt_res:end_time;
                
                xValues = ones(size(self.t_mesh_s'))*(self.x_mesh_m.(linkStr));
                tValues = self.t_mesh_s' * ones(size(self.x_mesh_m.(linkStr)));
                
                %===========================================
                % compute Moskowitz
                result = pbEnv.explSol(tValues,xValues);
                
                self.N.(linkStr) = result{1};
                activeComp = result{2};
                self.k.(linkStr) = pbEnv.density(tValues,xValues,activeComp);
            end
            
        end
        
        
        %===============================================================
        % plot the time-space diagram for the link
        % input:
        %       links: a column vector of link labels that we want to plot
        function plotLinks(self, links)
            
            if strcmp(links,'all')
                links = self.net.link_labels;
            end
            
            
            
            for link = links'
                
                linkStr = sprintf('link_%d',link);
                T_us_cum = self.net.network_hwy.(linkStr).T_us_cum;
                T_ds_cum = self.net.network_hwy.(linkStr).T_ds_cum;
                len_link = self.net.network_hwy.(linkStr).para_postkm*1000;
                
                %===========================================
                %Transformation for better color presentation
                %kc=>0.5km, km=>km
                k_c_tmp = self.net.network_hwy.(linkStr).para_kc;
                k_m_tmp = self.net.network_hwy.(linkStr).para_km;
                k_trans = mapping(self.k.(linkStr), [0 k_c_tmp; k_c_tmp k_m_tmp],...
                    [0 0.5*k_m_tmp; 0.5*k_m_tmp k_m_tmp]);
                
                %===========================================
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 1 scrsz(3) scrsz(4)]);
                title(sprintf('Link No. %d',link),'fontsize',24);
                
                colormap jet
                
                [~, ~] = LH_plot2D(self.t_mesh_s, self.x_mesh_m.(linkStr),...
                    self.N.(linkStr),k_trans, self.fd.(linkStr));
                
                hold on
                % Plot the discretization markers at the upstream
                for i = 1:length(T_us_cum)
                    plot([T_us_cum(i) T_us_cum(i)],...
                         [0 len_link/50],...
                        'k','LineWidth',2);
                end
                
                % plot the discretization markers at the downstream
                for i = 1:length(T_ds_cum)
                    plot([T_ds_cum(i) T_ds_cum(i)],...
                        [len_link-len_link/50 len_link],...
                        'k','LineWidth',2);
                end
                hold off

                set(gca,'fontsize',20)
                xlabel({'time (s)'},'fontsize',24);
                ylabel({'space (m)'},'fontsize',24);
                
            end
               
        end
        
        
        %===============================================================
        % plot the time space diagram for the links at this junction
        % connection: this function concatenate two links 
        % merge: this function plots a 1x2 subplots: 1->3, 2->3
        % diverge: this function plots a 1x2 subplots: 1->2, 1->3
        % input:
        %       juncs: the column vector of junction labels we want to
        %           plot, each junction will be plot in one seperate figure
        function plotJuncs(self, juncs)
            global INDEX_DOWN INDEX_UP
            
            if strcmp(juncs,'all')
                juncs = self.net.junc_labels;
            end
            
            for junc = juncs'
                
                juncStr = sprintf('junc_%d',junc);
                
                if strcmp(self.net.network_junc.(juncStr).type_junc,'merge')
                    %===========================================
                    % compute the total through flow and penalty for
                    % not following the priority rule
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).inlabel(1));
                    q_s1 = self.x(self.dv_index.(linkStr)(1,INDEX_DOWN):...
                                 self.dv_index.(linkStr)(2,INDEX_DOWN));
                             
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).inlabel(2));
                    q_s2 = self.x(self.dv_index.(linkStr)(1,INDEX_DOWN):...
                                 self.dv_index.(linkStr)(2,INDEX_DOWN)); 
                    
                    % the total throughput flow at the junction
                    tt_flow = sum((q_s1+q_s2).*self.net.network_junc.(juncStr).T);
                    
                    % L1 norm penalty for not following the priority rule
                    tt_pen = sum( abs(q_s1 - q_s2*self.net.network_junc.(juncStr).ratio(1)/...
                        self.net.network_junc.(juncStr).ratio(2)).*self.net.network_junc.(juncStr).T );
                    
                    %===========================================
                    % Normalize density
                    % kc=>0.5, km=>1
                    link1 = sprintf('link_%d',self.net.network_junc.(juncStr).inlabel(1));
                    k_c_tmp = self.net.network_hwy.(link1).para_kc;
                    k_m_tmp = self.net.network_hwy.(link1).para_km;
                    kNorm1 = mapping(self.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link2 = sprintf('link_%d',self.net.network_junc.(juncStr).inlabel(2));
                    k_c_tmp = self.net.network_hwy.(link2).para_kc;
                    k_m_tmp = self.net.network_hwy.(link2).para_km;
                    kNorm2 = mapping(self.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link3 = sprintf('link_%d',self.net.network_junc.(juncStr).outlabel);
                    k_c_tmp = self.net.network_hwy.(link3).para_kc;
                    k_m_tmp = self.net.network_hwy.(link3).para_km;
                    kNorm3 = mapping(self.k.(link3),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    kLeft = [kNorm1 kNorm3];
                    kRight = [kNorm2 kNorm3];
                    
                    NLeft = [self.N.(link1) self.N.(link3)];
                    NRight = [self.N.(link2) self.N.(link3)];
                    
                    xScaleLeft = [self.x_mesh_m.(link1)...
                        self.x_mesh_m.(link3) + max(self.x_mesh_m.(link1)) + self.dx_res/1000 ];
                    xScaleRight = [self.x_mesh_m.(link2)...
                        self.x_mesh_m.(link3) + max(self.x_mesh_m.(link2)) + self.dx_res/1000];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    colormap jet
                    
                    [~, ~] = LH_plot2D_junc([self.net.network_hwy.(link1).para_postkm*1000,...
                        self.net.network_hwy.(link2).para_postkm*1000],...
                        self.net.network_junc.(juncStr).T_cum',...
                        length(self.net.network_junc.(juncStr).T), ...
                        self.t_mesh_s, xScaleLeft, xScaleRight,...
                        NLeft, kLeft, NRight, kRight, self.net.network_junc.(juncStr));
                    h = suptitle(sprintf('Total Flow: %f\n Number of steps: %d \n Penalty: %f',...
                        tt_flow, length(self.net.network_junc.(juncStr).T), tt_pen));
                    set(h,'FontSize',24)
                    
                    
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'diverge')
                    %===========================================
                    % compute the total through flow and penalty for
                    % not following the distribution rule
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).outlabel(1));
                    q_r1 = self.x(self.dv_index.(linkStr)(1,INDEX_UP):...
                                 self.dv_index.(linkStr)(2,INDEX_UP));
                             
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).outlabel(2));
                    q_r2 = self.x(self.dv_index.(linkStr)(1,INDEX_UP):...
                                 self.dv_index.(linkStr)(2,INDEX_UP)); 
                                       
                    tt_flow = sum((q_r1+q_r2).*self.net.network_junc.(juncStr).T);
                    % L1 norm penalty
                    tt_pen = sum( abs(q_r1 - q_r2*self.net.network_junc.(juncStr).ratio(1)/...
                        self.net.network_junc.(juncStr).ratio(2)).*self.net.network_junc.(juncStr).T );
                    
                    
                    %===========================================
                    %Normalize density
                    %kc=>0.5, km=>1
                    link1 = sprintf('link_%d',self.net.network_junc.(juncStr).inlabel);
                    k_c_tmp = self.net.network_hwy.(link1).para_kc;
                    k_m_tmp = self.net.network_hwy.(link1).para_km;
                    kNorm1 = mapping(self.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link2 = sprintf('link_%d',self.net.network_junc.(juncStr).outlabel(1));
                    k_c_tmp = self.net.network_hwy.(link2).para_kc;
                    k_m_tmp = self.net.network_hwy.(link2).para_km;
                    kNorm2 = mapping(self.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link3 = sprintf('link_%d',self.net.network_junc.(juncStr).outlabel(2));
                    k_c_tmp = self.net.network_hwy.(link3).para_kc;
                    k_m_tmp = self.net.network_hwy.(link3).para_km;
                    kNorm3 = mapping(self.k.(link3),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    kLeft = [kNorm1 kNorm2];
                    kRight = [kNorm1 kNorm3];
                    
                    NLeft = [self.N.(link1) self.N.(link2)];
                    NRight = [self.N.(link1) self.N.(link3)];
                    
                    xScaleLeft = [self.x_mesh_m.(link1), ...
                        self.x_mesh_m.(link2) + max(self.x_mesh_m.(link1)) + self.dx_res/1000 ];
                    xScaleRight = [self.x_mesh_m.(link1), ...
                        self.x_mesh_m.(link3) + max(self.x_mesh_m.(link1)) + self.dx_res/1000 ];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    colormap jet
                    
                    [~, ~] = LH_plot2D_junc([self.net.network_hwy.(link1).para_postkm*1000,...
                        self.net.network_hwy.(link1).para_postkm*1000],...
                        self.net.network_junc.(juncStr).T_cum',...
                        length(self.net.network_junc.(juncStr).T), ...
                        self.t_mesh_s, xScaleLeft, xScaleRight,...
                        NLeft, kLeft, NRight, kRight,self.net.network_junc.(juncStr));
                    h = suptitle(sprintf('Total Flow: %f\n Number of steps: %d\nPenalty: %f\n',...
                        tt_flow, length(self.net.network_junc.(juncStr).T)), tt_pen);
                    set(h,'FontSize',24)
                    
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'connection')
                    %===========================================
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).outlabel);
                    q_thru = self.x( self.dv_index.(linkStr)(1,INDEX_UP):...
                                    self.dv_index.(linkStr)(2,INDEX_UP));
                    
                    tt_flow = sum((q_thru).*self.net.network_junc.(juncStr).T);
                    
                    %===========================================
                    %Normalize density
                    %kc=>0.5, km=>1
                    link1 = sprintf('link_%d',self.net.network_junc.(juncStr).inlabel);
                    k_c_tmp = self.net.network_hwy.(link1).para_kc;
                    k_m_tmp = self.net.network_hwy.(link1).para_km;
                    kNorm1 = mapping(self.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link2 = sprintf('link_%d',self.net.network_junc.(juncStr).outlabel);
                    k_c_tmp = self.net.network_hwy.(link2).para_kc;
                    k_m_tmp = self.net.network_hwy.(link2).para_km;
                    kNorm2 = mapping(self.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    kComb = [kNorm1 kNorm2];
                    
                    NComb = [self.N.(link1) self.N.(link2)];
                    
                    xScaleComb = [self.x_mesh_m.(link1)...
                        self.x_mesh_m.(link2) + max(self.x_mesh_m.(link1)) + self.dx_res/1000];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    colormap jet
                    
                    [~, ~] = LH_plot2D_junc(self.net.network_hwy.(link1).para_postkm*1000,...
                        self.net.network_junc.(juncStr).T_cum',...
                        length(self.net.network_junc.(juncStr).T), ... 
                        self.t_mesh_s, xScaleComb,...
                        NComb, kComb,self.net.network_junc.(juncStr));
                    h = suptitle(sprintf('Total Flow: %f\n Number of steps: %d\n',...
                        tt_flow, length(self.net.network_junc.(juncStr).T)));
                    set(h,'FontSize',24)
                end
                
            end
            
        end     %end function
        
        
        %===============================================================
        % Analytically check if the solution is entropy solution, that is,
        % the throughput flow is maximized at all time steps.
        % NOTE: This function is meant to find the nonentropic solution due
        % to the boundary discretization. If the entropic objective is not
        % added to the objective function, then the optimal solution may be
        % non-entropic at some steps. In this case, this function will
        % simply return a warning.
        % Intuition: using Moskowitz function, we can analytically compute
        % the sending and receiving function at each boundary
        % discretization point, then compare if we hold back any flow
        % input: 
        %       net: the network object 
        %       entropy_tol: cars, the tolerance for entropy, for instance, 
        %           if the error to entropic solution is only 1 car, we
        %           consider this solution as an entropic solution and stop
        %           updating the discretization grid.
        % output: 
        %       TF: true or false, this solution is an entropic solution?
        %       steps: struct, .(juncStr), the id of the first non entropic
        %           step
        function [TF, steps] = checkEntropy(self, entropyTol)
            
            global INDEX_UP INDEX_DOWN
            
            TF = true;
            steps = struct;
            
            % check each junction
            for junc = self.net.junc_labels'
                
                juncStr = sprintf('junc_%d',junc);
                
                % initialize as empty
                steps.(juncStr) = [];
                
                % the boudnary discretization grid at this junction;
                % they are also saved in the links
                T_grid = self.net.network_junc.(juncStr).T;
                T_cum_grid = self.net.network_junc.(juncStr).T_cum;
                
                % Here consider connection case
                if strcmp(self.net.network_junc.(juncStr).type_junc,'connection')
                    
                    outlink = self.net.network_junc.(juncStr).outlabel;
                    linkStr = sprintf('link_%d', outlink);
                    
                    % extract the thru flow value which is downstream
                    % inflow in the connection case
                    q_thru = self.x(self.dv_index.(linkStr)(1,INDEX_UP):...
                                   self.dv_index.(linkStr)(2,INDEX_UP));
                    
                    % check each step
                    % break once found the first non-entropic step
                    for i = 1: length(T_grid)
                        
                        % set the reference point of this step
                        self.t_ref = T_cum_grid(i);   
                        t_start = T_cum_grid(i);
                        t_end = T_cum_grid(i+1); % end time of this point
                        
                        % Sample the start and end point of this step;
                        % compute the min of sending and receiving function
                        % which should be entropic solution
                        d_M = self.samplePointsJunc(t_end, junc);
                        
                        if abs(q_thru(i)*T_grid(i)-d_M) > entropyTol
                            % difference greater than the tolerance
                            % This may due to the discretization or we
                            % simply did not add entropic component in the
                            % objective function
                            d_M_SR = [self.samplePointsLink(t_end, inlink, 'downstream');
                                      self.samplePointsLink(t_end, outlink, 'upstream')];
                            t_C = (t_start+t_end)/2;
                            d_M_C = [self.samplePointsLink(t_C, inlink, 'downstream');
                                     self.samplePointsLink(t_C, outlink, 'upstream')];
                            
                            if self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(1), d_M_SR(1)]') && ...
                               self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(2), d_M_SR(2)]')     
                                % meaning caused by not setting entropic
                                % solution or traffic control
                                warning('WARNING: Step %d is not entropic\n', i);
                                continue
                            else
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the entropic solution. In this case, the
                            % affected domain of the discretization error
                            % is small. we treat this as entropic
                            continue
                        end
                    end     %end each step
                    
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'merge')
                    
                    % extract the outflow of the two incoming links
                    inlinks = self.net.network_junc.(juncStr).inlabel;
                    outlink = self.net.network_junc.(juncStr).outlabel;
                    linkStr = sprintf('link_%d',inlinks(1));
                    q_s1 = self.x(self.dv_index.(linkStr)(1,INDEX_DOWN):...
                                 self.dv_index.(linkStr)(2,INDEX_DOWN) );
                    linkStr = sprintf('link_%d',inlinks(2));
                    q_s2 = self.x(self.dv_index.(linkStr)(1,INDEX_DOWN):...
                                 self.dv_index.(linkStr)(2,INDEX_DOWN) );
                    
                    % check each step
                    for i = 1: length(T_grid)
                        
                        self.t_ref = T_cum_grid(i);   % start time of this point
                        t_start = T_cum_grid(i);
                        t_end = T_cum_grid(i+1); % end time of this point
                       
                        % a 1 x 2 vector true unique solution from the two incoming links
                        d_M = self.samplePointsJunc( t_end, junc);
                        
                        if abs( q_s1(i)*T_grid(i)-d_M(1) ) > entropyTol ||...
                                abs( q_s2(i)*T_grid(i)-d_M(2) ) > entropyTol
                            % does not match the unique solution from the
                            % sampling approach
                            % the sending and receiving functions on the
                            % links, just to check there is not
                            % intersection on any three links
                            d_M_SR = [self.samplePointsLink(t_end, inlinks(1), 'downstream');
                                      self.samplePointsLink(t_end, inlinks(2), 'downstream');
                                      self.samplePointsLink(t_end, outlink, 'upstream')];
                            t_C = (t_start+t_end)/2;
                            d_M_C = [self.samplePointsLink(t_C, inlinks(1), 'downstream');
                                     self.samplePointsLink(t_C, inlinks(2), 'downstream');
                                     self.samplePointsLink(t_C, outlink, 'upstream')];
                            
                            if self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(1), d_M_SR(1)]') &&...
                                  self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(2), d_M_SR(2)]') &&...
                                  self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(3), d_M_SR(3)]')
                                % meaning caused by not setting entropic
                                % solution or traffic control
                                warning('WARNING: Step %d is not entropic\n', i);
                                continue
                            else
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the entropic solution. In this case, the
                            % affected domain of the discretization error
                            % is small. we treat this as entropic
                            continue
                        end
                    end     % end each step
                    
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'diverge')
                                        
                    % extract the outflow of the two incoming links
                    outlinks = self.net.network_junc.(juncStr).outlabel;
                    inlink = self.net.network_junc.(juncStr).inlabel;
                    linkStr = sprintf('link_%d',outlinks(1));
                    q_r1 = self.x(self.dv_index.(linkStr)(1,INDEX_UP):...
                                 self.dv_index.(linkStr)(2,INDEX_UP) );
                    linkStr = sprintf('link_%d',outlinks(2));
                    q_r2 = self.x(self.dv_index.(linkStr)(1,INDEX_UP):...
                                 self.dv_index.(linkStr)(2,INDEX_UP) );
                    
                    % check each step
                    for i = 1: length(T_grid)
                        
                        self.t_ref = T_cum_grid(i);   % start time of this point
                        t_start = T_cum_grid(i);
                        t_end = T_cum_grid(i+1); % end time of this point
                        
                        % a 1 x 2 vector true unique solution from the two incoming links
                        d_M = self.samplePointsJunc( t_end, junc);
                        
                        if abs( q_r1(i)*T_grid(i)-d_M(1) ) > entropyTol ||...
                                abs( q_r2(i)*T_grid(i)-d_M(2) ) > entropyTol
                            % does not match the unique solution from the
                            % sampling approach
                            d_M_SR = [self.samplePointsLink(t_end, outlinks(1), 'upstream');
                                      self.samplePointsLink(t_end, outlinks(2), 'upstream');
                                      self.samplePointsLink(t_end, inlink, 'downstream')];
                            t_C = (t_start+t_end)/2;
                            d_M_C = [self.samplePointsLink(t_C, outlinks(1), 'upstream');
                                     self.samplePointsLink(t_C, outlinks(2), 'upstream');
                                     self.samplePointsLink(t_C, inlink, 'downstream');];
                            
                            if self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(1), d_M_SR(1)]') && ...
                               self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(2), d_M_SR(2)]') && ...
                               self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(3), d_M_SR(3)]')
                                % meaning caused by not setting entropic
                                % solution or traffic control
                                warning('WARNING: Step %d is not entropic\n', i);
                                continue
                            else
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the entropic solution. In this case, the
                            % affected domain of the discretization error
                            % is small. we treat this as entropic
                            continue
                        end
                    end     % end each step
                    
                end     % end diverge
                
                
            end % end each junction
            
            % steps = unique(steps);
            
        end % end function
        
        
        %===============================================================
        % analytically update the discretization by identify the 
        % intersection point of the shock waves at the boundaries. Those
        % points will be added to the discretization grid to eliminate the
        % discretization error.
        % input: 
        %       steps: struct, .(juncStr), the first non-entropic step. 
        %           Note, the update has to be done iteratively, 
        %           earlier steps first
        % output: 
        %       T_new_cum: struct, .(juncStr). The updated time grid. 
        function T_new_cum = updateTimeDiscretization(self, steps)
            
            global INDEX_UP INDEX_DOWN 
            
            % recursively search uses a binary cut method to locate the
            % intersection point. 
            % Set the following parameters to limit the number of
            % iterations: stop if we think the position is roughtly good.
            searchDepth = 0;
            d_t = 1.0e-3;    % 0.001 second
            
            % the new discretization grid to be returned.
            T_new_cum = struct;
            
            % check each junction
            for junc = self.net.junc_labels'
                
                juncStr = sprintf('junc_%d',junc);
                
                % get the link label at this junction
                inlink = self.net.network_junc.(juncStr).inlabel;
                outlink = self.net.network_junc.(juncStr).outlabel;
                
                T_tmp_cum = self.net.network_junc.(juncStr).T_cum;
                end_time = self.net.network_junc.(juncStr).T_cum(end);
                num_steps = length(self.net.network_junc.(juncStr).T);
                
                % find the start and end time of the next step of the
                % non-entropic step, since the non-entropy could be caused 
                % by the intersection from the following step.
                if ~isempty(steps.(juncStr))
                    
                    nonentropic_step = steps.(juncStr);
                    
                    self.t_ref = self.net.network_junc.(juncStr).T_cum(nonentropic_step);
                    t_left = self.net.network_junc.(juncStr).T_cum(nonentropic_step);
                
                    if nonentropic_step < num_steps
                        t_right = self.net.network_junc.(juncStr).T_cum(nonentropic_step + 2);
                    else
                        t_right = self.net.network_junc.(juncStr).T_cum(nonentropic_step+1);
                    end
                end
                    
                
                if strcmp( self.net.network_junc.(juncStr).type_junc, 'connection')
                    
                    if ~isempty(steps.(juncStr))
                        
                        % search the intersection of shockwaves on incoming links
                        % t_found = searchIntersection(obj, t_interval, M_interval, slope, ...
                        %                      link, bound, searchDepth, dt_tol)
                        t_found_s = searchIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink, 'downstream', searchDepth, d_t);
                        % search the intersection of shockwaves on outgoing links
                        t_found_r = searchIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    % This reduces the number of discretized intervals in
                    % the boudnary without causing discretization error
                    linkStr = sprintf('link_%d', outlink);
                    q_thru = self.x(self.dv_index.(linkStr)(1, INDEX_UP):...
                                   self.dv_index.(linkStr)(2, INDEX_UP));
                    tmp_group = groupSameElement(q_thru, inf);
                    T_tmp_cum = [T_tmp_cum(tmp_group(:,1));...
                                 end_time];
                    T_tmp_cum = unique(T_tmp_cum);
                                        
                elseif strcmp( self.net.network_junc.(juncStr).type_junc, 'merge')
                    
                    if ~isempty(steps.(juncStr))
                    
                        % corner points found from sending link 1, 2 and
                        % receiving link 3
                        t_found_s1 = searchIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink(1), 'downstream', searchDepth, d_t);
                       
                        t_found_s2 = searchIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink(2), 'downstream', searchDepth, d_t);
                                        
                        t_found_r = searchIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s1; t_found_s2; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', inlink(1));
                    q_s1 = self.x( self.dv_index.(linkStr)(1, INDEX_DOWN):...
                                  self.dv_index.(linkStr)(2, INDEX_DOWN) );
                    linkStr = sprintf('link_%d', inlink(2));
                    q_s2 = self.x( self.dv_index.(linkStr)(1, INDEX_DOWN):...
                                  self.dv_index.(linkStr)(2, INDEX_DOWN) );
                    
                    T_tmp_cum_s1 = T_tmp_cum;
                    tmp_group = groupSameElement(q_s1, inf);
                    T_tmp_cum_s1 = [ T_tmp_cum_s1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_s2 = T_tmp_cum;
                    tmp_group = groupSameElement(q_s2, inf);
                    T_tmp_cum_s2 = [ T_tmp_cum_s2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_s1; T_tmp_cum_s2]);
                                                  
                elseif strcmp( self.net.network_junc.(juncStr).type_junc, 'diverge')
                     
                    if ~isempty(steps.(juncStr))
                        
                        % corner points found from sending link 1 and
                        % receiving link 2, 3
                        t_found_s = searchIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink, 'downstream', searchDepth, d_t);
                        t_found_r1 = searchIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink(1), 'upstream', searchDepth, d_t);
                        t_found_r2 = searchIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink(1), 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r1; t_found_r2];
                    
                    else 
                        t_found = [];
                    end
                  
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', outlink(1));
                    q_r1 = self.x( self.dv_index.(linkStr)(1, INDEX_UP):...
                                  self.dv_index.(linkStr)(2, INDEX_UP) );
                    linkStr = sprintf('link_%d', outlink(2));
                    q_r2 = self.x( self.dv_index.(linkStr)(1, INDEX_UP):...
                                  self.dv_index.(linkStr)(2, INDEX_UP) );
                    
                    T_tmp_cum_r1 = T_tmp_cum;
                    tmp_group = groupSameElement(q_r1, inf);
                    T_tmp_cum_r1 = [ T_tmp_cum_r1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_r2 = T_tmp_cum;
                    tmp_group = groupSameElement(q_r2, inf);
                    T_tmp_cum_r2 = [ T_tmp_cum_r2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_r1; T_tmp_cum_r2]);
                
                else
                    t_found = [];
                end
                
                % add new found intersection point in the discretization grid
                T_tmp_cum = [T_tmp_cum; t_found(t_found<end_time)];
                
                % return sorted and unique discretization accumulative points
                T_new_cum.(juncStr) = unique_tol(T_tmp_cum, 1.0e-8);
                
            end
            
        end
        
        
        
        %===============================================================
        % This is a recursive function aiming at finding the corner point
        % of a piecewise linear nondecreasing function in an interval
        % In principle, you only need the time interval and the function f
        % for computing the values; the rest of the input parameters are
        % just to avoid recomputing some values and save time.
        % Note: for computational efficiency, this code may not work for a
        % piece wise line with too many corners and certain properties. But
        % this is very unlikely to happen in one sigle time step.
        % input: 
        %        t_interval, float vector, [t_start; t_end]
        %        M_interval, float vector, [M_start; M_end], NaN if not
        %           computed
        %        link, int, the link id which defines f
        %        bound, 'upstream', 'downstream', which defines f
        %        searchDepth, int, keeps track of the search depth
        %        dr_tol, float, the minimal resolution we will investigate
        % output: found time points 1 x n column
        function t_found = searchIntersection(self, t_interval, M_interval, slope, ...
                                              link, bound, searchDepth, dt_tol)
            
            t_found = [];                  
            % if search depth is greater than 3, 
            % or the time interval is shorter than dt_res, return []       
            if searchDepth > 3 || t_interval(2)-t_interval(1) <= dt_tol
                return
            end
            
            % First, compute the two boundary values and the center point
            t_L = t_interval(1);
            t_R = t_interval(2);
            t_C = (t_L+t_R)/2;
            % If value not computed yet
            if isnan(M_interval(1))
                M_L = self.samplePointsLink(t_L, link, bound);
            else
                M_L = M_interval(1);
            end
            if isnan(M_interval(2))
                M_R = self.samplePointsLink(t_R, link, bound);
            else
                M_R = M_interval(2);
            end
            M_C = self.samplePointsLink(t_C, link, bound);
            
            % Second, check if on a straight line, return []
            if self.onStraightLine([t_L; t_C; t_R], [M_L; M_C; M_R])
                % double check if on a straight line. Some times, the
                % segments are symmetric and those 3 points are on the same
                % line, but the sending function is not a line
                t_Q_1 = t_L + (t_R-t_L)/4;      % first quarter
                t_Q_3 = t_L + 3*(t_R-t_L)/4;    % third quarter
                
                M_Q_1 = self.samplePointsLink(t_Q_1, link, bound);
                M_Q_3 = self.samplePointsLink(t_Q_3, link, bound);
                
                if self.onStraightLine([t_L; t_Q_1; t_C], [M_L; M_Q_1; M_C]) &&...
                   self.onStraightLine([t_C; t_Q_3; t_R], [M_C; M_Q_3; M_R])
                    % sampled 5 points on the same line, then very likely
                    % to be a straight line
                    return
                end
            end
            
            % Third, find the left and right slope and compute the
            % intersection 
            if isnan(slope(1))
                s_L = self.findSlope([t_L; t_C], [M_L; M_C], 'left',...
                                 link, bound, dt_tol);
            else
                s_L = slope(1);
            end
            if isnan(slope(2))
                s_R = self.findSlope([t_C; t_R], [M_C; M_R], 'right',...
                            	 link, bound, dt_tol);
            else
                s_R = slope(2);
            end
            % intersection point 
            if abs(s_L - s_R) >= 1.0e-3
                % make sure two slopes are not parallel, otherwise directly
                % split them.
                t_insct = (M_R - M_L + s_L*t_L - s_R*t_R)/(s_L - s_R);
                M_insct = s_L*(t_insct - t_L) + M_L;
                
                if t_insct <= t_R && t_insct >= t_L
                    % if the intersection is in the time step, then
                    % validate, otherwise, move on to split
                    M_validate = self.samplePointsLink(t_insct, link, bound);
                    if abs(M_insct - M_validate) <= 1e-6
                        % found the intersection
                        t_found = [t_found; t_insct];
                        return
                    end
                end
            end
            
            % Forth, recursion, split the interval if there are more than 
            % one corners based on: 
            % 1. slopes are parallel; 2. intersection outside of
            % the interval; 3. cornerpoint not valid in sample
            t_found_left = self.searchIntersection([t_L; t_C], [M_L; M_C],...
                [s_L; NaN], link, bound, searchDepth+1, dt_tol);
            t_found_right = self.searchIntersection([t_C; t_R], [M_C; M_R],...
                [NaN; s_R], link, bound, searchDepth+1, dt_tol);
            
            t_found = [t_found; t_found_left; t_found_right];
                        
        end
        
        
        %===============================================================
        % This is a recursive function that finds the slope at the boundary
        % point of a piecewise linear line. A utility function for finding
        % the intersections
        % input: 
        %        t_interval, float vector, [t_start; t_end]
        %        M_interval, float vector, [M_start; M_end], NaN if not
        %           computed
        %        slope_side, 'left','right', which side slope
        %        link, int, the link id which defines f
        %        bound, 'upstream', 'downstream', which defines f
        %        dr_tol, float, the minimal resolution we will investigate
        % output: found time points 1 x n column
        function slope = findSlope(self, t_interval, M_interval, slope_side, ...
                                   link, bound, dt_tol)
                              
            % if the time interval is too small, than stop computing slope
            % since this may give large numerical error
            if t_interval(2) - t_interval(1) <= dt_tol
                slope = NaN;
                return
            end
            
            % First, sample boundary points values and the center point
            t_L = t_interval(1);
            t_R = t_interval(2);
            t_C = (t_L+t_R)/2;
            % If value not computed yet
            if isnan(M_interval(1))
                M_L = self.samplePointsLink(t_L, link, bound);
            else
                M_L = M_interval(1);
            end
            if isnan(M_interval(2))
                M_R = self.samplePointsLink(t_R, link, bound);
            else
                M_R = M_interval(2);
            end
            
            M_C = self.samplePointsLink(t_C, link, bound);
            
            % Second, check if three points are on the same line
            if self.onStraightLine([t_L; t_C; t_R],[M_L; M_C; M_R])
                % if true, return the slope
                slope = (M_R-M_L)/(t_R-t_L);
            else
                if strcmp(slope_side, 'left')
                    % keep the left half to find the slope
                    slope = self.findSlope([t_L; t_C], [M_L; M_C], slope_side,...
                                         link, bound, dt_tol);
                elseif strcmp(slope_side, 'right')
                    % keep the right half to find the slope
                    slope = self.findSlope([t_C; t_R], [M_C; M_R], slope_side,...
                                         link, bound, dt_tol);
                end
            end
                              
                              
        end
        
        
        
        %===============================================================
        % This function tells whether three points are on the same straight
        % line
        % input:
        %       x: 3x1 vector, float
        %       y: 3x1 vector, float
        % output:
        %       TF: true if three points on one straight line with small
        %           tolerance; otherwise false
        function TF = onStraightLine(self, x, y)
            
            % if a horizontal line
            if abs(x(3)-x(2)) <= 1.0e-10 && abs(x(2)-x(1)) <= 1.0e-6
                TF = true;
                return
            end
            
            % if a vertical line
            if abs(y(3)-y(2)) <= 1.0e-10 && abs(y(2)-y(1)) <= 1.0e-6
                TF = true;
                return
            end
            
            % otherwise, compute slope
            slope = (y(3)-y(1))/(x(3)-x(1));
            if abs(slope*(x(2)-x(1)) + y(1) - y(2)) <= 1.0e-6
                TF = true;
            else
                TF = false;
            end
            
        end
        
        
        %===============================================================
        % sample points at boundaries at a junction
        % This function basically samples the specified point at all links,
        % and then compute the entropic value based on the sending and
        % receiving function.
        % Note the returned value is the vehicle id with reference to the
        % time point self.t_ref
        % input: 
        %        t_sample: float, the time point to be sampled
        %        junc: the label of the junction to be sampled
        % output: 
        %       d_M: for connection, float, the entropic car ID at t_sample
        %            for merge/diverge, 2x1 float, the entropic car ID at
        %            two incoming/outgoing links
        function d_M = samplePointsJunc(self, t_sample, junc)
            
            % t_ref is the reference point where M = 0
            if t_sample == self.t_ref
                d_M = 0;
                return
            end
            
            juncStr = sprintf('junc_%d',junc);
            
            if strcmp(self.net.network_junc.(juncStr).type_junc,'connection')
                
                inLink = self.net.network_junc.(juncStr).inlabel;
                outLink = self.net.network_junc.(juncStr).outlabel;
                
                % simply take the minimum of the sending and receiving of
                % two links
                
                d_M = min(self.samplePointsLink(t_sample, inLink,'downstream'),...
                          self.samplePointsLink(t_sample, outLink,'upstream'));
                
            elseif strcmp(self.net.network_junc.(juncStr).type_junc,'merge')
                
                inLinks = self.net.network_junc.(juncStr).inlabel;
                outLink = self.net.network_junc.(juncStr).outlabel;
                
                % slope of the priority with y: inLinks(2); and x:
                % inlinks(1)
                sPriority = self.net.network_junc.(juncStr).ratio(2)/...
                    self.net.network_junc.(juncStr).ratio(1);
                
                % here we extract the sending and receiving function
                M_R = self.samplePointsLink(t_sample, outLink, 'upstream');
                M_S1 = self.samplePointsLink(t_sample, inLinks(1), 'downstream');
                M_S2 = self.samplePointsLink(t_sample, inLinks(2), 'downstream');
                
                % case 1: M_S1 + M_S2 <= M_R, then through flow is
                % M_s1 + M_s2
                if M_S1 + M_S2 <= M_R
                    d_M = [M_S1; M_S2];
                
                % case 2: intersection point in side (0-M_S1, 0-M_S2) box
                elseif ( M_R/(1+sPriority) <= M_S1) &&...
                        (sPriority*M_R/(1+sPriority) <= M_S2)
                    
                    d_M = [M_R/(1+sPriority); sPriority*M_R/(1+sPriority)];                    
                
                % case 3: constrained by M_S1 flow
                elseif M_R/(1+sPriority) > M_S1
                    d_M = [M_S1, M_R-M_S1];
                
                % case 4: constrained by M_S2 flow
                elseif sPriority*M_R/(1+sPriority) > M_S2
                    d_M = [M_R-M_S2 , M_S2];
                end

            elseif strcmp(self.net.network_junc.(juncStr).type_junc,'diverge')
                
                inLink = self.net.network_junc.(juncStr).inlabel;
                outLinks = self.net.network_junc.(juncStr).outlabel;
                
                % slope of the distribution with y: outLinks(2); and x:
                % outlinks(1)
                sDistribution = self.net.network_junc.(juncStr).ratio(2)/...
                    self.net.network_junc.(juncStr).ratio(1);
                
                % here we extract the sending and receiving function
                M_R1 = self.samplePointsLink(t_sample, outLinks(1), 'upstream');
                M_R2 = self.samplePointsLink(t_sample, outLinks(2), 'upstream');
                M_S = self.samplePointsLink(t_samplw, inLink, 'downstream');
                
                % case 1: M_R1 + M_R2 <= M_S, then through flow is
                % M_R1+M_R2
                if M_R1 + M_R2 <= M_S
                    d_M = [M_R1, M_R2];
                    
                % case 2: intersection point in side (0-M_R1, 0-M_R2) box
                elseif ( M_S/(1+sDistribution) <= M_R1) &&...
                        (sDistribution*M_S/(1+sDistribution) <= M_R2)
                    d_M = [M_S/(1+sDistribution), sDistribution*M_S/(1+sDistribution)];
                    
                % case 3: constrained by M_R1 flow
                elseif M_S/(1+sDistribution) > M_R1
                    d_M = [M_R1, M_S-M_R1];
                    
                % case 4: constrained by M_R2 flow
                elseif sDistribution*M_S/(1+sDistribution) > M_R2
                    d_M = [M_S-M_R2 , M_R2];
                end

            end
            
        end
        
        
        %===============================================================
        % sample the boundary point of a link with respect to the reference
        % time self.t_ref at its corresponding bound.
        % input: 
        %        t_sample: float, the boundary time to be sampled
        %        link: int, the link ID to be sampled
        %        bound: string, 'upstream', 'downstream', specify which
        %               boundary we would like to sample
        % output: M, ( a vector of vehicle IDs ).
        function M = samplePointsLink(self, t_sample, link, bound)
            global INDEX_UP INDEX_DOWN
            
            % since t_ref is the reference point, so simply 0
            if t_sample == self.t_ref
                M = 0;
                return
            end
            
            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            M = ones(2,1)*NaN;
            
            % compute the M of t on this link
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.net.network_hwy.(linkStr).para_vf;
            w = self.net.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.net.network_hwy.(linkStr).para_kc;
            k_m = self.net.network_hwy.(linkStr).para_km;            
            
            % Get the initial number of vehicles 
            IC(:,1) = self.net.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.net.network_hwy.(linkStr).X_grid_cum(2:end);
            IC(:,3) = self.net.network_hwy.(linkStr).IC;
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end); %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end
            
            % select the boudnary conditions we need to use to compute the
            % vehicle id at the sampling point.
            % We only need to use the boundary conditions both upstream and
            % downstream up to the sampling point time.
            % the number of boundary steps we need to take into account
            num_us_pre_steps = sum(self.net.network_hwy.(linkStr).T_us_cum < t_sample);
            num_ds_pre_steps = sum(self.net.network_hwy.(linkStr).T_ds_cum < t_sample);
            
            
            % First extract the boundary conditions that to be considered
            % for computing the vehilce label at the sampling point
            BC_us(:,1) = self.net.network_hwy.(linkStr).T_us_cum(1:num_us_pre_steps);
            BC_us(:,2) = self.net.network_hwy.(linkStr).T_us_cum(2:num_us_pre_steps+1);
            BC_us(num_us_pre_steps, 2) = t_sample;
            BC_ds(:,1) = self.net.network_hwy.(linkStr).T_ds_cum(1:num_ds_pre_steps);
            BC_ds(:,2) = self.net.network_hwy.(linkStr).T_ds_cum(2:num_ds_pre_steps+1);
            BC_ds(num_ds_pre_steps, 2) = t_sample;
            
            % Use the computed solution values as the boundary condition
            BC_us(:,3) = self.x( self.dv_index.(linkStr)(1, INDEX_UP):...
                                self.dv_index.(linkStr)(1, INDEX_UP) + num_us_pre_steps-1);
            BC_ds(:,3) = self.x( self.dv_index.(linkStr)(1, INDEX_DOWN):...
                                self.dv_index.(linkStr)(1, INDEX_DOWN) + num_ds_pre_steps-1);
            
            % Second compute the absolute vehicle ID
            BC_us_num_veh = (BC_us(:,2) - BC_us(:,1)).*BC_us(:,3);
            BC_cum_us_M = [0; cumsum(BC_us_num_veh)];
            BC_cum_us_M(end) = [];
            
            BC_ds_num_veh = (BC_ds(:,2) - BC_ds(:,1)).*BC_ds(:,3);
            BC_cum_ds_M = [0; cumsum(BC_ds_num_veh) ] + IC_total_num_veh;
            BC_cum_ds_M(end) = [];
            
            % if the sample is on the upstream bound
            if strcmp(bound,'upstream')
                
                if sum( self.net.network_hwy.(linkStr).T_us_cum > self.t_ref &...
                        self.net.network_hwy.(linkStr).T_us_cum < t_sample) ~= 0
                    % remove the last two pieces of upstream conditions
                    BC_us( max(1, num_us_pre_steps-1) : num_us_pre_steps, 3) = NaN; 
                else
                    % remove the last piece of upstream condition
                    BC_us( num_us_pre_steps, 3) = NaN;
                end
                
            % if the sample is on the downstream bound
            elseif strcmp(bound,'downstream')
                
                if sum( self.net.network_hwy.(linkStr).T_ds_cum > self.t_ref &...
                        self.net.network_hwy.(linkStr).T_ds_cum < t_sample ) ~= 0
                    % remove the last two pieces of downstream conditions
                    BC_ds( max(1,num_ds_pre_steps-1):num_ds_pre_steps , 3) = NaN;
                else
                    % remove the last piece of downstream condition
                    BC_ds( num_ds_pre_steps , 3) = NaN;
                end
            else
                error('ERROR: You have to specify which bound you would like to sample');
            end
            
            % points coordinates
            len_link = self.net.network_hwy.(linkStr).para_postkm*1000;
            if strcmp(bound,'upstream')
                position = 0;
            elseif strcmp(bound,'downstream')
                position = len_link;
            end
            
            t_array = [self.t_ref; t_sample];
            for i = 1:length(t_array)
                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*t_array(i) <= position &...
                                    IC(:,2)+v_f*t_array(i) >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*t_array(i) > position &...
                              IC(:,1)+w*t_array(i) <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*t_array(i) - ...
                    IC(activeCharacterDomain,1)) );
                M(i) = minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*t_array(i) - ...
                    IC(activeFanDomain,1)  )  );
                M(i) = minNonEmpty(M(i), tmp_M);
                
                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*t_array(i) <= position &...
                                    IC(:,2)+w*t_array(i) >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*t_array(i) > position &...
                              IC(:,2)+w*t_array(i) <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*t_array(i) - ...
                    IC(activeCharacterDomain,1)) -k_m*t_array(i)*w);
                M(i) = minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*t_array(i) - IC(activeFanDomain,2)) ...
                    - k_m*t_array(i)*w);
                M(i) = minNonEmpty(M(i), tmp_M);
                
                %UPSTREAM==================================================
                % now check the upstream conditions, only apply those
                % boundary conditions that are not NaN
                notNaN = ~isnan(BC_us(:,3));
                
                inCharacterDomain = BC_us(:,1)+ position/v_f <= t_array(i) &...
                                    BC_us(:,2)+ position/v_f >= t_array(i); 
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_us(:,2)+ position/v_f < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_us_M(activeCharacterDomain,1) +...
                    BC_us(activeCharacterDomain,3).*(t_array(i) - position/v_f - ...
                    BC_us(activeCharacterDomain,1)) );
                M(i) = minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_us_M(activeFanDomain,1) +...
                    BC_us_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - position/v_f - ...
                    BC_us(activeFanDomain,2)) );
                M(i) = minNonEmpty(M(i), tmp_M);
                
                
                %DOWNSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(BC_ds(:,3));
                
                inCharacterDomain = BC_ds(:,1) + (position - len_link)/w <= t_array(i) &...
                                    BC_ds(:,2)+ (position - len_link)/w >= t_array(i);
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_ds(:,2)+ (position - len_link)/w < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_ds_M(activeCharacterDomain,1) +...
                    BC_ds(activeCharacterDomain,3).*(t_array(i) -...
                    (position - len_link)/w - ...
                    BC_ds(activeCharacterDomain,1))...
                    - k_m*(position - len_link));
                M(i) = minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_ds_M(activeFanDomain,1) +...
                    BC_ds_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - (position-len_link)/v_f - ...
                    BC_ds(activeFanDomain,2)) );
                
                M(i) = minNonEmpty(M(i), tmp_M);
                
                
            end % end for each t (here we just have (t_ref and t) )
            
            if all(~isnan(M)) && all(~isempty(M))
                
                % use car ID at t_ref as a reference
                M = M(2) - M(1);
                
            else
                error('failed to sample points');
            end
            
            
        end 
        
        
    end     %end methods
    
end






















