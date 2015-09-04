% Yanning Li, Sep 01, 2015
% This class handles the construction of a network object.
% The network object is defiend as a struct with methods to add highway 
% links and junctions.
% network:
% ---- hwy: a struct, containing information for each hwy link
% ---- junc: a struct, containing the informatino for each junc, including
%            the time grid at the junction
% ---- initial and boundary conditions
% ---- other auxiliary properties

classdef initNetwork < handle
    
    properties
        
        % defined as struct, containing information associated with highway
        % links and junctions
        % Initial and boundary conditions are saved in the links
        network_junc;
        network_hwy;
        
        % Auxiliary properties
        % a vector array containing the labels of juncs and links
        junc_labels;
        link_labels;
        num_juncs;
        num_links;
        
    end
    
    methods
        %===============================================================
        function self = initNetwork(~)
                        
            self.network_junc = struct;
            self.network_hwy = struct;
            
            self.junc_labels = [];
            self.link_labels = [];
            self.num_links = 0;
            self.num_juncs = 0;
            
        end
        
        %===============================================================
        % Create a junction. Only support connection, merge, and diverge
        % input: 
        %       junc: int, assign a junction id
        %       inlabel/outlabel: numerical vectors containing the labels of links
        %       type_junc: strings 'connection', 'merge' or 'diverge'
        %       ratio: 2x1 numerical array in same order to inlabel/outlabel.
        %           e.g. [0.4; 0.6]
        %       T: num_steps x 1 numerical array; time durations for each time step
        function addJunc(self, junc, inlabel, outlabel, type_junc, ratio, T)
            
            % make sure no duplicated label
            if any(self.junc_labels == junc)
                error('ERROR: Duplicated junc label: %d.\n', junc)
            end
            
            juncStr = sprintf('junc_%d', junc);
            self.network_junc.(juncStr).inlabel = inlabel;
            self.network_junc.(juncStr).outlabel = outlabel;
            self.network_junc.(juncStr).type_junc = type_junc;
            self.network_junc.(juncStr).ratio = ratio;
            
            % set the time duration and the time grid
            self.network_junc.(juncStr).T = T;
            self.network_junc.(juncStr).T_cum = [0; cumsum(T)];
            
            self.junc_labels = [self.junc_labels; junc];
            self.num_juncs = self.num_juncs + 1;
        end
        
        %===============================================================
        % Add a link to the network.
        % Link is refered as a segment on the highway. The highway should
        % be divided to links at locations where sensors are deployed, or
        % the property of the road chagnes (drops of lanes).
        % input:
        %       link: numerical label
        %       para: a struct containing road properties
        %       num_lanes: int, the number of lanes
        %       lengthKM: float, length of the link in km
        %       linkType: string, 'freeway', 'onramp', 'offramp'
        function addLink(self, link, para, num_lanes, lengthKM, linkType)
            
            % make sure no duplicated label
            if any(self.link_labels == link)
                error('ERROR: Duplicated link label.\n')
            end
            
            linkStr = sprintf('link_%d',link);
            
            % set the road parameters
            self.network_hwy.(linkStr).para_vf = para.v;
            self.network_hwy.(linkStr).para_w = para.w;
            self.network_hwy.(linkStr).para_kc = num_lanes*para.k_c_pl; %kc per lane
            self.network_hwy.(linkStr).para_km = num_lanes*para.k_m_pl;
            self.network_hwy.(linkStr).para_qmax = num_lanes*para.q_max_pl;
            self.network_hwy.(linkStr).para_vmin = para.v_min;
            self.network_hwy.(linkStr).para_vmax = para.v_max;
            
            self.network_hwy.(linkStr).para_postkm = lengthKM;
            self.network_hwy.(linkStr).para_linktype = linkType;
            
            % save the label
            self.link_labels = [self.link_labels; link];
            self.num_links = self.num_links + 1;
            
        end
        
        %===============================================================
        % Set initial condition of each link, assuming even discretization 
        % over space
        % input: 
        %       link: int, link label
        %       rho_ini: num_segments x 1 float, normalized to rho_c
        function setInitialCon(self, link, rho_ini)
            
            if ~any(self.link_labels == link)
                error('ERROR: Link %d not defined in the network\n', link);
            end
            
            num_seg = length(rho_ini); 
            linkStr = sprintf('link_%d',link);
            
            % Each link saves its own initial condition as absolute values,
            % as well as the space grid, 
            % length(grid vector) = length(rho_ini)+1
            tmp_dx = self.network_hwy.(linkStr).para_postkm*1000/num_seg;
            
            self.network_hwy.(linkStr).X_grid_cum = (0:num_seg)'*tmp_dx;
            self.network_hwy.(linkStr).IC = self.columnize(rho_ini)*...
                self.network_hwy.(linkStr).para_kc;

        end
        

        %===============================================================
        % Set boundary condition of each link
        % input:
        %       link: int, the link label
        %       q_in: num_in_steps x 1 float, inflow normalized to q_max
        %       q_out: num_out_steps x 1 float, outflow normalized to q_max
        %       T_in: num_in_steps x 1 float durations of each step
        %       T_out: num_out_steps x 1 float durations of each step
        % If q_in or q_out is [], then it is decision variable to be
        % estimated. It will be set as NaN. 
        % Note T_in T_out must be set. One of it may be the standard 
        % discretization at the boundary of the network. The other one may 
        % be the dynamically set. 
        function setBoundaryCon(self, link, q_in, q_out, T_in, T_out)
                        
            if ~any(self.link_labels == link)
                error('ERROR: Link %d not defined in the network.\n', link);
            end
            
            % the flow data and time must have the same length
            if (~isempty(q_in) && length(q_in)~=length(T_in)) ||...
                    (~isempty(q_out) && length(q_out)~=length(T_out))
                error('Data length must be equal to its time discretization');
            end
            
            num_steps_in = length(T_in);
            num_steps_out = length(T_out);
            
            linkStr = sprintf('link_%d',link);
            
            % save the boundary discritization grid in struct
            self.network_hwy.(linkStr).T_us = T_in;
            self.network_hwy.(linkStr).T_us_cum = [0; cumsum(T_in)];
            
            self.network_hwy.(linkStr).T_ds = T_out;
            self.network_hwy.(linkStr).T_ds_cum = [0; cumsum(T_out)];
            
            % set the boundary conditions in absolute values veh/hr
            if ~isempty(q_in)
                self.network_hwy.(linkStr).BC_us = self.columnize(q_in)*...
                    self.network_hwy.(linkStr).para_qmax;
            else
                self.network_hwy.(linkStr).BC_us = ones(num_steps_in,1)*NaN;
            end
            
            if ~isempty(q_out)
                self.network_hwy.(linkStr).BC_ds = self.columnize(q_out)*...
                                self.network_hwy.(linkStr).para_qmax;
            else
                self.network_hwy.(linkStr).BC_ds = ones(num_steps_out,1)*NaN;
            end
           
        end
        
        
        % Some useful functions
        function [co] = columnize(self, v)
            if iscolumn(v)
                co = v;
            else
                co = v';
            end
        end
        
     
    end
    
end

    
    






