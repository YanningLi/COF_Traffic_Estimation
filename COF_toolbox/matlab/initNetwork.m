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
        network_junc;
        network_hwy;
        
        % Initial and boundary conditions
        Con_initial;
        q_upstream;
        q_downstream;
        
        % Auxiliary properties
        % a vector array containing the labels of juncs and links
        junc_labels;
        link_labels;
        num_juncs;
        num_links;
        
    end
    
    methods
        %===============================================================
        function obj = initNetwork(~)
            
            global DENS FLUX;
            
            obj.network_junc = struct;
            obj.network_hwy = struct;
            obj.Con_initial = zeros(0,DENS);
            obj.q_upstream = zeros(0,FLUX);
            obj.q_downstream = zeros(0,FLUX);
           
            obj.junc_labels = [];
            obj.link_labels = [];
            obj.num_links = 0;
            obj.num_juncs = 0;
            
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
        function addJunc(obj, junc, inlabel, outlabel, type_junc, ratio, T)
            
            % make sure no duplicated label
            if any(obj.junc_labels == junc)
                error('ERROR: Duplicated junc label: %d.\n', junc)
            end
            
            juncStr = sprintf('junc_%d', junc);
            obj.network_junc.(juncStr).inlabel = inlabel;
            obj.network_junc.(juncStr).outlabel = outlabel;
            obj.network_junc.(juncStr).type_junc = type_junc;
            obj.network_junc.(juncStr).ratio = ratio;
            
            % set the time duration and the time grid
            obj.network_junc.(juncStr).T = T;
            obj.network_junc.(juncStr).T_cum = [0; cumsum(T)];
            
            obj.junc_labels = [obj.junc_labels; junc];
            obj.num_junc = obj.num_junc + 1;
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
        function addLink(obj, link, para, num_lanes, lengthKM, linkType)
            global V_PARA W_PARA K_C_PARA K_M_PARA Q_MAX_PARA V_MIN_PARA V_MAX_PARA;
            global POSTKM LINKTYPE;
            
            % make sure no duplicated label
            if any(obj.link_labels == link)
                error('ERROR: Duplicated link label.\n')
            end
            
            linkStr = sprintf('link_%d',link);
            obj.network_hwy.(linkStr).paras(1,V_PARA) = para.v;
            obj.network_hwy.(linkStr).paras(1,W_PARA) = para.w;
            obj.network_hwy.(linkStr).paras(1,K_C_PARA) = num_lanes*para.k_c_pl;
            obj.network_hwy.(linkStr).paras(1,K_M_PARA) = num_lanes*para.k_m_pl;
            obj.network_hwy.(linkStr).paras(1,Q_MAX_PARA) = num_lanes*para.q_max_pl;
            obj.network_hwy.(linkStr).paras(1,V_MIN_PARA) = para.v_min;
            obj.network_hwy.(linkStr).paras(1,V_MAX_PARA) = para.v_max;
            
            obj.network_hwy.(linkStr).profile(1,POSTKM) = lengthKM;
            obj.network_hwy.(linkStr).profile(1,LINKTYPE) = linkType;
            
            % save the label
            obj.link_labels = [obj.link_labels; link];
            obj.num_links = obj.num_link + 1;
            
        end
        
        %===============================================================
        % Set initial condition of each link, assuming even discretization 
        % over space
        % input: 
        %       link: int, link label
        %       rho_ini: num_segments x 1 float, normalized to rho_c
        function setInitialCon(obj, link, rho_ini)
            global LINKLABEL TIME P_MIN P_MAX DENS K_C_PARA POSTKM
            
            if ~any(obj.link_labels == link)
                error('ERROR: Link %d not defined in the network\n', link);
            end
            
            num_seg = length(rho_ini); 
            linkStr = sprintf('link_%d',link);
            
            tmp_link = zeros(num_seg,DENS);  
            tmp_link(:,LINKLABEL) = link;
            tmp_link(:,TIME)= 0;
            
            tmp_dx = obj.network_hwy.(linkStr).profile(1,POSTKM)*1000/num_seg;
            tmp_link(:,P_MIN) = (0:num_seg-1)'*tmp_dx;
            tmp_link(:,P_MAX) = (1:num_seg)'*tmp_dx;
            
            tmp_link(:,DENS) = obj.columnize(rho_ini)*obj.network_hwy.(linkStr).paras(1,K_C_PARA);
            
            % save in the format defined in the toolbox.
            obj.Con_initial = [obj.Con_initial; tmp_link];
        end
        

        %===============================================================
        % Set boundary condition of each link
        % input:
        %       link: int, the link label
        %       q_in: num_in_steps x 1 float, inflow normalized to q_max
        %       q_out: num_out_steps x 1 float, outflow normalized to q_max
        %       T_in: num_in_steps x 1 float durations of each step
        %       T_out: num_out_steps x 1 float durations of each step
        function setBoundaryCon(obj, link, q_in, q_out, T_in, T_out)
            global Q_MAX_PARA LINKLABEL T_MIN T_MAX FLUX;
            
            if ~any(obj.link_labels == link)
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
            
            q_upstreamTmp(1:num_steps_in,LINKLABEL) = link;
            q_downstreamTmp(1:num_steps_out,LINKLABEL) = link;
            
            % set the boundary conditions t_min and t_max
            T_in_cum = [0; cumsum(T_in)];
            T_out_cum = [0; cumsum(T_out)];
            
            q_upstreamTmp(1:num_steps_in,T_MIN) = T_in_cum(1:end-1);
            q_upstreamTmp(1:num_steps_in,T_MAX) = T_in_cum(2:end);
            q_downstreamTmp(1:num_steps_out,T_MIN) = T_out_cum(1:end-1);
            q_downstreamTmp(1:num_steps_out,T_MAX) = T_out_cum(2:end);
            
            if ~isempty(q_in)
                q_upstreamTmp(1:num_steps_in,FLUX) = ...
                    obj.columnize(q_in)*obj.network_hwy.(linkStr).paras(1,Q_MAX_PARA);
            else
                q_upstreamTmp(1:num_steps_in,FLUX) = NaN;
            end
            
            if ~isempty(q_out)
                q_downstreamTmp(1:num_steps_out,FLUX) = ...
                    obj.columnize(q_out)*obj.network_hwy.(linkStr).paras(1,Q_MAX_PARA);
            else
                q_downstreamTmp(1:num_steps_out,FLUX) = NaN;
            end
            
            obj.q_upstream = [obj.q_upstream; q_upstreamTmp];
            obj.q_downstream = [obj.q_downstream; q_downstreamTmp];
            
        end
        
        
        % Some useful functions
        function [co] = columnize(obj,v)
            if iscolumn(v)
                co = v;
            else
                co = v';
            end
        end
        
     
    end
    
end

    
    






