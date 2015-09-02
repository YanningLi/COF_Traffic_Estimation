% New class for treating solution
% including compute Moskowitz function, plot figures etc
% ::postSolution_dis
% overwrite the methods: checkEntropy, updateDiscretization

classdef postSolution_dis_ana < handle
    
    properties
        x;  % solution
        
        k;  %structure for each link
        N;
        postmScale; %structure for each link
        tScale; %array assuming all links have same time discretization
        fd;
        
        dx_res;
        dt_res
    end
    
    methods
        
        %===============================================================
        % compute the Moskowitz function
        % dx,dt are the resolution for computing Moskowitz
        function obj=postSolution_dis_ana(x,network,program,dx_res,dt_res)
            global NUM_UP NUM_DOWN NUM_INI
            global V_PARA W_PARA K_M_PARA POSTKM
            global TIME LINKLABEL P_MIN P_MAX
            
            obj.x = x;
            
            for link = 1:network.num_link
                linkStr = sprintf('link_%d',link);
                
                q_in =x(program.Dec(link,NUM_UP,1):program.Dec(link,NUM_UP,2),1);
                q_out =x(program.Dec(link,NUM_DOWN,1):program.Dec(link,NUM_DOWN,2),1);
                p_ini =x(program.Dec(link,NUM_INI,1):program.Dec(link,NUM_INI,2),1);
                
                %L = []; % for internal condition, assume none
                %r = [];
                
                obj.fd.(linkStr) = LH_Tfd(network.network_hwy.(linkStr).paras(V_PARA),...
                    network.network_hwy.(linkStr).paras(W_PARA),...
                    network.network_hwy.(linkStr).paras(K_M_PARA));
                xi = 0;
                chi = network.network_hwy.(linkStr).profile(1,POSTKM)*1000;
                pbEnv = LH_general(obj.fd.(linkStr),xi,chi);
                
                %===========================================
                % extract initial segment vector
                Con_ini = network.Con_initial(network.Con_initial(:,TIME)==program.start_time &...
                    network.Con_initial(:,LINKLABEL) == link,:);
                ini_seg = unique([Con_ini(:,P_MIN),Con_ini(:,P_MAX)]);
                
                if size(ini_seg,1)~=1
                    ini_seg = ini_seg';
                end
                if size(p_ini,1)~=1
                    p_ini = p_ini';
                end
                
                pbEnv.setIniDens(ini_seg,p_ini);
                
                a=program.T_cum;
                if ~isrow(a)
                    a = a';
                end
                if ~isrow(q_in)
                    q_in = q_in';
                end
                if ~isrow(q_out)
                    q_out = q_out';
                end
                pbEnv.setUsFlows(a,q_in);
                pbEnv.setDsFlows(a,q_out);
                
                %===========================================
                % specify resolution
                nx = floor((chi)/dx_res);
                dx=(chi)/nx;
                obj.dt_res = dt_res;
                obj.dx_res = dx_res;
                
                obj.postmScale.(linkStr) = 0:dx:chi;
                obj.tScale = 0:dt_res:program.end_time;
                
                xValues = ones(size(obj.tScale'))*(obj.postmScale.(linkStr));
                tValues = obj.tScale' * ones(size(obj.postmScale.(linkStr)));
                
                %===========================================
                % compute Moskowitz
                tic
                result = pbEnv.explSol(tValues,xValues);
                
                obj.N.(linkStr) = result{1};
                activeComp = result{2};
                obj.k.(linkStr) = pbEnv.density(tValues,xValues,activeComp);
                toc
            end
            
        end
        
        
        %===============================================================
        % plot link
        function plotLinks(obj,network,links)
            global K_C_PARA K_M_PARA
            
            if strcmp(links,'all')
                
                links = 1:network.num_link;
                
            end
            
            for i = 1:network.num_link
                linkStr = sprintf('link_%d',links(i));
                
                %===========================================
                %Transformation for better color presentation
                %kc=>0.5km, km=>km
                k_c_tmp = network.network_hwy.(linkStr).paras(K_C_PARA);
                k_m_tmp = network.network_hwy.(linkStr).paras(K_M_PARA);
                k_trans = mapping(obj.k.(linkStr), [0 k_c_tmp; k_c_tmp k_m_tmp],...
                    [0 0.5*k_m_tmp; 0.5*k_m_tmp k_m_tmp]);
                
                %===========================================
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 1 scrsz(3) scrsz(4)]);
                title(sprintf('Link No. %d',i),'fontsize',24);
                
                [~, ~] = LH_plot2D(obj.tScale, obj.postmScale.(linkStr),...
                    obj.N.(linkStr),k_trans, obj.fd.(linkStr));
                set(gca,'fontsize',20)
                xlabel({'time (s)'},'fontsize',24);
                ylabel({'x (m)'},'fontsize',24);
                
            end
               
        end
        
        
        %===============================================================
        % plot junc
        % middle point of the junction
        % middleX = [1700m, 1700m; 1400m 1400m]
        function plotJuncs(obj,network,program,juncs,T,num_steps)
            global K_C_PARA K_M_PARA NUM_DOWN NUM_UP
            
            if strcmp(juncs,'all')
                juncs = 1:network.num_junc;
            end
            
            for i = 1:network.num_junc
                juncStr = sprintf('junc_%d',juncs(i));
                
                if strcmp(network.network_junc.(juncStr).type_junc,'merge')
                    %===========================================
                    % compute the total through flow and penalty for
                    % title
                    inLinks = network.network_junc.(juncStr).inlabel;
                    
                    q_s1 = obj.x(program.Dec(inLinks(1),NUM_DOWN,1):...
                        program.Dec(inLinks(1),NUM_DOWN,2) );
                    q_s2 = obj.x(program.Dec(inLinks(2),NUM_DOWN,1):...
                        program.Dec(inLinks(2),NUM_DOWN,2) );
                    
                    tt_flow = sum((q_s1+q_s2).*T);
                    % L1 norm penalty
                    tt_pen = sum( abs(q_s1 - q_s2*network.network_junc.(juncStr).ratio(1)/...
                        network.network_junc.(juncStr).ratio(2)) );
                    
                    %===========================================
                    %Normalize density
                    %kc=>0.5, km=>1
                    link1 = sprintf('link_%d',network.network_junc.(juncStr).inlabel(1));
                    k_c_tmp = network.network_hwy.(link1).paras(K_C_PARA);
                    k_m_tmp = network.network_hwy.(link1).paras(K_M_PARA);
                    kNorm1 = mapping(obj.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    link2 = sprintf('link_%d',network.network_junc.(juncStr).inlabel(2));
                    k_c_tmp = network.network_hwy.(link2).paras(K_C_PARA);
                    k_m_tmp = network.network_hwy.(link2).paras(K_M_PARA);
                    kNorm2 = mapping(obj.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    link3 = sprintf('link_%d',network.network_junc.(juncStr).outlabel);
                    k_c_tmp = network.network_hwy.(link3).paras(K_C_PARA);
                    k_m_tmp = network.network_hwy.(link3).paras(K_M_PARA);
                    kNorm3 = mapping(obj.k.(link3),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    kLeft = [kNorm1 kNorm3];
                    kRight = [kNorm2 kNorm3];
                    
                    NLeft = [obj.N.(link1) obj.N.(link3)];
                    NRight = [obj.N.(link2) obj.N.(link3)];
                    
                    xScaleLeft = [obj.postmScale.(link1)...
                        obj.postmScale.(link3) + max(obj.postmScale.(link1)) + obj.dx_res/1000 ];
                    xScaleRight = [obj.postmScale.(link2)...
                        obj.postmScale.(link3) + max(obj.postmScale.(link2)) + obj.dx_res/1000];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    [~, ~] = LH_plot2D_junc_dis([network.network_hwy.(link1).profile(1)*1000,...
                        network.network_hwy.(link2).profile(1)*1000],...
                        [0, cumsum(T')], num_steps, obj.tScale, xScaleLeft, xScaleRight,...
                        NLeft, kLeft, NRight, kRight, network.network_junc.(juncStr));
                    h = suptitle(sprintf('Total Flow: %f\n Number of steps: %d',tt_flow, length(T)));
                    set(h,'FontSize',24)
                    
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                    %===========================================
                    % compute the total through flow and penalty for
                    % title
                    outLinks = network.network_junc.(juncStr).outlabel;
                    
                    % extract the inflow of the two outcoming links
                    q_r1 = obj.x(program.Dec(outLinks(1),NUM_UP,1):...
                        program.Dec(outLinks(1),NUM_UP,2) );
                    q_r2 = obj.x(program.Dec(outLinks(2),NUM_UP,1):...
                        program.Dec(outLinks(2),NUM_UP,2) );
                    
                    tt_flow = sum((q_r1+q_r2).*T);
                    % L1 norm penalty
                    tt_pen = sum( abs(q_r1 - q_r2*network.network_junc.(juncStr).ratio(1)/...
                        network.network_junc.(juncStr).ratio(2)) );
                    
                    
                    %===========================================
                    %Normalize density
                    %kc=>0.5, km=>1
                    link1 = sprintf('link_%d',network.network_junc.(juncStr).inlabel);
                    k_c_tmp = network.network_hwy.(link1).paras(K_C_PARA);
                    k_m_tmp = network.network_hwy.(link1).paras(K_M_PARA);
                    kNorm1 = mapping(obj.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    link2 = sprintf('link_%d',network.network_junc.(juncStr).outlabel(1));
                    k_c_tmp = network.network_hwy.(link2).paras(K_C_PARA);
                    k_m_tmp = network.network_hwy.(link2).paras(K_M_PARA);
                    kNorm2 = mapping(obj.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    link3 = sprintf('link_%d',network.network_junc.(juncStr).outlabel(2));
                    k_c_tmp = network.network_hwy.(link3).paras(K_C_PARA);
                    k_m_tmp = network.network_hwy.(link3).paras(K_M_PARA);
                    kNorm3 = mapping(obj.k.(link3),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    kLeft = [kNorm1 kNorm2];
                    kRight = [kNorm1 kNorm3];
                    
                    NLeft = [obj.N.(link1) obj.N.(link2)];
                    NRight = [obj.N.(link1) obj.N.(link3)];
                    
                    xScaleLeft = [obj.postmScale.(link1)...
                        obj.postmScale.(link2) + max(obj.postmScale.(link1))];
                    xScaleRight = [obj.postmScale.(link1)...
                        obj.postmScale.(link3) + max(obj.postmScale.(link1))];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    [~, ~] = LH_plot2D_junc_dis([network.network_hwy.(link1).profile(1)*1000,...
                        network.network_hwy.(link1).profile(1)*1000],...
                        T, num_steps, obj.tScale, xScaleLeft, xScaleRight,...
                        NLeft, kLeft, NRight, kRight,network.network_junc.(juncStr));
                    h = suptitle(sprintf('Total Flow: %f\n Number of steps: %d',tt_flow, length(T)));
                    set(h,'FontSize',24)
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'connection')
                    %===========================================
                    outLink = network.network_junc.(juncStr).outlabel;
                    q_thru = obj.x(program.Dec(outLink,NUM_UP,1):...
                        program.Dec(outLink,NUM_UP,2) );
                    tt_flow = sum((q_thru).*T);
                    
                    %===========================================
                    %Normalize density
                    %kc=>0.5, km=>1
                    link1 = sprintf('link_%d',network.network_junc.(juncStr).inlabel);
                    k_c_tmp = network.network_hwy.(link1).paras(K_C_PARA);
                    k_m_tmp = network.network_hwy.(link1).paras(K_M_PARA);
                    kNorm1 = mapping(obj.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    link2 = sprintf('link_%d',network.network_junc.(juncStr).outlabel);
                    k_c_tmp = network.network_hwy.(link2).paras(K_C_PARA);
                    k_m_tmp = network.network_hwy.(link2).paras(K_M_PARA);
                    kNorm2 = mapping(obj.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    kComb = [kNorm1 kNorm2];
                    
                    NComb = [obj.N.(link1) obj.N.(link2)];
                    
                    xScaleComb = [obj.postmScale.(link1)...
                        obj.postmScale.(link2) + max(obj.postmScale.(link1))];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    [~, ~] = LH_plot2D_junc_dis(network.network_hwy.(link1).profile(1)*1000,...
                        T, num_steps, obj.tScale, xScaleComb,...
                        NComb, kComb,network.network_junc.(juncStr));
                    h = suptitle(sprintf('Total Flow: %f \n Number of steps: %d',tt_flow, length(T)));
                    set(h,'FontSize',24)
                end
                
            end
            
        end     %end function
        
        
        %===============================================================
        % analyze junction, stitch two incoming (or out going) links
        % together to check the distribution rule of priority rule
        function checkJunc(obj,x,network,program,juncs)
            global K_C_PARA K_M_PARA NUM_DOWN NUM_UP
            if strcmp(juncs,'all')
                
                for i = 1:network.num_junc
                    juncStr = sprintf('junc_%d',i);
                    
                    if strcmp(network.network_junc.(juncStr).type_junc,'merge')
                        
                        %===========================================
                        %Normalize density
                        %kc=>0.5, km=>1
                        link1 = sprintf('link_%d',network.network_junc.(juncStr).inlabel(1));
                        k_c_tmp = network.network_hwy.(link1).paras(K_C_PARA);
                        k_m_tmp = network.network_hwy.(link1).paras(K_M_PARA);
                        kNorm1 = mapping(obj.k.(link1),...
                            [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                            [0 0.5; 0.6 1]);
                        
                        link2 = sprintf('link_%d',network.network_junc.(juncStr).inlabel(2));
                        k_c_tmp = network.network_hwy.(link2).paras(K_C_PARA);
                        k_m_tmp = network.network_hwy.(link2).paras(K_M_PARA);
                        kNorm2 = mapping(obj.k.(link2),...
                            [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                            [0 0.5; 0.6 1]);
                        
                        % Flip the density and trace for link 2
                        kNorm2(:,1:size(kNorm2,2)) = kNorm2(:,size(kNorm2,2):-1:1);
                        Ntmp = obj.N.(link2)(:, size(obj.N.(link2),2):-1:1);
                        
                        kComb = [kNorm1 kNorm2];
                        
                        NComb = [obj.N.(link1) Ntmp];
                        
                        xScaleComb = [obj.postmScale.(link1)...
                            obj.postmScale.(link2) + max(obj.postmScale.(link1))];
                        
                        %===========================================
                        scrsz = get(0,'ScreenSize');
                        figure('Position',[1 1 scrsz(3) scrsz(4)]);
                        %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                        
                        [C, h] = LH_plot2D_junc(obj.tScale, xScaleComb,...
                            NComb, kComb, network.network_junc.(juncStr));
                        %                        set(gca,'fontsize',20)
                        %                        xlabel({'time (s)'},'fontsize',24);
                        %                        ylabel({'x (m)'},'fontsize',24);
                        
                        disp('      q1      q2      q3      q1/q3       q2/q3');
                        disp([x(program.Dec(1,NUM_DOWN,1):program.Dec(1,NUM_DOWN,2)) ...
                            x(program.Dec(2,NUM_DOWN,1):program.Dec(2,NUM_DOWN,2)) ...
                            x(program.Dec(3,NUM_UP,1):program.Dec(3,NUM_UP,2)) ...
                            x(program.Dec(1,NUM_DOWN,1):program.Dec(1,NUM_DOWN,2))./...
                            x(program.Dec(3,NUM_UP,1):program.Dec(3,NUM_UP,2)) ...
                            x(program.Dec(2,NUM_DOWN,1):program.Dec(2,NUM_DOWN,2))./...
                            x(program.Dec(3,NUM_UP,1):program.Dec(3,NUM_UP,2)) ]);
                        
                        
                    elseif strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                        
                        %===========================================
                        %Normalize density
                        %kc=>0.5, km=>1
                        link1 = sprintf('link_%d',network.network_junc.(juncStr).outlabel(1));
                        k_c_tmp = network.network_hwy.(link1).paras(K_C_PARA);
                        k_m_tmp = network.network_hwy.(link1).paras(K_M_PARA);
                        kNorm1 = mapping(obj.k.(link1),...
                            [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                            [0 0.5; 0.6 1]);
                        
                        link2 = sprintf('link_%d',network.network_junc.(juncStr).outlabel(2));
                        k_c_tmp = network.network_hwy.(link2).paras(K_C_PARA);
                        k_m_tmp = network.network_hwy.(link2).paras(K_M_PARA);
                        kNorm2 = mapping(obj.k.(link2),...
                            [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                            [0 0.5; 0.6 1]);
                        
                        % Flip the density and trace for link 1
                        kNorm1(:,1:size(kNorm1,2)) = kNorm1(:,size(kNorm1,2):-1:1);
                        Ntmp = obj.N.(link1)(:, size(obj.N.(link1),2):-1:1);
                        
                        kComb = [kNorm1 kNorm2];
                        
                        NComb = [Ntmp obj.N.(link2)];
                        
                        xScaleComb = [obj.postmScale.(link1)...
                            obj.postmScale.(link2) + max(obj.postmScale.(link1))];
                        
                        %===========================================
                        scrsz = get(0,'ScreenSize');
                        figure('Position',[1 1 scrsz(3) scrsz(4)]);
                        %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                        
                        [C, h] = LH_plot2D_junc(obj.tScale, xScaleComb,...
                            NComb, kComb, network.network_junc.(juncStr));
                        %                        set(gca,'fontsize',20)
                        %                        xlabel({'time (s)'},'fontsize',24);
                        %                        ylabel({'x (m)'},'fontsize',24);
                        
                    end
                    
                end     %end for
                
            else
                
                for i = 1:length(juncs)
                    juncStr = sprintf('junc_%d',juncs(i));
                    
                    if strcmp(network.network_junc.(juncStr).type_junc,'merge')
                        
                        %===========================================
                        %Normalize density
                        %kc=>0.5, km=>1
                        link1 = sprintf('link_%d',network.network_junc.(juncStr).inlabel(1));
                        k_c_tmp = network.network_hwy.(link1).paras(K_C_PARA);
                        k_m_tmp = network.network_hwy.(link1).paras(K_M_PARA);
                        kNorm1 = mapping(obj.k.(link1),...
                            [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                            [0 0.5; 0.6 1]);
                        
                        link2 = sprintf('link_%d',network.network_junc.(juncStr).inlabel(2));
                        k_c_tmp = network.network_hwy.(link2).paras(K_C_PARA);
                        k_m_tmp = network.network_hwy.(link2).paras(K_M_PARA);
                        kNorm2 = mapping(obj.k.(link2),...
                            [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                            [0 0.5; 0.6 1]);
                        
                        % Flip the density and trace for link 2
                        kNorm2(:,1:size(kNorm2,2)) = kNorm2(:,size(kNorm2,2):-1:1);
                        Ntmp = obj.N.(link2)(:, size(obj.N.(link2),2):-1:1);
                        
                        kComb = [kNorm1 kNorm2];
                        
                        NComb = [obj.N.(link1) Ntmp];
                        
                        xScaleComb = [obj.postmScale.(link1)...
                            obj.postmScale.(link2) + max(obj.postmScale.(link1))];
                        
                        %===========================================
                        scrsz = get(0,'ScreenSize');
                        figure('Position',[1 1 scrsz(3) scrsz(4)]);
                        %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                        
                        [C, h] = LH_plot2D_junc(obj.tScale, xScaleComb,...
                            NComb, kComb, network.network_junc.(juncStr));
                        %                        set(gca,'fontsize',20)
                        %                        xlabel({'time (s)'},'fontsize',24);
                        %                        ylabel({'x (m)'},'fontsize',24);
                        
                        disp('      q1      q2      q3      q1/q3       q2/q3');
                        disp([x(program.Dec(1,NUM_DOWN,1):program.Dec(1,NUM_DOWN,2)) ...
                            x(program.Dec(2,NUM_DOWN,1):program.Dec(2,NUM_DOWN,2)) ...
                            x(program.Dec(3,NUM_UP,1):program.Dec(3,NUM_UP,2)) ...
                            x(program.Dec(1,NUM_DOWN,1):program.Dec(1,NUM_DOWN,2))./...
                            x(program.Dec(3,NUM_UP,1):program.Dec(3,NUM_UP,2)) ...
                            x(program.Dec(2,NUM_DOWN,1):program.Dec(2,NUM_DOWN,2))./...
                            x(program.Dec(3,NUM_UP,1):program.Dec(3,NUM_UP,2)) ]);
                        
                    elseif strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                        
                        %===========================================
                        %Normalize density
                        %kc=>0.5, km=>1
                        link1 = sprintf('link_%d',network.network_junc.(juncStr).outlabel(1));
                        k_c_tmp = network.network_hwy.(link1).paras(K_C_PARA);
                        k_m_tmp = network.network_hwy.(link1).paras(K_M_PARA);
                        kNorm1 = mapping(obj.k.(link1),...
                            [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                            [0 0.5; 0.6 1]);
                        
                        link2 = sprintf('link_%d',network.network_junc.(juncStr).outlabel(2));
                        k_c_tmp = network.network_hwy.(link2).paras(K_C_PARA);
                        k_m_tmp = network.network_hwy.(link2).paras(K_M_PARA);
                        kNorm2 = mapping(obj.k.(link2),...
                            [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                            [0 0.5; 0.6 1]);
                        
                        % Flip the density and trace for link 1
                        kNorm1(:,1:size(kNorm1,2)) = kNorm1(:,size(kNorm1,2):-1:1);
                        Ntmp = obj.N.(link1)(:, size(obj.N.(link1),2):-1:1);
                        
                        kComb = [kNorm1 kNorm2];
                        
                        NComb = [Ntmp obj.N.(link2)];
                        
                        xScaleComb = [obj.postmScale.(link1)...
                            obj.postmScale.(link2) + max(obj.postmScale.(link1))];
                        
                        %===========================================
                        scrsz = get(0,'ScreenSize');
                        figure('Position',[1 1 scrsz(3) scrsz(4)]);
                        %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                        
                        [C, h] = LH_plot2D_junc(obj.tScale, xScaleComb,...
                            NComb, kComb, network.network_junc.(juncStr));
                        %                        set(gca,'fontsize',20)
                        %                        xlabel({'time (s)'},'fontsize',24);
                        %                        ylabel({'x (m)'},'fontsize',24);
                        
                    end
                    
                end     %end for
                
            end
            
        end     %end function
        
        
        %===============================================================
        % analytically check if the solution is entropy sample those points
        % and compare with LP solutions
        % input: x, (solution)
        %        network,
        %        program,
        %        tolerance for entrpy (e.g 10 cars)
        % output: [TF, steps], (true, false), (steps that are not entropy)
        function [TF, steps] = checkEntropy(obj, x, network, program, entropyTol)
            
            global NUM_UP NUM_DOWN
            
            TF = true;
            steps = [];
            
            % check each junction
            for junc = 1:network.num_junc
                juncStr = sprintf('junc_%d',junc);
                
                % Here consider connection case
                if strcmp(network.network_junc.(juncStr).type_junc,'connection')
                    
                    outLink = network.network_junc.(juncStr).outlabel;
                    
                    % extract the thru flow value which is downstream
                    % inflow in the connection case
                    q_thru = x(program.Dec(outLink,NUM_UP,1):...
                        program.Dec(outLink,NUM_UP,2) );
                    
                    % check each step
                    for i = 1: program.num_steps
                        
                        t_ref = program.T_cum(i);   % start time of this point
                        t_end = program.T_cum(i+1); % end time of this point
                        
                        % Here the sampling identifier is 0, meaning return
                        % the true entropy solution
                        d_M = obj.samplePointsJunc(t_ref, t_end, junc, 0, x, network, program.Dec);
                        
                        if abs(q_thru(i)*program.T(i)-d_M) > entropyTol
                            % difference greater than the tolerance
                            TF = false;
                            steps = [ steps; i];
                        else
                            % double check the middle point
                            d_M_middle = obj.samplePointsJunc(t_ref,...
                                t_ref+(t_end-t_ref)/2, junc, 0, x, network, program.Dec);
                            
                            if abs(q_thru(i)*program.T(i)/2 - d_M_middle) > entropyTol/2 ||...
                                    abs(q_thru(i)*program.T(i)/2 - (d_M - d_M_middle) ) > entropyTol/2
                                % not entropy
                                TF = false;
                                steps = [steps; i];
                            end
                        end
                    end     %end each step
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'merge')
                    
                    inLinks = network.network_junc.(juncStr).inlabel;
                    
                    % extract the outflow of the two incoming links
                    q_s1 = x(program.Dec(inLinks(1),NUM_DOWN,1):...
                        program.Dec(inLinks(1),NUM_DOWN,2) );
                    q_s2 = x(program.Dec(inLinks(2),NUM_DOWN,1):...
                        program.Dec(inLinks(2),NUM_DOWN,2) );
                    
                    % check each step
                    for i = 1: program.num_steps
                        
                        t_ref = program.T_cum(i);   % start time of this point
                        t_end = program.T_cum(i+1); % end time of this point
                        
                        % a 1 x 2 vector true unique solution from the two incoming links
                        d_M = obj.samplePointsJunc(t_ref, t_end, junc, 0, x, network, program.Dec);
                        
                        if abs( q_s1(i)*program.T(i)-d_M(1) ) > entropyTol ||...
                                abs( q_s2(i)*program.T(i)-d_M(2) ) > entropyTol
                            % does not match the unique solution from the
                            % sampling approach
                            TF = false;
                            steps = [ steps; i];
                            
                        else
                            % doubel check the middle points
                            d_M_middle = obj.samplePointsJunc(t_ref,...
                                t_ref+(t_end-t_ref)/2, junc, 0, x, network, program.Dec);
                            
                            if abs( q_s1(i)*program.T(i)/2 - d_M_middle(1) ) > entropyTol/2 ||...
                                abs( q_s2(i)*program.T(i)/2 - d_M_middle(2) ) > entropyTol/2
                            TF = false;
                                steps = [steps; i];
                            end
                            
                        end
                    end     % end each step
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                    
                    outLinks = network.network_junc.(juncStr).outlabel;
                    
                    % extract the inflow of the two outcoming links
                    q_r1 = x(program.Dec(outLinks(1),NUM_UP,1):...
                        program.Dec(outLinks(1),NUM_UP,2) );
                    q_r2 = x(program.Dec(outLinks(2),NUM_UP,1):...
                        program.Dec(outLinks(2),NUM_UP,2) );
                    
                    % check each step
                    for i = 1: program.num_steps
                        
                        t_ref = program.T_cum(i);   % start time of this point
                        t_end = program.T_cum(i+1); % end time of this point
                        
                        % a 1 x 2 vector true unique solution from the two incoming links
                        d_M = obj.samplePointsJunc(t_ref, t_end, junc, 0, x, network, program.Dec);
                        
                        if abs( q_r1(i)*program.T(i)-d_M(1) ) > entropyTol ||...
                                abs( q_r2(i)*program.T(i)-d_M(2) ) > entropyTol
                            % does not match the unique solution from the
                            % sampling approach
                            TF = false;
                            steps = [ steps; i];
                            
                        else
                            % doubel check the middle points
                            d_M_middle = obj.samplePointsJunc(t_ref,...
                                t_ref+(t_end-t_ref)/2, junc, 0, x, network, program.Dec);
                            
                            if abs( q_r1(i)*program.T(i)/2 - d_M_middle(1) ) > entropyTol/2 ||...
                                abs( q_r2(i)*program.T(i)/2 - d_M_middle(2) ) > entropyTol/2
                            TF = false;
                                steps = [steps; i];
                            end
                            
                        end
                    end     % end each step
                    
                end     % end diverge
                
                
            end % end each junction
            
            steps = unique(steps);
            
        end % end function
        
        
        %===============================================================
        % analytically update the discretization by identify flow corners.
        % return T_cum,  with T(1) = 0 and last(T) = length Simulation Time.
        function T_cum = updateTimeDiscretization(obj, x, steps, network, program)
            
            searchDepth = 0;
            d_t = 1.0e-3;    % 0.001 second
            T_new = program.T_cum;
            
            
            % check each junction
            for junc = 1:network.num_junc
                juncStr = sprintf('junc_%d',junc);
                
                if strcmp( network.network_junc.(juncStr).type_junc, 'connection')
                    % find the interval from the begining of this nonentropy step
                    % to the end of the following step. (Since the nonentropy could
                    % be caused by the intersection from following step)
                    t_left = program.T_cum(steps(1));
                    if steps(1) < program.num_steps
                        t_right = program.T_cum(steps(1)+2);
                    else
                        t_right = program.T_cum(steps(1)+1);
                    end
                    
                    % query the values of those two points with reference to t_left
                    % sending and receiving
                    M_left = 0;
                    M_right = obj.samplePointsJunc(t_left, t_right, junc, 4, x, network, program.Dec);
                    
                    t_found_s = searchIntersection(obj, [t_left, t_right], [M_left, M_right(1)], ...
                        searchDepth, -1, -1, d_t, t_left, junc, 1, x, network, program.Dec);
                    t_found_r = searchIntersection(obj, [t_left, t_right], [M_left, M_right(2)], ...
                        searchDepth, -1, -1, d_t, t_left, junc, 2, x, network, program.Dec);
                    
                    t_found = [t_found_s; t_found_r];
                    
                    for t=1:length(t_found)
                        
                        % add to T_cum
                        if t_found(t) < program.end_time
                            T_new( length(T_new)+1 ) = t_found(t);
                        end
                        
                    end
                    
                    % return sorted and unique discretization accumulative points
                    T_cum = unique_tol(T_new, 1.0e-8);
                    
                elseif strcmp( network.network_junc.(juncStr).type_junc, 'merge')
                    
                    % The merge is the same as the connection case
                    t_left = program.T_cum(steps(1));
                    if steps(1) < program.num_steps
                        t_right = program.T_cum(steps(1)+2);
                    else
                        t_right = program.T_cum(steps(1)+1);
                    end
                    
                    % query the values of those two points with reference to t_left
                    M_left = 0;
                    
                    % Here M_right is the NOT the entropy solution; but the
                    % sending M of incoming links
                    M_right = obj.samplePointsJunc(t_left, t_right, junc, 4, x, network, program.Dec);
                    
                    % corner points found from sending link 1, 2 and
                    % receiving link 3
                    t_found_s1 = searchIntersection(obj, [t_left, t_right], [M_left, M_right(1)], ...
                        searchDepth, -1, -1, d_t, t_left, junc, 1, x, network, program.Dec);
                    t_found_s2 = searchIntersection(obj, [t_left, t_right], [M_left, M_right(2)], ...
                        searchDepth, -1, -1, d_t, t_left, junc, 2, x, network, program.Dec);
                    t_found_r = searchIntersection(obj, [t_left, t_right], [M_left, M_right(3)], ...
                        searchDepth, -1, -1, d_t, t_left, junc, 3, x, network, program.Dec);
                   
                    
                    t_found = [t_found_s1; t_found_s2; t_found_r];
                    for t=1:length(t_found)
                        
                        % add to T_cum
                        if t_found(t) < program.end_time
                            T_new( length(T_new)+1 ) = t_found(t);
                        end
                        
                    end
                    
                    % return sorted and unique discretization accumulative points
                    T_cum = unique_tol(T_new, 1.0e-8);
                    
                elseif strcmp( network.network_junc.(juncStr).type_junc, 'diverge')
                     % Find the potential range for the intersection
                    t_left = program.T_cum(steps(1));
                    if steps(1) < program.num_steps
                        t_right = program.T_cum(steps(1)+2);
                    else
                        t_right = program.T_cum(steps(1)+1);
                    end
                    
                    % query the values of those two points with reference to t_left
                    M_left = 0;
                    
                    % Here M_right is the NOT the entropy solution; but the
                    % sending and receiving M of incoming and outgoing
                    % links
                    M_right = obj.samplePointsJunc(t_left, t_right, junc, 4, x, network, program.Dec);
                    
                    % corner points found from sending link 1 and
                    % receiving link 2, 3
                    t_found_s = searchIntersection(obj, [t_left, t_right], [M_left, M_right(1)], ...
                        searchDepth, -1, -1, d_t, t_left, junc, 3, x, network, program.Dec);
                    t_found_r1 = searchIntersection(obj, [t_left, t_right], [M_left, M_right(2)], ...
                        searchDepth, -1, -1, d_t, t_left, junc, 1, x, network, program.Dec);
                    t_found_r2 = searchIntersection(obj, [t_left, t_right], [M_left, M_right(3)], ...
                        searchDepth, -1, -1, d_t, t_left, junc, 2, x, network, program.Dec);
                   
                    
                    t_found = [t_found_s; t_found_r1; t_found_r2];
                    for t=1:length(t_found)
                        
                        % add to T_cum
                        if t_found(t) < program.end_time
                            T_new( length(T_new)+1 ) = t_found(t);
                        end
                        
                    end
                    
                    % return sorted and unique discretization accumulative points
                    T_cum = unique_tol(T_new, 1.0e-8);
                end
                
                
            end
            
        end
        
        
        
        %===============================================================
        % This is a recursive function aiming at finding the corner point
        % of a piecewise linear nondecreasing function in an interval
        % input: t_interval, (time interval)
        %        M_interval, (function value of start and end point)
        %        searchDepth, (maximal recursion depth
        %        s_left/ s_right, (left/right derivative of the function)
        %        d_t, (maximal time resolution we need)
        %        (t_ref, junc, x, network, Dec) are used to sample points
        %        linkIndex tells which value to returnwhen sampling
        % output: found time points 1 x n column
        function t_found = searchIntersection(obj, t_interval, M_interval, ...
                searchDepth,s_left, s_right, d_t, t_ref, junc, linkIndex, x, network, Dec)
            
            % if "recursion depth >=0" and "is one straight line"
            % return [] immediatly if it seems like to be a straight line
            if searchDepth >= 0 &&...
                    max(s_left,s_right) >= 0 &&...
                    abs(max(s_left,s_right) -...
                    (M_interval(2)-M_interval(1))/(t_interval(2)-t_interval(1)))...
                    <=1.0e-3
                % finished checking if it is just a straight line and checked 2 depths
                t_found = [];
                return
            end
            
            % if "not a straight line" and ("searchDepth >=4" or "len_interval<=d_t")
            % return middle point as an approximation
            if ( max(s_left,s_right) >= 0 &&...
                    abs(max(s_left,s_right) -...
                    (M_interval(2)-M_interval(1))/(t_interval(2)-t_interval(1)))...
                    >=1.0e-3) &&...
                    ( searchDepth >= 8   ||...
                    (t_interval(2)-t_interval(1))<=d_t )
                
                t_found = mean(t_interval);
                return
            end
            
            % left and right points
            t_A = t_interval(1);
            t_B = t_interval(2);
            M_A = M_interval(1);
            M_B = M_interval(2);
            
            % algorithm starts by pick 3 points
            t_L1 = t_A + (t_B-t_A)/4;
            t_L2 = t_A + (t_B-t_A)/2;
            t_R1 = t_A + 3*(t_B-t_A)/4;
            t_R2 = t_A + (t_B-t_A)/2;
            
            %figure(2)
            %scatter(t_L2, obj.samplePointsJunc(t_ref,t_L2,junc,x,network,Dec), 'b*');
            leftFound = false;
            rightFound = false;
            
            t_found = [];
            
            % Use binary search to find the left slope
            if s_left ~= -1
                sL = s_left;
            else
                % get leftside slope
                while leftFound == false
                    %scatter(t_L1, obj.samplePointsJunc(t_ref, t_L1, junc, x, network, Dec), 'b*');
                    
                    if abs((obj.samplePointsJunc(t_ref, t_L1, junc, linkIndex, x, network, Dec)-M_A)/(t_L1-t_A) -...
                            (obj.samplePointsJunc(t_ref, t_L2, junc, linkIndex, x, network, Dec)-M_A)/(t_L2-t_A)) <= 1.0e-8
                        % if find the exact slope
                        leftFound = true;
                        sL = (obj.samplePointsJunc(t_ref, t_L2, junc, linkIndex, x, network, Dec)-M_A)/(t_L2-t_A);
                        
                    elseif (t_L1-t_A <= d_t)
                        % one more trial
                        t_L0 = (t_A+t_L1)/2;
                        if abs( (obj.samplePointsJunc(t_ref, t_L0, junc, linkIndex, x, network, Dec)...
                                -obj.samplePointsJunc(t_ref, t_A, junc, linkIndex, x, network, Dec))/(t_L0-t_A) -...
                                (obj.samplePointsJunc(t_ref, t_L1, junc, linkIndex, x, network, Dec)...
                                -obj.samplePointsJunc(t_ref, t_A, junc, linkIndex, x, network, Dec))/(t_L1-t_A) ) < 1.0e-8
                            
                            % found the left slope
                            leftFound = true;
                            sL = (obj.samplePointsJunc(t_ref, t_L1, junc, linkIndex, x, network, Dec)...
                                -obj.samplePointsJunc(t_ref, t_A, junc, linkIndex, x, network, Dec))/(t_L1-t_A);
                        else
                            % the coner point is between t_A and t_L1
                            % approximate the coner point
                            leftFound = true;
                            % push approximated point to array, use sign to mark this is an
                            % approximated point
                            t_found = [t_found; (t_A+t_L1)/2];
                            sL = (obj.samplePointsJunc(t_ref, t_L2, junc, linkIndex, x, network, Dec)...
                                -obj.samplePointsJunc(t_ref, t_L1, junc, linkIndex, x, network, Dec))/(t_L2-t_L1);
                            t_A = t_L1; % adjust t_A
                        end
                    else
                        t_L2 = t_L1;
                        t_L1 = t_A+(t_L1-t_A)/2;
                    end
                end
            end
            
            % Use binary search to find the right slope
            if s_right ~= -1
                sR = s_right;
            else
                while rightFound == false
                    %scatter(t_R1, obj.samplePointsJunc(t_ref, t_R1, junc, x, network, Dec), 'b*');
                    
                    if abs((obj.samplePointsJunc(t_ref, t_R1, junc, linkIndex, x, network, Dec)-M_B)/(t_R1-t_B) -...
                            (obj.samplePointsJunc(t_ref, t_R2, junc, linkIndex, x, network, Dec)-M_B)/(t_R2-t_B)) <= 1.0e-8
                        
                        rightFound = true;
                        sR = (obj.samplePointsJunc(t_ref, t_R2, junc, linkIndex, x, network, Dec)-M_B)/(t_R2-t_B);
                    elseif t_B-t_R1 <= d_t
                        
                        % one more trial
                        t_R0 = (t_B+t_R1)/2;
                        if abs( (obj.samplePointsJunc(t_ref, t_R0, junc, linkIndex, x, network, Dec)...
                                -obj.samplePointsJunc(t_ref, t_B, junc, linkIndex, x, network, Dec))/(t_R0-t_B) -...
                                (obj.samplePointsJunc(t_ref, t_R1, junc, linkIndex, x, network, Dec)...
                                -obj.samplePointsJunc(t_ref, t_B, junc, linkIndex, x, network, Dec))/(t_R1-t_B) ) < 1.0e-8
                            
                            % found the left slope
                            rightFound = true;
                            sR = (obj.samplePointsJunc(t_ref, t_R1, junc, linkIndex, x, network, Dec)...
                                -obj.samplePointsJunc(t_ref, t_B, junc, linkIndex, x, network, Dec))/(t_R1-t_B);
                        else
                            % if the point is too close to the bound point
                            % approximate
                            rightFound = true;
                            % use sign to mark this is an approximated corner point
                            t_found = [t_found; (t_B+t_R1)/2];
                            sR = (obj.samplePointsJunc(t_ref, t_R2, junc, linkIndex, x, network, Dec)...
                                -obj.samplePointsJunc(t_ref, t_R1, junc, linkIndex, x, network, Dec))/(t_R2-t_R1);
                            t_B = t_R1; % update right boundary.
                        end
                    else
                        t_R2 = t_R1;
                        t_R1 = t_R1+(t_B-t_R1)/2;
                    end
                end
            end
            
            % define the new left and right bounds for potential split
            if s_left ~= -1
                t_newLeft = t_A;
                t_splitLeft = t_A;
            else
                t_newLeft = t_L1;
                t_splitLeft = t_L2;
            end
            
            if s_right ~= -1
                t_newRight = t_B;
                t_splitRight = t_B;
            else
                t_newRight = t_R1;
                t_splitRight = t_R2;
            end
            
            % compute the intersection point
            % make sure not devided by 0
            if (abs(sL - sR) >= 1.0e-3)
                t_intersection = ((M_B-sR*t_B)-(M_A-sL*t_A))/(sL-sR);
            end
            
            % if "lines are not identical/parallel" and
            % "t_intersection is in the interval" and
            % "t_intersection-d_t and t_intersection has slope sL" and
            % "t_intersection+d_t and t_intersection has slope sR" and
            if ( (abs(sL - sR) >= 1.0e-3)&&...
                    ( ( t_intersection >= t_newLeft ) && ( t_intersection <= t_newRight) ) &&...
                    (abs(obj.samplePointsJunc(t_ref, t_intersection, junc, linkIndex, x, network, Dec) - ...
                    (sL*t_intersection+M_A-sL*t_A)  ) <=1.0e-8 ) &&...
                    (abs( ( obj.samplePointsJunc(t_ref, t_intersection, junc, linkIndex, x, network, Dec) -...
                    obj.samplePointsJunc(t_ref, t_intersection-d_t, junc, linkIndex, x, network, Dec) )/d_t - sL) <= 1.0e-8 ) &&...
                    (abs( ( obj.samplePointsJunc(t_ref, t_intersection+d_t, junc, linkIndex, x, network, Dec) -...
                    obj.samplePointsJunc(t_ref, t_intersection, junc, linkIndex, x, network, Dec) )/d_t - sR) <= 1.0e-8) )
                
                t_found = [t_found; ((M_B-sR*t_B)-(M_A-sL*t_A))/(sL-sR)];
                
                
                % there are two cases: it is one straight line, there are
                % corner points in the interval; Maybe we can use this
                % knowledge and choose the split point more wisely
                %             elseif (abs(sL - sR) <= 1.0e-3) &&...
                %                     ( abs( (M_B-M_A)/(t_B-t_A) -sL) <=1.0e-3 )
            else
                
                t_split = (t_splitLeft+t_splitRight)/2;
                
                % split the interval
                %scatter( t_newLeft, obj.samplePointsJunc(t_ref, t_newLeft, junc, x, network, Dec), 'rs');
                %scatter( t_split, obj.samplePointsJunc(t_ref, t_split, junc, x, network, Dec), 'rs');
                %scatter( t_newRight, obj.samplePointsJunc(t_ref, t_newRight, junc, x, network, Dec), 'rs');
                
                t_found = [t_found;...
                    obj.searchIntersection([t_newLeft, t_split],...
                    [obj.samplePointsJunc(t_ref, t_newLeft, junc, linkIndex, x, network, Dec),...
                    obj.samplePointsJunc(t_ref, t_split, junc, linkIndex, x, network, Dec)],...
                    searchDepth+1, sL, -1, d_t,...
                    t_ref, junc, linkIndex, x, network, Dec)  ];
                
                t_found = [t_found;...
                    obj.searchIntersection([t_split, t_newRight],...
                    [obj.samplePointsJunc(t_ref, t_split, junc, linkIndex, x, network, Dec),...
                    obj.samplePointsJunc(t_ref, t_newRight, junc, linkIndex, x, network, Dec)],...
                    searchDepth+1, -1, sR, d_t,...
                    t_ref, junc, linkIndex, x, network, Dec)  ];
                
            end
            
            
        end
        
        
        
        %===============================================================
        % sample points at boundaries at a junction
        % input: t,
        %        t_ref, (the starting time of the nonentropy step, all vehicle
        %        labels are computed with reference to this point)
        %        junc, x, network ( for merges and diverges, maybe need
        %                              priority and distribution factor)
        %        rIndex: 4- return 1x3 vector for merge and diverge
        %                1- return 1st sending or receiving for
        %                merge/diverge
        %                2- return 2nd sending or receiving for
        %                merge/diverge
        %                3- return 3rd link receiving/sending
        %                0- return the entropy solution
        %        
        % output: d_M, ( number of vehicles that can travel from
        %         the begaining of the nonentropy step to time t,
        %         Could be a 1 x 2 vector for the merge and diverge )
        function d_M = samplePointsJunc(obj, t_ref, t, junc, rIndex, x, network, Dec)
            
            juncStr = sprintf('junc_%d',junc);
            
            if strcmp(network.network_junc.(juncStr).type_junc,'connection')
                
                inLink = network.network_junc.(juncStr).inlabel;
                outLink = network.network_junc.(juncStr).outlabel;
                
                tmp_d_M = [obj.samplePointsLink(t_ref, t, inLink,'dn',x,network, Dec) ,...
                    obj.samplePointsLink(t_ref, t, outLink,'up',x,network, Dec)]
                
                % min of supply-demand
                if rIndex == 0
                    d_M = min( tmp_d_M );
                elseif rIndex == 1
                    d_M = tmp_d_M(1);
                elseif rIndex == 2
                    d_M = tmp_d_M(2);
                elseif rIndex == 4
                    d_M = tmp_d_M;
                end
                
                
            elseif strcmp(network.network_junc.(juncStr).type_junc,'merge')
                
                inLinks = network.network_junc.(juncStr).inlabel;
                outLink = network.network_junc.(juncStr).outlabel;
                
                % slope of the priority with y: inLinks(2); and x:
                % inlinks(1)
                sPriority = network.network_junc.(juncStr).ratio(2)/...
                    network.network_junc.(juncStr).ratio(1);
                
                % here we extract the sending and receiving function
                M_R = obj.samplePointsLink(t_ref, t, outLink, 'up', x, network, Dec);
                M_S1 = obj.samplePointsLink(t_ref, t, inLinks(1), 'dn', x, network, Dec);
                M_S2 = obj.samplePointsLink(t_ref, t, inLinks(2), 'dn', x, network, Dec);
                
                % case 1: M_S1 + M_S2 <= M_R, then through flow is
                % M_s1+M_s2
                if M_S1 + M_S2 <= M_R
                    d_M = [M_S1, M_S2];
                    % case 2: intersection point in side (0-M_S1, 0-M_S2) box
                elseif ( M_R/(1+sPriority) <= M_S1) &&...
                        (sPriority*M_R/(1+sPriority) <= M_S2)
                    d_M = [M_R/(1+sPriority), sPriority*M_R/(1+sPriority)];                    
                    % case 3: constrained by M_S1 flow
                elseif M_R/(1+sPriority) > M_S1
                    d_M = [M_S1, M_R-M_S1];
                    % case 4: constrained by M_S2 flow
                elseif sPriority*M_R/(1+sPriority) > M_S2
                    d_M = [M_R-M_S2 , M_S2];
                end
                
                % This is the sending function sampled from
                % samplePointsLink function.
                d_M_sendRev = [M_S1, M_S2, M_R];
                
                % when rIndex == 0, called to check entropy
                % otherwise, just a call to samplePointsLink
                if rIndex == 1
                    d_M = d_M_sendRev(1);
                elseif rIndex == 2
                    d_M = d_M_sendRev(2);
                elseif rIndex == 3
                    d_M = d_M_sendRev(3);
                elseif rIndex == 4
                    d_M = d_M_sendRev;
                end
                
            elseif strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                % TODO, to fix the return index
                
                inLink = network.network_junc.(juncStr).inlabel;
                outLinks = network.network_junc.(juncStr).outlabel;
                
                % slope of the distribution with y: outLinks(2); and x:
                % outlinks(1)
                sDistribution = network.network_junc.(juncStr).ratio(2)/...
                    network.network_junc.(juncStr).ratio(1);
                
                % here we extract the sending and receiving function
                M_R1 = obj.samplePointsLink(t_ref, t, outLinks(1), 'up', x, network, Dec);
                M_R2 = obj.samplePointsLink(t_ref, t, outLinks(2), 'up', x, network, Dec);
                M_S = obj.samplePointsLink(t_ref, t, inLink, 'dn', x, network, Dec);
                
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
                
                d_M_sendRev = [M_R1, M_R2, M_S];
                
                % check which value to return
                if rIndex == 1
                    d_M = d_M_sendRev(1);
                elseif rIndex == 2
                    d_M = d_M_sendRev(2);
                elseif rIndex == 3
                    d_M = d_M_sendRev(3);
                elseif rIndex == 4
                    d_M = d_M_sendRev;
                end
                
            end
            
        end
        
        
        %===============================================================
        % sample points at a link boundary
        % input: t,
        %        t_ref, (the starting time of the nonentropy step, all vehicle
        %        labels are computed with reference to this point)( a vector to be sampled)
        %        link,  ( on which link to be sampled)
        %        'up'/'dn' ( upstream or downstream bound
        %        x, (solution from program)
        %        network, ( which contains information of the link)
        % output: M, ( a vector of vehicle IDs ).
        function d_M = samplePointsLink(obj, t_ref, t, link, bound, x, network, Dec)
            global NUM_UP NUM_DOWN
            global V_PARA W_PARA K_M_PARA K_C_PARA POSTKM
            global LINKLABEL P_MIN DENS
            
            % save the absolute vehicle id on this link in M = [ 0 0 ]
            M = ones(1,length(t)+1)*NaN;
            
            % compute the M of [t_ref t] on this link
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = network.network_hwy.(linkStr).paras(V_PARA);
            w = network.network_hwy.(linkStr).paras(W_PARA); % <0
            k_c = network.network_hwy.(linkStr).paras(K_C_PARA);
            k_m = network.network_hwy.(linkStr).paras(K_M_PARA);
            
            % determine the number of boundary conditions that need to be
            % used to sample the time point
            num_inspSteps = sum(network.T_cum < t); % t \in (t_ref, t_ref+interval]
            
            % save the initial condition in format [x_begin x_end density]
            Con_ini = network.Con_initial( network.Con_initial(:,LINKLABEL) == link,...
                P_MIN: DENS);
            Con_inidM = -(Con_ini(:,2)-Con_ini(:,1)).*Con_ini(:,3);
            Con_iniCumM(1) = 0;
            for i =1:length(Con_inidM)
                Con_iniCumM(i+1) = sum(Con_inidM(1:i));
            end
            Con_iniTotalM = Con_iniCumM( length(Con_iniCumM) ); % <0
            Con_iniCumM( length(Con_iniCumM) ) = [];
            if ~iscolumn(Con_iniCumM)
                Con_iniCumM = Con_iniCumM';
            end
            
            
            % if in domain of initial conditions, update the solutions at
            % those two points
            
            % save the upstream/downstream condition in format [t_begin t_end flow]
            Con_us(:,1) = network.T_cum(1:num_inspSteps);
            Con_us(:,2) = network.T_cum(2:num_inspSteps+1);
            Con_us(num_inspSteps,2) = t;    % update the end time
            Con_ds(:,1) = network.T_cum(1:num_inspSteps);
            Con_ds(:,2) = network.T_cum(2:num_inspSteps+1);
            Con_ds(num_inspSteps,2) = t;
            
            Con_us(:,3) = x( Dec(link,NUM_UP,1):Dec(link,NUM_UP,1) + num_inspSteps - 1 );
            Con_ds(:,3) = x( Dec(link,NUM_DOWN,1):Dec(link,NUM_DOWN,1) + num_inspSteps - 1 );
            
            % get the absolute vehicle label at those step points
            Con_usdM = (Con_us(:,2)-Con_us(:,1)).*Con_us(:,3);
            Con_usCumM(1) = 0;
            for i = 1:length(Con_usdM)
                Con_usCumM(1+i) = sum( Con_usdM(1:i) );
            end
            Con_usCumM( length(Con_usCumM) ) = [];
            if ~iscolumn(Con_usCumM)
                Con_usCumM = Con_usCumM';
            end
            
            Con_dsdM = (Con_ds(:,2)-Con_ds(:,1)).*Con_ds(:,3);
            Con_dsCumM(1) = Con_iniTotalM;
            for i = 1:length(Con_dsdM)
                Con_dsCumM(1+i) = Con_iniTotalM + sum( Con_dsdM(1:i) );
            end
            Con_dsCumM( length(Con_dsCumM) ) = [];
            if ~iscolumn(Con_dsCumM)
                Con_dsCumM = Con_dsCumM';
            end
            
            if strcmp(bound,'up')
                
                if sum(network.T_cum > t_ref & network.T_cum < t) ~= 0
                    % remove the last two pieces of upstream conditions
                    Con_us( max(1,num_inspSteps-1):num_inspSteps , 3) = NaN;
                else
                    % remove the last piece of upstream condition
                    Con_us( num_inspSteps ,3 ) = NaN;
                end
                
            elseif strcmp(bound,'dn')
                
                if sum( network.T_cum > t_ref & network.T_cum < t ) ~= 0
                    % remove the last two pieces of downstream conditions
                    Con_ds( max(1,num_inspSteps-1):num_inspSteps , 3) = NaN;
                else
                    % remove the last piece of downstream condition
                    Con_ds( num_inspSteps , 3) = NaN;
                end
            else
                error('Check bound in samplePointsLink');
            end
            
            
            % points coordinates
            len_link = network.network_hwy.(linkStr).profile(POSTKM)*1000;
            if strcmp(bound,'up')
                position = 0;
            elseif strcmp(bound,'dn')
                position = len_link;
            end
            
            t_array = [ t_ref, t];
            for i = 1:length(t_array)
                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (Con_ini(:,3) <= k_c);
                greaterThanKc = (Con_ini(:,3) > k_c);
                
                %INITIAL====<k_c===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = ( Con_ini(lowerThanKc,1)+v_f*t_array(i) <= position &...
                    Con_ini(lowerThanKc,2)+v_f*t_array(i) >= position);
                
                inFanDomain = ( Con_ini(lowerThanKc,1)+v_f*t_array(i) > position &...
                    Con_ini(lowerThanKc,1)+w*t_array(i) <= position);
                
                % solution in characteristic domain
                tmp_M = min( Con_iniCumM(inCharacterDomain,1) -...
                    Con_ini(inCharacterDomain,3).*(position - v_f*t_array(i) - ...
                    Con_ini(inCharacterDomain,1)) );
                M(i) = minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( Con_iniCumM(inFanDomain,1) -...
                    k_c.*(position - v_f*t_array(i) - ...
                    Con_ini(inFanDomain,1)  )  );
                M(i) = minNonEmpty(M(i), tmp_M);
                
                %INITIAL====>k_c===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = ( Con_ini(greaterThanKc,1)+w*t_array(i) <= position &...
                    Con_ini(greaterThanKc,2)+w*t_array(i) >= position);
                
                inFanDomain = ( Con_ini(greaterThanKc,2)+v_f*t_array(i) > position &...
                    Con_ini(greaterThanKc,2)+w*t_array(i) <= position);
                
                % solution in characteristic domain
                tmp_M = min( Con_iniCumM(inCharacterDomain,1) -...
                    Con_ini(inCharacterDomain,3).*(position - w*t_array(i) - ...
                    Con_ini(inCharacterDomain,1)) -k_m*t_array(i)*w);
                M(i) = minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( Con_iniCumM(inFanDomain,1)+ Con_inidM(inFanDomain,1)...
                    - k_c.*(position - w*t_array(i) - Con_ini(inFanDomain,2)) ...
                    - k_m*t_array(i)*w);
                M(i) = minNonEmpty(M(i), tmp_M);
                
                
                %UPSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(Con_us(:,3));
                inCharacterDomain = ( Con_us(notNaN,1)+ position/v_f <= t_array(i) &...
                    Con_us(notNaN,2)+ position/v_f >= t_array(i));
                
                inFanDomain = ( Con_us(notNaN,2)+ position/v_f < t_array(i) );
                
                % solution in characteristic domain
                tmp_M = min( Con_usCumM(inCharacterDomain,1) +...
                    Con_us(inCharacterDomain,3).*(t_array(i) - position/v_f - ...
                    Con_us(inCharacterDomain,1)) );
                M(i) = minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( Con_usCumM(inFanDomain,1) +...
                    Con_usdM(inFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - position/v_f - ...
                    Con_us(inFanDomain,2)) );
                M(i) = minNonEmpty(M(i), tmp_M);
                
                
                %DOWNSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(Con_ds(:,3));
                inCharacterDomain = ( Con_ds(notNaN,1) + (position - len_link)/w <= t_array(i) &...
                    Con_ds(notNaN,2)+ (position - len_link)/w >= t_array(i));
                
                inFanDomain = ( Con_ds(notNaN,2)+ (position - len_link)/w < t_array(i) );
                
                % solution in characteristic domain
                tmp_M = min( Con_dsCumM(inCharacterDomain,1) +...
                    Con_ds(inCharacterDomain,3).*(t_array(i) -...
                    (position - len_link)/w - ...
                    Con_us(inCharacterDomain,1))...
                    - k_m*(position - len_link));
                
                M(i) = minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( Con_dsCumM(inFanDomain,1) +...
                    Con_dsdM(inFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - (position-len_link)/v_f - ...
                    Con_us(inFanDomain,2)) );
                
                M(i) = minNonEmpty(M(i), tmp_M);
                
                
            end % end for each t (here we just have (t_ref and t) )
            
            if all(~isnan(M)) && all(~isempty(M))
                
                % use car ID at t_ref as a reference
                d_M = M(2) - M(1);
                
            else
                error('failed to sample points');
            end
            
            
        end % end samplePointsLink function
        
        
    end     %end methods
    
end






















