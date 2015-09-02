% initEnv
% initialize the simulation environment, including parameters, index


%% Path
clear global path path_data path_profiles path_result;
global path path_data path_profiles path_result;
path = [pwd '/'];   
path_data = [path 'data/'];    % place to input data
path_profiles = [path 'profiles/'];
path_result = [path 'result/'];
if ~exist(path_profiles, 'dir')
    mkdir(path_profiles);
end
if ~exist(path_result, 'dir')
    mkdir(path_result);
end

%% Global index
clear global POSTKM LINKTYPE; 
clear global V_PARA W_PARA K_C_PARA K_M_PARA Q_MAX_PARA V_MIN_PARA V_MAX_PARA LAST_PARA; 
global POSTKM LINKTYPE; 
global V_PARA W_PARA K_C_PARA K_M_PARA Q_MAX_PARA V_MIN_PARA V_MAX_PARA LAST_PARA; 

% Declare global constants Dec.hwy_#.counts index
clear global NUM_UP NUM_DOWN NUM_INI NUM_LAST;
global NUM_UP NUM_DOWN NUM_INI NUM_LAST;

% Declare global constants for value condition index
clear global LINKLABEL TIME P_MIN P_MAX DENS FLUX;
global LINKLABEL TIME P_MIN P_MAX DENS FLUX;

%% network_hwy index
% network_hwy.hwy_#.topoPoints
% LON = 1; LAT = 2; POSTKM = 3;

% Index for network_hwy.link_#; ramps 
POSTKM = 1;
LINKTYPE = 2;   %0-mainline; 1-onramp; only used for onramp control

% Index for network_hwy.hwy_#; paras
V_PARA = 1;
W_PARA = 2;
K_C_PARA = 3;
K_M_PARA = 4;
Q_MAX_PARA = 5;
V_MIN_PARA = 6;
V_MAX_PARA = 7;
LAST_PARA = V_MAX_PARA;

%% value conditions index
% initial
LINKLABEL = 1;
TIME = 2; 
P_MIN = 3;  
P_MAX = 4;
DENS = 5;

% boundary
% HWYLABEL = 1; LINKLABEL = 2; T_MIN = 3; 
FLUX = 3;


%% Dec 3d matrix containing index of variabels in the decision variable
% Decision variable locator, s atruct
        % Dec.(linkStr) = [start of upstream, start of downstream, start of initial;
        %                  end of upstream, end of downstream, end of initial]
INDEX_UP = 1;
INDEX_DOWN = 2;
INDEX_INI = 3;
INDEX_LAST = NUM_INI;


%% DecBool matrix for index of boolean variables modeling meters 
% num_link x 2 matrix [start of index, end of index]
% [101, 105]    % link 1 has onramp, define 5(tail_steps) boolean variable
% [105, 105]    % not onramp
% [106, 110]    % on ramp
% [110, 110]    % total length of decision variable


%% default_para
default_para = struct;
default_para.beta_off = 0.2;
default_para.v = 65*1609/3600;  %65 miles/hr
default_para.w = -(7.5349)*1609/3600;     %m/s calibrated from corsim
%default_para.w = -(8.125)*1609/3600;     %m/s Here v/w = 8, just to resolve the discritization error
default_para.k_c_pl = 34.6569/1609;          %veh/m calibrated from corsim
default_para.q_max_pl = default_para.k_c_pl*default_para.v;  %veh/s
default_para.q_on_max = default_para.k_c_pl*default_para.v; 
default_para.q_off_max = default_para.k_c_pl*default_para.v;
default_para.k_m_pl = default_para.k_c_pl*(default_para.w-default_para.v)/default_para.w;
default_para.v_min = 0*1609/3600;
default_para.v_max = 65*1609/3600;
default_para.e_max = 0.0;









