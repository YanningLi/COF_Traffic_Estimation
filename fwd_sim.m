% This script is a simple forward simulation for debugging the COF toolbox. 
% For other applications, write a script file based on this one to define
% the optimization program.
% Yanning, Sep 01, 2015

% TODO:
%   1. print computation time. Auto log output.

% Remark: 
%   1. This version does not handle ramp junctions yet.


tic
clearvars -except dbg

profile on
%===============================================================
% configuration paramters
runCtrl = 'no_control'; % This parameter can be later used for optimal control 
entropyTolerance = 1; % number of vehicles that we would consider as entropic solution

sim_steps = 7;
step_length = 30;   % seconds
T_durations = ones(sim_steps,1)*step_length;
T_BC_cum = [0; cumsum(T_durations)];


start_time = 0;
end_time = sum(T_durations);

%===============================================================
% This is the time discretization grid. Will be iteratively updated to
% eliminate the discretization error.
T_junc= T_durations;

% initialize environment, parameters, index.
initEnv;    

%===============================================================
if strcmp(runCtrl,'no_control') == 1
    % here we define an example
    % A merge with three links: link 1, 2 merge to link 3
    
    ratio = 1.5817; % split ratio for the merge. Need to be calibrated.

    % test length
    len_link1 = 1.2;    %km
    len_link2 = 1.1;
    len_link3 = 1.0;
    
    % Initial traffic density is initialized as segments with equal length
    % Normalized to rho_c
    Ini_1 = [1 1 1 1];
    Ini_2 = [3 0 3 3];
    Ini_3 = [3 0 3 3];
    
    % Boundary condition at the two entrances and the one exit
    % Normalized to q_max
    q1_us_data = ones(sim_steps,1);
    q2_us_data = ones(sim_steps,1);
    q3_ds_data = ones(sim_steps,1);
    
    % adjust flows for feaibility and different scenarios
    q1_us_data(10:sim_steps,1) = 0;
    q2_us_data([3 4 5],1) = 0.2;
    q3_ds_data(1:sim_steps,1) = 0;
    
elseif strcmp(runCtrl, 'corsim') == 1
    % Initial density are all 0
    % total simulation time 6000s, 0-10min:0; 10-90min:large; 90-100min: 0
    % 10 segments for initial conditions
    Ini_1 = zeros(2,1);
    Ini_2 = zeros(2,1);
    Ini_3 = zeros(5,1);
    % add 1 car in the 4th and 5th segment to avoid infeasibility
    Ini_3(4:5) = (1/(len_link3*1000/5))/default_para.k_c_pl;
    
    % here qus and down are supposed to read from files
    data = dlmread('./data/dataLog_cfg5.csv',';');
    
    %normalized flow veh/s to the maximal flow on each link
    q1_us_data = (data(data(:,1)==1,2)/30)/default_para.q_max_pl;
    q2_us_data = (data(data(:,1)==4,2)/30)/default_para.q_max_pl;
    q3_ds_data = (data(data(:,1)==2,3)/30)/default_para.q_max_pl;
    
    % only use the sim_steps data
    q1_us_data = q1_us_data(1:sim_steps);
    q2_us_data = q2_us_data(1:sim_steps);
    q3_ds_data = q3_ds_data(1:sim_steps);
    
    % Due to the fluctuation of the flow data and we take the mean max as the
    % maximal flow, some boundary values here are greater than 1.
    % Truncate them to 1 to avoid infeasibility.
    q1_us_data(q1_us_data > 1) = 1;
    q2_us_data(q2_us_data > 1) = 1;
    q3_ds_data(q3_ds_data > 1) = 1;
end

% represent the boundary condition with the time grid information
q1_us_TF = [T_BC_cum(1:sim_steps), T_BC_cum(2:sim_steps+1), q1_us_data];
q2_us_TF = [T_BC_cum(1:sim_steps), T_BC_cum(2:sim_steps+1), q2_us_data];
q3_ds_TF = [T_BC_cum(1:sim_steps), T_BC_cum(2:sim_steps+1), q3_ds_data];


% This is the while loop for iteratively regridding and eventually getting
% the entropy condition
% Flags for iteratively updating time discretization
getEntropy = false;
loopCounter = 0;
while getEntropy == false && loopCounter <=10
    
    loopCounter = loopCounter+1;
    
    % Keep track of the objective funciton values, the solution, and the
    % split ratio
    f = [];
    Sol = [];
    R1 = [];
    R2 = [];
    
    %===============================================================
    % Define the network
    % initLink(obj, link, para, num_lanes, lengthKM, linkType)
    net = initNetwork;
    net.addLink(1, default_para, 1, len_link1, 'freeway');
    net.addLink(2, default_para, 1, len_link2, 'freeway');
    net.addLink(3, default_para, 1, len_link2, 'freeway');
    net.addJunc(1, [1; 2], 3, 'merge', [ratio; 1], T_durations);
    
%     %===============================================================
%     % Regrid the boundary data
%     q1_us_data = regridBoundaryCondition(q1_us_TF,T_new_cum);
%     q2_us_data = regridBoundaryCondition(q2_us_TF,T_new_cum);
%     q3_ds_data = regridBoundaryCondition(q3_ds_TF,T_new_cum);
    
    %===============================================================                         
    % set initial conditoins
    net.setInitialCon(1,Ini_1);
    net.setInitialCon(2,Ini_2);
    net.setInitialCon(3,Ini_3);
    
    % setBoundaryCon(obj, link, q_in, q_out, T_in, T_out)
    net.setBoundaryCon(1, q1_us_data, [], T_durations, T_junc);
    net.setBoundaryCon(2, q2_us_data, [], T_durations, T_junc);
    net.setBoundaryCon(3, [], q3_ds_data, T_junc, T_durations);

    %===============================================================
    % define and solve optimization program
    LP = optProgram;
    LP.setConfig(start_time, end_time);
    LP.setConstraints(net, default_para.e_max);
    
    %===============================================================
    % Add objective functions
    LP.addEntropy(net, 1);
    %LP.addOnRampControl(net,[1]);
    %LP.maxUpflow([1 2 4]);
    %LP.maxDownflow(3);
    %===============================================================    
    toc
    
    tic
    %solve program
    [x, fval, exitflag, output] = LP.solveProgram;
    
    %===============================================================
    % Post process the data
    Mos = postSolution(x, net, LP.dv_index, LP.end_time, 2, 2);
    Mos.plotJuncs(net, 'all')
    
    [getEntropy, steps] = Mos.checkEntropy(entropyTolerance);
    
    if getEntropy == false
        hold on
        T_junc_cum = Mos.updateTimeDiscretization(steps);
    end

    % updated T_junc
    % a quick hack, only works for this particular one junc scenario
    T_junc = T_junc_cum.junc_1(2:end)-T_junc_cum.junc_1(1:end-1);
    
end



toc


% profile viewer
% p = profile('info');
% profsave(p,'./result/matlab_profiler_results');





