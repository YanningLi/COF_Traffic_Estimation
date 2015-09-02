% Yanning Li, Sep 01, 2015

% This class sets the equality constraints for the optimization program.


classdef setEqConstraints
    
    properties
        
       EqMatrix;  
        
    end
    
    
    methods
        
        function si = setEqConstraints_grid(network,Dec)
            global NUM_UP NUM_DOWN
            
            si.EqMatrix = zeros(0,Dec(size(Dec,1),size(Dec,2),2));
            
            for junc = 1:network.num_junc
                
                juncStr = sprintf('junc_%d',junc);
                num_steps = length(network.network_junc.(juncStr).T);
                
                array = zeros(num_steps,Dec(size(Dec,1),size(Dec,2),2));
                
                if strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                    
                    array(1:num_steps,Dec(network.network_junc.(juncStr).inlabel,NUM_DOWN,1):...
                        Dec(network.network_junc.(juncStr).inlabel,NUM_DOWN,2)) = diag(ones(num_steps,1));
                    array(1:num_steps,Dec(network.network_junc.(juncStr).outlabel(1),NUM_UP,1):...
                        Dec(network.network_junc.(juncStr).outlabel(1),NUM_UP,2)) = -diag(ones(num_steps,1));
                    array(1:num_steps,Dec(network.network_junc.(juncStr).outlabel(2),NUM_UP,1):...
                        Dec(network.network_junc.(juncStr).outlabel(2),NUM_UP,2)) = -diag(ones(num_steps,1));
                    
                    si.EqMatrix = [si.EqMatrix; array];
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'merge')
                    
                    array(1:num_steps,Dec(network.network_junc.(juncStr).inlabel(1),NUM_DOWN,1):...
                        Dec(network.network_junc.(juncStr).inlabel(1),NUM_DOWN,2)) = diag(ones(num_steps,1));
                    array(1:num_steps,Dec(network.network_junc.(juncStr).inlabel(2),NUM_DOWN,1):...
                        Dec(network.network_junc.(juncStr).inlabel(2),NUM_DOWN,2)) = diag(ones(num_steps,1));
                    array(1:num_steps,Dec(network.network_junc.(juncStr).outlabel,NUM_UP,1):...
                        Dec(network.network_junc.(juncStr).outlabel,NUM_UP,2)) = -diag(ones(num_steps,1));
                    
                    si.EqMatrix = [si.EqMatrix; array];
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'onramp')
                    
                    array(1:num_steps,Dec(network.network_junc.(juncStr).inlabel(1),NUM_DOWN,1):...
                        Dec(network.network_junc.(juncStr).inlabel(1),NUM_DOWN,2)) = diag(ones(num_steps,1));
                    array(1:num_steps,Dec(network.network_junc.(juncStr).inlabel(2),NUM_DOWN,1):...
                        Dec(network.network_junc.(juncStr).inlabel(2),NUM_DOWN,2)) = diag(ones(num_steps,1));
                    array(1:num_steps,Dec(network.network_junc.(juncStr).outlabel,NUM_UP,1):...
                        Dec(network.network_junc.(juncStr).outlabel,NUM_UP,2)) = -diag(ones(num_steps,1));
                    
                    si.EqMatrix = [si.EqMatrix; array];
                
                
                elseif strcmp(network.network_junc.(juncStr).type_junc,'connection')
                    
                    array(1:num_steps,Dec(network.network_junc.(juncStr).inlabel,NUM_DOWN,1):...
                        Dec(network.network_junc.(juncStr).inlabel,NUM_DOWN,2)) = diag(ones(num_steps,1));
                    array(1:num_steps,Dec(network.network_junc.(juncStr).outlabel,NUM_UP,1):...
                        Dec(network.network_junc.(juncStr).outlabel,NUM_UP,2)) = -diag(ones(num_steps,1));
                    
                    si.EqMatrix = [si.EqMatrix; array];
                    
                end
                
            end
            
        end
        
        
        
    end
    
end
    
    
