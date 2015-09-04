% Yanning Li, Sep 01, 2015

% This class sets the equality constraints for the optimization program.


classdef setEqConstraints
    
    properties
        
       EqMatrix;  
        
    end
    
    
    methods
        
        function si = setEqConstraints(network, dev_index, dev_index_max)
            global INDEX_UP INDEX_DOWN
            
            si.EqMatrix = zeros(0,dev_index_max);
            
            for junc = network.junc_labels'
                
                juncStr = sprintf('junc_%d',junc);
                
                num_steps = length(network.network_junc.(juncStr).T);
                
                array = zeros(num_steps, dev_index_max );
                
                if strcmp(network.network_junc.(juncStr).type_junc,'diverge')
                    
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).inlabel);
                    array(1:num_steps,...
                          dev_index.(linkStr)(1, INDEX_DOWN):dev_index.(linkStr)(2, INDEX_DOWN)) =...
                          diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).outlabel(1));
                    array(1:num_steps,...
                          dev_index.(linkStr)(1, INDEX_UP):dev_index.(linkStr)(2, INDEX_UP)) =...
                          -diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).outlabel(2));
                    array(1:num_steps,...
                          dev_index.(linkStr)(1, INDEX_UP):dev_index.(linkStr)(2, INDEX_UP)) =...
                          -diag(ones(num_steps,1));

                    si.EqMatrix = [si.EqMatrix; array];
                    
                elseif strcmp(network.network_junc.(juncStr).type_junc,'merge')
                    
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).inlabel(1));
                    array(1:num_steps,...
                          dev_index.(linkStr)(1, INDEX_DOWN):dev_index.(linkStr)(2, INDEX_DOWN)) =...
                          diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).inlabel(2));
                    array(1:num_steps,...
                          dev_index.(linkStr)(1, INDEX_DOWN):dev_index.(linkStr)(2, INDEX_DOWN)) =...
                          diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).outlabel);
                    array(1:num_steps,...
                          dev_index.(linkStr)(1, INDEX_UP):dev_index.(linkStr)(2, INDEX_UP)) =...
                          -diag(ones(num_steps,1));
 
                    si.EqMatrix = [si.EqMatrix; array];
  
                elseif strcmp(network.network_junc.(juncStr).type_junc,'connection')
                    
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).inlabel);
                    array(1:num_steps,...
                          dev_index.(linkStr)(1, INDEX_DOWN):dev_index.(linkStr)(2, INDEX_DOWN)) =...
                          diag(ones(num_steps,1));
                      
                    linkStr = sprintf('link_%d',network.network_junc.(juncStr).outlabel);
                    array(1:num_steps,...
                          dev_index.(linkStr)(1, INDEX_UP):dev_index.(linkStr)(2, INDEX_UP)) =...
                          -diag(ones(num_steps,1));
                                          
                    si.EqMatrix = [si.EqMatrix; array];
                    
                end
                
            end
            
        end
        
        
        
    end
    
end
    
    
