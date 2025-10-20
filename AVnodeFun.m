function [sA, sV, state, t_last] = AVnodeFun(AE,VE,t,state,t_last,delayA,delayR,stim_d1,stim_d2,ERP)
% AVnodeFun implements the same finite state machine of PKJnodeFun for just
% and prints in the command window AVN activation times.

    if(state == 3 && (t-t_last) > ERP)
        state = 0; %REST
    end

    if (state == 0 && AE) 
        
        state = 1;  %anterograde conduction
        t_last = t;
        disp(strcat("AVN activated on the atrial side at: ",num2str(t)));
    elseif (state == 0 && VE)
        
        state = 2;  %retrograde conduction
        t_last = t;
        disp(strcat("AVN activated on the ventricular side at: ",num2str(t)));
    end
    
    if(state == 1 && (t-t_last) >= delayA) %anterograde stimulation

        sA = false;
        sV = true;

        if((t-t_last-delayA) >= stim_d1)
            
            state = 3; %ERP
        end

    elseif(state == 2 && (t-t_last) >= delayR) %retrograde stimualtion
    
        sA = true;
        sV = false;

        if((t-t_last-delayR) >= stim_d2)
            
            state = 3; %ERP
        end
    else
        
        sA = false;
        sV = false;
    end



end