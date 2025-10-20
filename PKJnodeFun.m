function [s1,s2,state,t_last]=PKJnodeFun(A1,A2,t,state,t_last,delay1,delay2,d1,d2,ERP)
% PKJnodeFun implements a finite state machine for the Purkinje muscular
% junctions. A1 is anterograde activation and A2 retrograde
% activation. State=1 indicates anterograde conduction, state=2 indicates
% retrogarde conduction, state=3 is the refractory period, state=0 is the
% rest condition. delay1 and delay2 are anterograde and retrograde delays.
% d1, d2 are duration of stimulation in the two ways. ERP is
% the effective refractory period. See runSimulation for details about
% parameter values. PKJnodeFun returns the stimulus flag s1 and s2, the
% actual state and the time of the last activation (which also serves as
% input to the function).
if isscalar(delay1)
    delay1=delay1*ones(size(A1));
end

if isscalar(delay2)
    delay1=delay2*ones(size(A1));
end

if isscalar(d1)
    delay1=d1*ones(size(A1));
end


if isscalar(d2)
    delay1=d2*ones(size(A1));
end

if isscalar(ERP)
    ERP=ERP*ones(size(A1));
end

[s1,s2,state,t_last] = arrayfun(@fun,A1,A2,state,t_last,delay1,delay2,d1,d2,ERP);


    function  [s1,s2,state,t_last]=fun(A1,A2,state,t_last,delay1,delay2,d1,d2,ERP)

        %ERP
        if(state == 3 && (t-t_last) > ERP)

            state = 0;

        end

        if (state == 0 && A1)

            state = 1;
            t_last = t;
            
        elseif (state == 0 && A2)

            state = 2;
            t_last = t;
        end

        if(state == 1 && (t-t_last) >= delay1)

            s1 = true;
            s2 = false;

            if((t-t_last-delay1) >= d1)

                state = 3;
            end

        elseif(state == 2 && (t-t_last) >= delay2)

            s1 = false;
            s2 = true;

            if((t-t_last-delay2) >= d2)

                state = 3;
            end

        else

            s1 = false;
            s2 = false;

        end

    end
end
