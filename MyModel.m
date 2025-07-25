function [Vm,state]=MyModel(Vm,state,Ie,dv2dt,domain,dt,scalar)
% This function serves as a template for the definition of ionic model
% function to be used with CardioMat.

if nargin==0
    % MODEL DEFAULT INITIAL CONDITIONS
    % Vm=....;
    % state_1=...;
    % state_2=...;
    % ............
    % state_N=....
else
    % Run model function 
    % [Vm,state_1, state_2,...., state_N]=arrayfun(@MyModel_step,Vm,state{1},state{2},........,state{N},Ie,dv2dt,domain,scalar);
end

% state{1}=state_1;
% state{2}=state_2;
% .....
% state{N}=state_N;

      %% MODEL FUNCTION
    % function [Vm,state_1,state_2,....,state_N]=MyModel_step(Vm,state_1,state_2,....,state_N,Ie,dv2dt,region,scalar)
    % 
    %     .........
    %     .........
    %     .........
    % 
    % 
    %     %update membrane potential
    %     Vm=Vm+dv2dt*dt;
    %     
    %     % update gate variables
    %     state_1=state_1+dstate_12dt*dt; OR use Rush-Larsen (recomended
    %     for most models)
    %     state_1=state_1_inf+(state_1-state_1_inf)*exp(-dt/state_1_tau);
    % end

end