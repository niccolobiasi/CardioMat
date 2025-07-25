function [Vm_array, state_array]=bench_model(model,Ie,dt,domain,scalar)
% bench_model is thought to test an ionic model Matlab function before
% using it into CardioMat. 
% bench_model(model) executes a single cell simulation by using the ionic
% model function specified by model string. The default simulation is
% 500 ms long and consists of a single stimulation at time t=50 ms with
% current= 40 mV/ms. Optionally, the user can input an Ie array whose
% length defines the simulation duration and the entry value the
% stimulation current for each time instant.
% Other optional inputs are the temporal resolution (dt), the "domain" and
% scalar parameters to be used in the model function. 
% The function returns the computed membrane potential and states
% in time as numeric arrays.



if nargin<5 || isempty(scalar)
    scalar=1;
end

if nargin<4 || isempty(domain)
    domain=1;
end

if nargin<3 || isempty(dt)
    dt=0.02;
end

if nargin<2 || isempty(Ie)
    Nt=round(500/dt)+1;
    Ie=zeros(Nt,1);
    tt=round(50/dt)+1;
    Npulses=1;
    dT=round(1/dt);
    BCL=1000;
    for k=1:Npulses
        Ie(tt:tt+dT)=40;
        tt=tt+round(BCL/dt);
    end
else
    Nt=length(Ie);
end

Vm_array=NaN(size(Ie));

[Vm,state]=feval(model);
Vm_array(1)=Vm;
state_array=NaN(Nt, length(state));
state_array(1,:)=cell2mat(state);

for i=1:Nt-1
    [Vm,state]=feval(model,Vm,state,Ie(i),0,domain,dt,scalar);
    Vm_array(i+1)=Vm;
    state_array(i+1,:)=cell2mat(state);
end




