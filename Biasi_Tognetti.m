function [Vm,state,dv2dt]=Biasi_Tognetti(Vm,state,Ie,dv2dt,domain,dt,scalar)
% This function implements the Biasi-Tognetti model for the CardioMat 
% toolbox. 7 different parameters sets are reported below:
%   -Epicardium (domain=1)
%   -Endocardium (domain=2)
%   -Midcardium (domain=3)
%   -Atria (domain=4)
%   -Brugada coved-type (domain=5)
%   -Brugada saddleback-type (domain=6)
%
% The parameters for each region can be modified directly inside the
% function. Additional parameter sets associated to different domain
% numbers can be added. Further details on this model are given in:
% https://doi.org/10.1371/journal.pone.0259066
% https://doi.org/10.1109/ACCESS.2022.3222830
% https://doi.org/10.1038/s41598-022-12239-9


if nargin==0
    Vm=-85; %if called without arguments returns the default inital conditions
    u=0;
    w=0;
else
    [Vm,u,w,dv2dt]=arrayfun(@Biasi_Tognetti_step,Vm,state{1},state{2},Ie,dv2dt,domain,scalar);
end

state{1}=u;
state{2}=w;

    function [Vm,u,w,dv2dt]=Biasi_Tognetti_step(Vm,u,w,Ie,dv2dt,region,scalar)
        
        %MODEL PARAMETERS
        if region ==2
            %ENDOCARDIUM
            A0 = 125;
            A1 = 60;
            e110 = 0.001;
            e120 = 0.0133;
            e21 = 0.0005;
            e22 = 0.0331;
            gamma = 9.25;
            gamma1 = 23.125;
            ew0 = 0.025;
            dw0 = 5;
            k = 1;
            c1 = 2.6;
            c2 = 1;
            c3 = 0.5;
            a = 0.18;
            B=-85; %[mV]
            alpha = 15;
            thetau = 0.2;
            g0 = 0.1;
            um = 0.58;

        elseif region ==3
            %MIDCARDIUM
            A0 = 134;
            A1 = 0;
            e110 = 6.03e-5;
            e120 = 0.013118;
            e21 = 0.00093397;
            e22 = 0.03;
            gamma = 11.875;
            gamma1 = 28.91;
            ew0 = 0.03;
            dw0 = 3;
            k = 1;
            c1 = 2.6;
            c2 = 1;
            c3 = 0.5;
            a = 0.18;
            B=-85; %[mV]
            alpha = 15;
            thetau = 0.2;
            g0 = 0.1;
            um = 1;

        elseif region ==4
            %ATRIA
            A0 = 100;
            A1 = 0;
            e110 = 0.0015;
            e120 = 0.0272;
            e21 = 0.0118;
            e22 = 0.0132;
            gamma = 6;
            gamma1 = 33;
            ew0 = 0.04;
            dw0 = 0.73;
            k = 1;
            c1 = 2.6;
            c2 = 1;
            c3 = 0.5;
            a = 0.18;
            B=-75; %[mV]
            alpha = 25;
            thetau = 0.1;
            g0 = 0.07;
            um = 0.8;

        elseif region==5
            %BRUGADA COVED
            A0 = 90;
            A1 = 500;
            e110 = 0.0059;
            e120 = 0;
            e21 = 0.015;
            e22 = 0;
            gamma = 3;
            gamma1 = 7.5;
            ew0 = 0.06;
            dw0 = 0.6;
            k = 1;
            c1 = 2.6;
            c2 = 1;
            c3 = 0.5;
            a = 0.18;
            B=-85; %[mV]
            alpha = 15;
            thetau = 0.2;
            g0 = 0.1;
            um = 0.58;

        elseif region==6
            %Brugada LD
            A0 = 90;
            A1 = 500;
            e110 = 0.0059;
            e120 = 0;
            e21 = 0.015;
            e22 = 0;
            gamma = 3;
            gamma1 = 7.5;
            ew0 = 0.06;
            dw0 = 0.35;
            k = 1;
            c1 = 2.6;
            c2 = 1;
            c3 = 0.5;
            a = 0.18;
            B=-85; %[mV]
            alpha = 15;
            thetau = 0.2;
            g0 = 0.1;
            um = 0.58;

        % elseif region==7
        %     NEW PARAMETER SET
        % 
        % 

        else
            %EPICARDIUM
            A0 = 135;
            A1 = 0;
            e110 = 0.0059;
            e120 = 0;
            e21 = 0.015;
            e22 = 0;
            gamma = 8;
            gamma1 = 20;
            ew0 = 0.04;
            dw0 = 0.6;
            k = 1;
            c1 = 2.6;
            c2 = 1;
            c3 = 0.5;
            a = 0.18;
            B=-85; %[mV]
            alpha = 15;
            thetau = 0.2;
            g0 = 0.1;
            um = 0.58;

        end
        fact=1-scalar/4;
        A0=A0*fact;
        A1=A1*fact;
        e110=e110*fact;
        e120=e120*fact;
        ew0=ew0*fact;
        dw0=dw0/fact;
        gamma=gamma*fact;
        gamma1=gamma1*fact;




        %compute parameters
        e0 = e110 + e120 *abs(u);
        e1 = e21 + e22*abs(u);
        A= A0+ A1*u^2;
        g=(gamma +gamma1*u)*(-tanh(alpha *(u-thetau))+1)/2+g0;
        vr=(Vm-B);
        v=vr/A;
        su=(u-um)^2/um^2;

        %Compute currents
        Iexc=c1*vr*(1-v)*(a-v);
        Irec=c2*u*vr;
        Ito=c3*w*vr*su;

        %compute derivatives
        dv2dt=dv2dt-k*g*(Iexc+Irec+Ito)+Ie;
        dw=dw0/su;
        if dw>1000
            dw=1000;
        end
        dw2dt=k*ew0*g*(v-dw*w);
        du2dt=v-u;

        %compute e
        if du2dt>=0
            if (dv2dt>=0 && dw2dt>=0)
                e=e0*g;
            else
                e=e0;
            end
        else
            e=e1;
        end


        %update variables
        Vm=Vm+dv2dt*dt;
        w=w+dw2dt*dt;
        u=u+k*e*du2dt*dt;
    end
end