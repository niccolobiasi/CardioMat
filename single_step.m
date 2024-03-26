function [Vm_array, u_array, w_array]=single_step(Vm,u,w,Ie,dv2dt,domain,dt,scalar)
% This function defines the ionic model to use, that in our case is the 
% Biasi and Tognetti model. Currently, the possibility of using other models
% is not implemented, but it will be implemented soon. 
% Vm,u, and w are the state variable of the model. Ie is the external
% current and dv2dt is the diffusive current. domain defines which
% paramater set should be used:
%1->epicardial
%2->endocardial
%3->midcardial
%4->atrial
%5->Brugada coved type AP
%6->Brugada lost dome type AP
% any other -> epicardial
%
%dt is the temporal resolution. teh scalar input represent the normalized
%pixelintensity that could be used to modulate EP properties of the model.
% The function runs on GPU.

[Vm_array, u_array,w_array]=arrayfun(@single_step_par,Vm,u,w,Ie,dv2dt,domain,scalar);


    function [Vm,u,w]=single_step_par(Vm,u,w,Ie,dv2dt,region,scalar)
       
        if isnan(Vm)
            Vm = NaN;
            u  = NaN;
            w  = NaN;
            return
        end

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
        fact=1-scalar/2;
        fAct=1-scalar/4;
        A0=A0*fAct;
        A1=A1*fAct;
        e110=e110*fact;
        e120=e120*fact;
        e21=e21*fact;
        e22=e22*fact;
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
        dv2dt=dv2dt-k*g.*(Iexc+Irec+Ito)+Ie;
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
