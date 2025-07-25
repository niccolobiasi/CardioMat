function [Vm,state]=CRN_1998(Vm,state,Ie,dv2dt,domain,dt,scalar)
% This function implements the Courtemanche-Ramirez-Nattel (CRN) model  
% (https://doi.org/10.1152/ajpheart.1998.275.1.H301) for the CardioMat 
% toolbox. Normal (domain=1) and chronic AF (domain=2) are included 
% (chronic AF parameters were taken from 
% https://doi.org/10.1016/s0008-6363(99)00034-6).

if nargin==0
    % MODEL DEFAULT INITIAL CONDITIONS
    Vm=-81.2;
    h=0.965;
    d=1.37e-4;
    xr=3.29e-5;
    Na_i=11.2;
    K_i=1.39e2;
    Ca_rel=1.49;
    oi=0.999;
    ui=0.999;
    Ca_i=1.02e-4;
    v=1;
    m=2.91e-3;
    j=9.78e-1;
    f=0.999;
    xs=1.87e-2;
    Ca_up=1.49;
    oa=3.04e-2;
    ua=4.96e-3;
    fCa=7.75e-1;
    u=0;
    w=0.999;
else
    % Run model function
    [Vm,h,d,xr,Na_i,K_i,Ca_rel,oi,ui,Ca_i,v,m,j,f,xs,Ca_up,oa,ua,fCa,u,w]=...
        arrayfun(@CRN_step,Vm,state{1},state{2},state{3},state{4},state{5},...
        state{6},state{7},state{8},state{9},state{10},state{11},state{12},state{13}...
        ,state{14},state{15},state{16},state{17},state{18},state{19},...
        state{20},Ie,dv2dt,domain,scalar);
end


state{1}=h;
state{2}=d;
state{3}=xr;
state{4}=Na_i;
state{5}=K_i;
state{6}=Ca_rel;
state{7}=oi;
state{8}=ui;
state{9}=Ca_i;
state{10}=v;
state{11}=m;
state{12}=j;
state{13}=f;
state{14}=xs;
state{15}=Ca_up;
state{16}=oa;
state{17}=ua;
state{18}=fCa;
state{19}=u;
state{20}=w;


%% MODEL FUNCTION
    function [Vm,h,d,xr,Na_i,K_i,Ca_rel,oi,ui,Ca_i,v,m,j,f,xs,Ca_up,oa,ua,fCa,u,w]=...
            CRN_step(Vm,h,d,xr,Na_i,K_i,Ca_rel,oi,ui,Ca_i,v,m,j,f,xs,Ca_up,oa,...
            ua,fCa,u,w,Ie,dv2dt,region,scalar)

        % CONSTANTS
        R = 8.3143;% [J/(K*mol)] gas constant
        T = 310; % [K] temperature
        F = 96.4867; % [C/mmol] Faraday constant
        Cm = 100; % [pF] membrane capacitance

        V_i = 13668;% [um^3] intrecellular volume
        V_up = 1109.52;% [um^3] SR uptake compartment volume
        V_rel = 96.48;% [um^3] SR release compartment volume

        K_o = 5.4; % [mM] extracellular potassium concentration
        Na_o=140;% [mM] extracellular sodium concentration
        Ca_o = 1.8;% [mM] extracellular calcium concentration


        GNa = 7.8; %[nS/pF] maximal sodium conductance
        G_K1 = 0.09;%[nS/pF] maximal K1 conductance
        G_to = 0.1652;%[nS/pF] maximal transient outward conductance
        G_Kr = 0.0294;%[nS/pF] maximal Kr conductance
        G_Ks = 0.129;%[nS/pF] maximal Ks conductance
        G_Ca_L = 0.1238;%[nS/pF] maximal Ca,L conductance
        g_B_Ca = 0.00113;%[nS/pF] maximal b,Ca conductance
        g_B_Na = 0.000674;%[nS/pF] maximal b,Na conductance
        if region==2
            G_to=0.5*G_to;
            G_Ca_L=G_Ca_L*0.3;
        end
        i_NaK_max = 0.6; %[pA/pF] mamximal NaK current
        I_NaCa_max = 1600;% [pA/pF] mamximal NaCa current
        I_CaP_max = 0.275;%[pA/pF] mamximal p,Ca current
        I_up_max = 0.005;%[mM/ms] mamximal reuptake current


        K_Q10 = 3; %scaling factor for temperature effects

        gamma = 0.35; %voltage dependence parameterfor INaCa

        Km_Na_i = 10;%[mM] Na_i half-saturation constant for INaK
        Km_K_o = 1.5;%[mM] K_o half-saturation constant for INaK
        K_mNa = 87.5;%[mM] Na_o half-saturation constant for INaCa
        K_mCa = 1.38;%[mM] Ca_o half-saturation constant for INaCa
        K_sat = 0.1;% Saturation factor for INaCa
        K_rel = 30;% [ms^-1] Maximal release rate for Irel
        K_up_max = 0.00092;%[mM] Ca_i half saturation constant for Iup

        Ca_up_max = 15;% [mM] maximal Ca concentration in uptake compartment
        CMDN_max = 0.05;% [mM] total calmodulin concentration in myoplasm
        TRPN_max = 0.07;% [mM] total troponin concentration in myoplasm
        CSQN_max = 10;% [mM] total calsequestrin concentration in SR release compartment
        Km_CMDN = 0.00238;%[mM] Ca_i half-saturation constant for calmodulin
        Km_TRPN = 0.0005;% [mM] Ca_i half-saturation constant for troponin
        Km_CSQN = 0.8;% [mM] Ca_rel half-saturation constant for Iup



        tau_tr = 180; %[ms] time constant Ca transfer
        tau_fCa = 2; %[ms] time constant fCa
        sigma =  (1/7)*(exp(Na_o/67.3)-1);
        tau_u = 8;%[ms] time constant gate variable u


        % Nerst Potentials
        E_K= R*T*log(K_o/K_i)/F;
        E_Na= R*T*log(Na_o/Na_i)/F;
        E_Ca= R*T*log(Ca_o/Ca_i)/(2*F);


        %Currents computation
        INa = Cm*GNa*m^3*h*j*(Vm-E_Na);
        IK1 = Cm*G_K1*(Vm-E_K)/(1+exp(0.07*(Vm+80)));
        Ito = Cm*G_to*oa^3*oi*(Vm-E_K);
        GKur = 0.005 + 0.05/(1+exp((-Vm+15)/13));
        if region==2
            GKur=GKur*0.5;
        end
        IKur = Cm*GKur*ua^3*ui*(Vm-E_K);
        IKr = Cm*G_Kr*xr*(Vm-E_K)/(1+exp((Vm+15)/22.4));
        IKs = Cm*G_Ks*xs^2*(Vm-E_K);
        ICaL = Cm*G_Ca_L*d*f*fCa*(Vm-65);
        fNaK=(1+0.1245*exp(-0.1*F*Vm/(R*T))+0.0365*sigma*exp(-F*Vm/(R*T)))^(-1);
        INaK = Cm*i_NaK_max*fNaK*K_o/((1+(Km_Na_i/Na_i)^1.5)*(K_o+Km_K_o));
        INaCa = Cm*I_NaCa_max*(exp(gamma*Vm*F/(R*T))*Na_i^3*Ca_o-exp((gamma-1)*Vm*F/(R*T))*Na_o^3*Ca_i)/ ...
            ((K_mNa^3+Na_o^3)*(K_mCa+Ca_o)*(1+K_sat*exp((gamma-1)*Vm*F/(R*T))));
        IbCa = Cm*g_B_Ca*(Vm-E_Ca);
        IbNa = Cm*g_B_Na*(Vm-E_Na);
        IpCa = Cm*I_CaP_max*(Ca_i)/(Ca_i+0.0005);

        Irel = K_rel*u^2*v*w*(Ca_rel - Ca_i);
        Itr = (Ca_up-Ca_rel)/tau_tr;
        Iup = I_up_max/(1+(K_up_max/Ca_i));
        Iupleak = Ca_up*I_up_max/Ca_up_max;

        Iion = INa+IK1+Ito+IKur+IKr+IKs+ICaL+IpCa+INaK+INaCa+IbNa+IbCa;

        %alpha and beta computation
        if abs(Vm-47.13)<1e-10
            am=3.2;
        else
            am = 0.32*(Vm+47.13)/(1-exp(-0.1*(Vm+47.13)));
        end

        bm = 0.08*exp(-Vm/11);

        if Vm >= -40
            ah = 0;
            bh = (0.13*(1+exp(-(Vm+10.66)/11.1)))^(-1);
            aj = 0;
            bj = 0.3*(exp(-2.535*10^(-7)*Vm))/(1+exp(-0.1*(Vm+32)));
        else
            ah = 0.135*exp(-(Vm+80)/6.8);
            bh = 3.56*exp(0.079*Vm)+3.1*10^5*exp(0.35*Vm);
            aj = (-127140*exp(0.2444*Vm)-3.474*10^(-5)* ...
                exp(-0.04391*Vm))*(Vm+37.78)/(1+exp(0.311*(Vm+79.23)));
            bj = 0.1212*(exp(-0.01052*Vm))/(1+exp(-0.1378*(Vm+40.14)));
        end

        aoa = 0.65*(exp(-(Vm+10)/8.5)+exp(-(Vm-30)/59))^(-1);
        boa = 0.65*(2.5+exp((Vm+82)/17))^(-1);
        aoi = (18.53 + exp((Vm+113.17)/10.59))^(-1);
        boi = (35.56 + exp(-(Vm+1.26)/7.44))^(-1);

        aua = 0.65*(exp(-(Vm+10)/8.5)+exp(-(Vm-30)/59))^(-1);
        bua = 0.65*(2.5+exp((Vm+82)/17))^(-1);
        aui = (21+exp(-(Vm-185)/28))^(-1);
        bui = exp((Vm-158)/16);

        if abs(Vm-14.1)<1e-10
            axr=0.0015;
        else
            axr = 0.0003*(Vm+14.1)/(1-exp(-(Vm+14.1)/5));
        end
        if abs(Vm-3.3328)<1e-10
            bxr= 3.7863e-4;
        else
            bxr = 7.3898*10^(-5)*(Vm-3.3328)/(-1+exp((Vm-3.3328)/5.1237));

        end

        if abs(Vm-19.9)<1e-10
            axs=6.8e-4;
            bxs=3.15e-4;
        else
            axs = 4*10^(-5)*(Vm-19.9)/(1-exp(-(Vm-19.9)/17));
            bxs = 3.5*10^(-5)*(Vm-19.9)/(exp((Vm-19.9)/9)-1);
        end

        % Tau values
        taum=(am+bm)^(-1);
        tauh=(ah+bh)^(-1);
        tauj=(aj+bj)^(-1);
        tauoa=((aoa+boa)^(-1))/K_Q10;
        tauoi=((aoi+boi)^(-1))/K_Q10;
        tauua=((aua+bua)^(-1))/K_Q10;
        tauui=((aui+bui)^(-1))/K_Q10;
        tauxr=(axr+bxr)^(-1);
        tauxs=(1/2)*(axs+bxs)^(-1);
        if abs(Vm+10)<1e-10
            taud=2.2894;
        else
            taud=(1-exp(-(Vm+10)/6.24))/(0.035*(Vm+10)*(1+exp(-(Vm+10)/6.24)));

        end
        Fn = 10^(-12)*V_rel*Irel-((5*10^(-13))/F)*(ICaL/2 - INaCa/5);

        tauf=9*(0.0197*exp(-0.0337^2*(Vm+10)^2)+0.02)^(-1);
        tauv=1.91+2.09*(1+exp(-(Fn-3.4175*10^(-13))/(13.67*10^(-16))))^(-1);


        if abs(Vm-7.9)<1e-10
            tauw=0.9231;
        else
            tauw=6.0*(1-exp(-(Vm-7.9)/(5)))/((1+0.3*exp(-(Vm-7.9)/(5)))*(Vm-7.9));

        end


        minf=am*taum;
        hinf=ah*tauh;
        jinf=aj*tauj;
        oainf=(1+exp(-(Vm+20.47)/17.54))^(-1);
        oiinf=(1+exp((Vm+43.1)/5.3))^(-1);
        uainf=(1+exp(-(Vm+30.3)/9.6))^(-1);
        uiinf=(1+exp((Vm-99.45)/27.48))^(-1);
        xrinf=(1+exp(-(Vm+14.1)/6.5))^(-1);
        xsinf=(1+exp(-(Vm-19.9)/12.7))^(-1/2);
        dinf=(1+exp(-(Vm+10)/8))^(-1);
        finf=(1+exp((Vm+28)/6.9))^(-1);
        fCainf=(1+Ca_i/0.00035)^(-1);
        uinf=(1+exp(-(Fn-3.4175*10^(-13))/(13.67*10^(-16))))^(-1);
        vinf=1-(1+exp(-(Fn-6.835*10^(-14))/(13.67*10^(-16))))^(-1);
        winf=1-(1+exp(-(Vm-40)/17))^(-1);

        %Concentrations update
        Na_i = Na_i+dt*(-3*INaK-3*INaCa-IbNa-INa)/(F*V_i);
        K_i = K_i+dt*(2*INaK-IK1-Ito-IKur-IKr-IKs)/(F*V_i); %trascuro IbK (background potassio)
        B1 = (2*INaCa-IpCa-ICaL-IbCa)/(2*F*V_i)+(V_up*(Iupleak-Iup)+Irel*V_rel)/V_i;
        B2 = 1 + (TRPN_max*Km_TRPN)/(Ca_i+Km_TRPN)^2 + (CMDN_max*Km_CMDN)/(Ca_i+Km_CMDN)^2;
        Ca_i = Ca_i+dt*B1/B2;
        Ca_up = Ca_up+dt*(Iup-Iupleak-Itr*V_rel/V_up);
        Ca_rel = Ca_rel+dt*((Itr-Irel)*(1+(CSQN_max*Km_CSQN)/(Ca_rel+Km_CSQN)^2)^(-1));





        % Membranepotential update
        Vm = Vm+dt*(-Iion/Cm+Ie+dv2dt);

        % Gating variabls update
        m=minf+(m-minf)*exp(-dt/taum);
        h=hinf+(h-hinf)*exp(-dt/tauh);
        j=jinf+(j-jinf)*exp(-dt/tauj);
        oa=oainf+(oa-oainf)*exp(-dt/tauoa);
        oi=oiinf+(oi-oiinf)*exp(-dt/tauoi);
        ua=uainf+(ua-uainf)*exp(-dt/tauua);
        ui=uiinf+(ui-uiinf)*exp(-dt/tauui);
        xr=xrinf+(xr-xrinf)*exp(-dt/tauxr);
        xs=xsinf+(xs-xsinf)*exp(-dt/tauxs);
        d=dinf+(d-dinf)*exp(-dt/taud);
        f=finf+(f-finf)*exp(-dt/tauf);
        fCa=fCainf+(fCa-fCainf)*exp(-dt/tau_fCa);
        u=uinf+(u-uinf)*exp(-dt/tau_u);
        v=vinf+(v-vinf)*exp(-dt/tauv);
        w=winf+(w-winf)*exp(-dt/tauw);


    end

end