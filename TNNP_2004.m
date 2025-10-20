function [Vm,state,dv2dt]=TNNP_2004(Vm,state,Ie,dv2dt,domain,dt,scalar)
% This function implements the TenTusscher-Noble-Noble-Panfilov (TNNP) 
% model (https://doi.org/10.1152/ajpheart.00794.2003) for the CardioMat 
% toolbox. Epicardial (domain=1), endocardial (domain=2), and midmyocardial
% (dmain=3) parameter sets are included.

if nargin==0
    Vm = -86.2;
    Nai = 11.6;
    Cai = 0.0002;
    CaSR = 0.2;
    Ki = 138.3;
    m = 0;
    h = 0.75;
    j = 0.75;
    d = 0;
    f= 1;
    fCa = 1;
    s = 1;
    r = 0;
    Xs = 0;
    Xr1 = 0;
    Xr2 = 1;
    g = 1;
else
    [Vm,Nai,Cai,CaSR,Ki,m,h,j,d,f,fCa,s,r,Xs,Xr1,Xr2,g,dv2dt]=...
        arrayfun(@TNNP_step,Vm,state{1},state{2},state{3},state{4},state{5},...
        state{6},state{7},state{8},state{9},state{10},state{11},state{12},...
        state{13},state{14},state{15},state{16},Ie,dv2dt,domain,scalar);
end


state{1}=Nai;
state{2}=Cai;
state{3}=CaSR;
state{4} = Ki;
state{5}=m;
state{6} = h;
state{7}=j;
state{8} = d;
state{9}=f;
state{10}=fCa;
state{11}=s;
state{12}=r;
state{13}=Xs;
state{14}=Xr1;
state{15}=Xr2;
state{16}=g;

    function [Vm,Nai,Cai,CaSR,Ki,m,h,j,d,f,fCa,s,r,Xs,Xr1,Xr2,g,dv2dt]=TNNP_step(Vm,Nai,Cai,CaSR,Ki,m,h,j,d,f,fCa,s,r,Xs,Xr1,Xr2,g,Ie,dv2dt,region,scalar)

        R=8314.3; %J*K^(-1)*Mol^(-1)
        T=310;	%K
        F=96486.7; %C/mmol
        Cm=0.185;  %uF               Membrane capacitance
        Nao=140.0; %mM,				 Extracellular sodium concentration
        Cao=2.0; %mM,				Extracellular calcium concentration
        Ko=5.4; %mM,				 Extracellular potassium concentration
        p_KNa=0.03; %dimensionless
        k_NaCa=1000; %pA/pF           Maximal INaCa
        tau_fCa=2.0; %ms,			 fCa gate time constant
        gamma=0.35; %dimensionless,		 Voltage dependence parameter of INaCa
        alpha=2.5; %dimensionless,		 Factor enhancing outward nature of INaCa
        K_mNai=87.5; %mM,			 Nai half saturation constant for INaCa
        K_mCa=1.38; %mM,			 Cai half saturation constant for INaCa
        K_sat=0.1; %dimensionless,		 Saturation factor for INaCa

        P_NaK=1.362; %pA/pF,	     Maximal INaK
        K_mK=1.0; %mM,				 Ko half saturation constant for INaK
        K_mNa=40.0; %mM,			 Nai half saturation constant for INaK

        K_pCa=0.0005; %mM            Cai half saturation constant of IpCa

        a_rel=0.016464; %mM/ms,			 Max CaSR-dependent Irel
        b_rel=0.25; %mM,			 CaSR half saturation constant of Irel
        c_rel=0.008232; %mM/ms,			Max CaSR-independent Irel

        Vmax_up=0.000425; %mM/ms,		Maximal Iup
        K_up=0.00025; %mM,			Half saturation constant of Iup
        V_leak=0.00008; %1/ms,	    Maximal Ileak
        tau_g=2.0; %ms,				 Ca dynamics g gate time constant
        Bufc=0.15; %mM,				 Total cytoplasmic buffer concentration
        Kbufc=0.001; %mM,			 Cai half saturation for cytoplasmic buffer
        Bufsr=10; %mM,				 Total SR Ca buffer concentration
        Kbufsr=0.3; %mM,			 CaSR half saturation for sarcoplasmic buffer
        Vc=0.016404; %nL,			 Cytoplasmic volume
        Vsr=0.001094; %nL,	     Sarcoplasmic reticulum volume

        %Conductances [nS/pF]
        GNa=14.838; %[cm^3*uF^(-1)*s^(-1)]
        GK1=5.405;
        if region==2 %endo
            Gto=0.073;
            GKs=0.245;
        elseif region==3 %mid
            Gto=0.294;
            GKs=0.062;
        else  %epi
            Gto=0.294;
            GKs=0.245;
        end

        GKr=0.096;
        GCaL=0.000175;
        Gpk=0.0146;
        GpCa=0.825;
        GbNa=0.00029;
        GbCa=0.000592;

        %Nerst Potentials
        EK=((R*T)/F)*log(Ko/Ki);
        EKs=((R*T)/F)*log((Ko+(p_KNa*Nao))/(Ki+(p_KNa*Nai)));
        ENa=((R*T)/F)*log(Nao/Nai);
        ECa=((R*T)/(2*F))*log(Cao/Cai);

        %fast sodium current
        INa=GNa*(m)^3*h*j*(Vm-ENa);
        am=1/(1+exp((-60-Vm)/5));
        bm=0.1/(1+exp((Vm+35)/5))+0.1/(1+exp((Vm-50)/200));
        m_inf=1/(1+exp((-56.86-Vm)/9.03))^2;
        tau_m=am*bm;
        if(Vm>=-40)
            ah=0;
            bh=0.77/(0.13*(1+exp((Vm+10.66)/(-11.1))));
            aj=0;
            bj=0.6*exp(0.057*Vm)/(1+exp(-0.1*(Vm+32)));
        else
            ah=0.057*exp(-(Vm+80)/6.8);
            bh=2.7*exp(0.079*Vm)+(3.1*10^5)*exp(0.3485*Vm);
            aj=((-25428*exp(0.2444*Vm)-6.948*10^(-6)*exp(-0.04391*Vm))*(Vm+37.78))/(1+exp(0.311*(Vm+79.23)));
            bj=(0.02424*exp(-0.01052*Vm))/(1+exp(-0.1378*(Vm+40.14)));
        end
        h_inf=1/(1+exp((71.55+Vm)/7.43))^2;
        tau_h=1/(ah+bh);
        j_inf=1/(1+exp((71.55+Vm)/7.43))^2;
        tau_j=1/(aj+bj);

        %L type Ca2+ current
        ICaL=(((GCaL*d*f*fCa*4*(Vm*F^2))/(R*T))*(Cai*exp(2*Vm*F/(R*T))-0.341*Cao))/(exp(2*Vm*F/(R*T))-1);

        ad=(1.4/(1+exp((-35-Vm)/13)))+0.25;
        bd=1.4/(1+exp((Vm+5)/5));
        yd=1/(1+exp((50-Vm)/20));
        tau_d=ad*bd+yd;
        d_inf=1/(1+exp((-5-Vm)/7.5));
        tau_f=1125*exp(-((Vm+27)^2)/240)+80+165/(1+exp((25-Vm)/10));
        f_inf=1/(1+exp((Vm+20)/7));
        afCa=1/(1+(Cai/0.000325)^8);
        bfCa=0.1/(1+exp((Cai-0.0005)/0.0001));
        yfCa=0.2/(1+exp((Cai-0.00075)/0.0008));
        fCa_inf=(afCa+bfCa+yfCa+0.23)/1.46;

        %transient outward current
        Ito=Gto*r*s*(Vm-EK);

        r_inf=1/(1+exp((20-Vm)/6));
        tau_r=9.5*exp(-((Vm+40)^2)/1800)+0.8;

        if region==2 %endo
            s_inf=1/(1+exp(Vm+28)/5); 
            tau_s=1000*exp(-(Vm+67)^2/1000)+8; 
        else  %epi,mid
            s_inf=1/(1+exp((Vm+20)/5));
            tau_s=85*exp(-((Vm+45)^2)/320)+(5/(1+exp((Vm-20)/5)))+3;
        end

        %slow delayed rectifier current
        IKs=GKs*(Xs^2)*(Vm-EKs);

        Xs_inf=1/(1+exp((-5-Vm)/14));
        aXs=1100/(1+exp((-10-Vm)/6))^0.5;
        bXs=1/(1+exp((Vm-60)/20));
        tau_Xs=aXs*bXs;


        %rapid delayed rectifier current
        IKr=GKr*((Ko/5.4)^0.5)*Xr1*Xr2*(Vm-EK);

        Xr1_inf=1/(1+exp((-26-Vm)/7));
        aXr1=450/(1+exp((-45-Vm)/10));
        bXr1=6/(1+exp((Vm+30)/11.5));
        tau_Xr1=aXr1*bXr1;
        Xr2_inf=1/(1+exp((Vm+88)/24));
        aXr2=3/(1+exp((-60-Vm)/20));
        bXr2=1.12/(1+exp((Vm-60)/20));
        tau_Xr2=aXr2*bXr2;


        %inward rectifier potassium current
        bK1=(3*exp(0.0002*(Vm-EK+100))+exp(0.1*(Vm-EK-10)))/(1+exp(-0.5*(Vm-EK)));
        aK1=0.1/(1+exp(0.06*(Vm-EK-200)));
        XK1_inf=aK1/(aK1+bK1);
        IK1=GK1*((Ko/5.4)^0.5)*XK1_inf*(Vm-EK);


        %Na+/Ca2+ exchanger current
        INaCa=k_NaCa*(exp(gamma*Vm*F/(R*T))*(Nai)^3*Cao-(exp((gamma-1)*Vm*F/(R*T))*Nao^3*Cai*alpha))/((K_mNai^3+Nao^3)*(K_mCa+Cao)*(1+K_sat*exp((gamma-1)*Vm*F/(R*T))));


        %Na+/K+ pump current
        INaK=P_NaK*Ko*Nai/((Ko+K_mK)*(Nai+K_mNa)*(1+0.1245*exp((-0.1*Vm*F)/(R*T))+0.0353*exp((-Vm*F)/(R*T))));
        IpCa=GpCa*Cai/(K_pCa+Cai);
        IpK=Gpk*(Vm-EK)/(1+exp((25-Vm)/5.98));
        IbNa=GbNa*(Vm-ENa);
        IbCa=GbCa*(Vm-ECa);

        Ileak=V_leak*(CaSR-Cai);
        Iup=Vmax_up/(1+((K_up)^2/(Cai)^2));
        Irel=(((a_rel*(CaSR)^2)/((b_rel)^2+(CaSR)^2))+c_rel)*d*g;
        if(Cai<0.00035)
            g_inf=1/(1+((Cai)^6/(0.00035)^6));
        else
            g_inf=1/(1+((Cai)^16/(0.00035)^16));
        end

        Cai=Cai+dt*(1/(1+Bufc*Kbufc/(Cai+Kbufc)^2)*(-(ICaL+IbCa+IpCa-2*INaCa)/(2*Vc*F)*Cm+Ileak-Iup+Irel));
        CaSR=CaSR+dt*(1/(1+Bufsr*Kbufsr/(CaSR+Kbufsr)^2)*(Vc/Vsr)*(-Ileak+Iup-Irel));


        %sodium and potassium dynamics
        Nai=Nai+dt*(-(INa+IbNa+3*INaK+3*INaCa)/(Vc*F)*Cm);
        Ki=Ki+dt*(-(IK1+Ito+IKr+IKs-2*INaK+IpK+Ie+dv2dt)/(Vc*F)*Cm);


        % Vm update
        Iion=INa+IK1+Ito+IKr+IKs+ICaL+INaCa+INaK+IpCa+IpK+IbCa+IbNa;
        dv2dt=dv2dt+Ie-Iion;
        Vm=dt*dv2dt+Vm;

        % Rush-Larsen gate variable updates
        m=m_inf+(m-m_inf)*exp(-dt/tau_m);
        h=h_inf+(h-h_inf)*exp(-dt/tau_h);
        j=j_inf+(j-j_inf)*exp(-dt/tau_j);
        d=d_inf+(d-d_inf)*exp(-dt/tau_d);
        f=f_inf+(f-f_inf)*exp(-dt/tau_f);
        if ~(fCa_inf>fCa && Vm>-60)
            fCa=fCa_inf+(fCa-fCa_inf)*exp(-dt/tau_fCa);
        end
        r=r_inf+(r-r_inf)*exp(-dt/tau_r);
        s=s_inf+(s-s_inf)*exp(-dt/tau_s);
        Xs=Xs_inf+(Xs-Xs_inf)*exp(-dt/tau_Xs);
        Xr1=Xr1_inf+(Xr1-Xr1_inf)*exp(-dt/tau_Xr1);
        Xr2=Xr2_inf+(Xr2-Xr2_inf)*exp(-dt/tau_Xr2);
        if ~(g_inf>g && Vm>-60)
            g=g_inf+(g-g_inf)*exp(-dt/tau_g);
        end
    			     
    end
end