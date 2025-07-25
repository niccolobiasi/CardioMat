function [Vm,state]=TP_2006(Vm,state,Ie,dv2dt,domain,dt,scalar)
% This function implements the TenTusscher-Panfilov (TP) 
% model ( https://doi.org/10.1152/ajpheart.00109.2006)for the CardioMat
%  toolbox. Epicardial (domain=1), endocardial (domain=2), and,
% midmyocardial (dmain=3) parameter sets are included.

if nargin==0
    % %standard initial conditions 
    % Vm = -86.2;
    % Nai = 7.67;
    % Cai = 0.00007;
    % CaSR = 1.3;
    % CaSS =0.00007;
    % Ki = 138.3;
    % m = 0;
    % h = 0.75;
    % j = 0.75;
    % d = 0;
    % f = 1;
    % f2 = 1; 
    % fcass = 1; 
    % s = 1;
    % r = 0;
    % Xs = 0;
    % Xr1 = 0;
    % Xr2 = 1;
    % R_m= 1;

    % CellML initial conditions
    Vm = -85.23;
    Ki = 136.89;
    Nai = 8.604;
    Cai = 0.000126;
    Xr1 = 0.00621;
    Xr2 = 0.4712;
    Xs = 0.0095;
    m = 0.00172;
    h = 0.7444;
    j = 0.7045;
    CaSS = 0.00036;
    d = 3.373e-5;
    f = 0.7888;
    f2 = 0.9755;
    fcass = 0.9953;
    s = 0.999998;
    r = 2.42e-8;
    CaSR = 3.64;
    R_m = 0.9073;
else
    [Vm,Nai,Cai,CaSR,CaSS,Ki,m,h,j,d,f,f2,fcass,s,r,Xs,Xr1,Xr2,R_m]=...
        arrayfun(@TP_step,Vm,state{1},state{2},state{3},state{4},state{5},...
        state{6},state{7},state{8},state{9},state{10},state{11},state{12},...
        state{13},state{14},state{15},state{16},state{17},state{18}, ...
        Ie,dv2dt,domain,scalar);
end

state{1}=Nai;
state{2}=Cai;
state{3}=CaSR;
state{4}=CaSS;
state{5}=Ki;
state{6}=m;
state{7}=h;
state{8}=j;
state{9}=d;
state{10}=f;
state{11}=f2;
state{12}=fcass;
state{13}=s;
state{14}=r;
state{15}=Xs;
state{16}=Xr1;
state{17}=Xr2;
state{18}=R_m;

    function [Vm,Nai,Cai,CaSR,CaSS,Ki,m,h,j,d,f,f2,fcass,s,r,Xs,Xr1,Xr2,R_m]=...
            TP_step(Vm,Nai,Cai,CaSR,CaSS,Ki,m,h,j,d,f,f2,fcass,s,r,Xs,Xr1,Xr2,R_m,Ie,dv2dt,region,scalar)

        R=8314.3; %J*K^(-1)*Mol^(-1)
        T=310;	%K
        F=96486.7; %C/mmol
        Cm=0.185;  %uF               Membrane capacitance
        Nao=140.0; %mM,				 Extracellular sodium concentration
        Cao=2.0; %mM,				 Extracellular calcium concentration
        Ko=5.4; %mM,				 Extracellular potassium concentration
        p_KNa=0.03; %dimensionless
        k_NaCa=1000; %pA/pF          Maximal INaCa
        gamma=0.35; %dimensionless,	 Voltage dependence parameter of INaCa (n)
        alpha=2.5; %dimensionless,	 Factor enhancing outward nature of INaCa
        K_mNai=87.5; %mM,			 Nai half saturation constant for INaCa
        K_mCa=1.38; %mM,			 Cai half saturation constant for INaCa
        K_sat=0.1; %dimensionless,	 Saturation factor for INaCa
        P_NaK=2.724; %pA/pF,	     Maximal INaK (knak)
        K_mK=1.0; %mM,				 Ko half saturation constant for INaK
        K_mNa=40.0; %mM,			 Nai half saturation constant for INaK
        K_pCa=0.0005; %mM            Cai half saturation constant of IpCa
        Vmax_up=0.006375; %mM/ms,	 Maximal Iup
        K_up=0.00025; %mM,			 Half saturation constant of Iup
        V_leak=0.00036; %1/ms,	     Maximal Ileak
        Vrel=0.102; %mM/ms      (wrong in publication)  Maximal Irel
        Vxfer=0.0038; %mM/ms         Maximal Ixfer
        k10=0.15; %mM^(-2)*ms^(-1)   R to O and RI to I Irel transition rate
        k20=0.045; %mM^(-1)*ms^(-1)  O to I and R to RI Irel transition rate
        k3=0.060; %ms^(-1)           O ro R and I to RI Irel transition rate
        k4=0.005; %ms^(-1) (wrong in publication) I to O and RI to I Irel transition rate       
        EC=1.5; %mM                  CaSR half-saturation constant of kcasr
        maxsr=2.5; %dimensionless    Maximum value of kcasr
        minsr=1;%dimensionless       Minimum value of kcasr
        Bufc=0.2; %mM,				 Total cytoplasmic buffer concentration
        Kbufc=0.001; %mM,			 Cai half saturation for cytoplasmic buffer
        Bufsr=10; %mM,				 Total SR Ca buffer concentration
        Kbufsr=0.3; %mM,			 CaSR half saturation for sarcoplasmic buffer
        BufSS=0.4; %mM               Total subspace buffer concentration
        KbufSS=0.00025; %mM          Cass half-saturation constant for subspace buffer
        Vc=0.016404; %nL,			 Cytoplasmic volume
        Vsr=0.001094; %nL,	     Sarcoplasmic reticulum volume
        Vss=0.00005468; %nL       Subspace volume

        %Conductances [nS/pF] 
        GNa=14.838; %[cm^3*uF^(-1)*s^(-1)]
        GK1=5.405;
        if region==2 %endo
            Gto=0.073;
            GKs=0.392;
        elseif region==3 %mid
            Gto=0.294;
            GKs=0.098;
        else  %epi
            Gto=0.294;
            GKs=0.392;
        end
        GKr=0.153;
        GCaL=0.0000398;
        Gpk=0.0146;
        GpCa=0.1238;
        GbNa=0.00029;
        GbCa=0.000592;

        %Nerst Potentials (no changes) 
        EK=((R*T)/F)*log(Ko/Ki);
        EKs=((R*T)/F)*log((Ko+(p_KNa*Nao))/(Ki+(p_KNa*Nai)));
        ENa=((R*T)/F)*log(Nao/Nai);
        ECa=((R*T)/(2*F))*log(Cao/Cai);

        %fast sodium current (no changes) 
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
        ICaL=GCaL*d*f*f2*fcass*4*((Vm-15)*F^2)/(R*T)*((0.25*CaSS*exp((2*(Vm-15)*F)/(R*T))-Cao)/(exp((2*(Vm-15)*F)/(R*T))-1));
        d_inf=1/(1+exp((-8-Vm)/7.5)); 
        ad=(1.4/(1+exp((-35-Vm)/13)))+0.25; 
        bd=1.4/(1+exp((Vm+5)/5)); 
        yd=1/(1+exp((50-Vm)/20)); 
        tau_d=ad*bd+yd; 
        f_inf=1/(1+exp((Vm+20)/7)); 
        af=1102.5*exp(-((Vm+27)/15)^2); 
        bf=200/(1+exp((13-Vm)/10)); 
        yf=180/(1+exp((Vm+30)/10))+20; 
        tau_f=af+bf+yf; 
        f2_inf=0.67/(1+exp((Vm+35)/7))+0.33; 
        af2=600*exp(-((Vm+25)^2/170)); 
        bf2=31/(1+exp((25-Vm)/10)); 
        yf2=16/(1+exp((Vm+30)/10)); 
        tau_f2=af2+bf2+yf2; 
        fcass_inf=0.6/(1+(CaSS/0.05)^2)+0.4; 
        tau_fcass=80/(1+(CaSS/0.05)^2)+2; 
     
        %transient outward current (no changes)
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
        aXs=1400/(1+exp((5-Vm)/6))^0.5;
        bXs=1/(1+exp((Vm-35)/15));
        tau_Xs=aXs*bXs+80;

        %rapid delayed rectifier current (no changes)
        IKr=GKr*((Ko/5.4)^0.5)*Xr1*Xr2*(Vm-EK);
        Xr1_inf=1/(1+exp((-26-Vm)/7));
        aXr1=450/(1+exp((-45-Vm)/10));
        bXr1=6/(1+exp((Vm+30)/11.5));
        tau_Xr1=aXr1*bXr1;
        Xr2_inf=1/(1+exp((Vm+88)/24));
        aXr2=3/(1+exp((-60-Vm)/20));
        bXr2=1.12/(1+exp((Vm-60)/20));
        tau_Xr2=aXr2*bXr2;

        %inward rectifier potassium current (no changes)
        bK1=(3*exp(0.0002*(Vm-EK+100))+exp(0.1*(Vm-EK-10)))/(1+exp(-0.5*(Vm-EK)));
        aK1=0.1/(1+exp(0.06*(Vm-EK-200)));
        XK1_inf=aK1/(aK1+bK1);
        IK1=GK1*((Ko/5.4)^0.5)*XK1_inf*(Vm-EK);

        %Na+/Ca2+ exchanger current (no changes)
        INaCa=k_NaCa*(exp(gamma*Vm*F/(R*T))*(Nai)^3*Cao-(exp((gamma-1)*Vm*F/(R*T))*Nao^3*Cai*alpha))/((K_mNai^3+Nao^3)*(K_mCa+Cao)*(1+K_sat*exp((gamma-1)*Vm*F/(R*T))));

        %Na+/K+ pump current (no changes)
        INaK=P_NaK*Ko*Nai/((Ko+K_mK)*(Nai+K_mNa)*(1+0.1245*exp((-0.1*Vm*F)/(R*T))+0.0353*exp((-Vm*F)/(R*T))));
        IpCa=GpCa*Cai/(K_pCa+Cai);
        IpK=Gpk*(Vm-EK)/(1+exp((25-Vm)/5.98));
        IbNa=GbNa*(Vm-ENa);
        IbCa=GbCa*(Vm-ECa);

        %calcium dynamics
        kcasr=maxsr-((maxsr-minsr)/(1+(EC/CaSR)^2));
        k1=k10/kcasr; 
        k2=k20*kcasr;
        R_m=R_m+dt*(-k2*CaSS*R_m+k4*(1-R_m));
        O=(k1*(CaSS)^2*R_m)/(k3+k1*(CaSS)^2);
        Irel=Vrel*O*(CaSR-CaSS); 
        Ileak=V_leak*(CaSR-Cai); 
        Iup=Vmax_up/(1+((K_up)^2/(Cai)^2)); 
        Ixfer=Vxfer*(CaSS-Cai);

        Cai=Cai+dt*(-((IbCa+IpCa-2*INaCa)/(2*Vc*F)*Cm)+(Vsr/Vc)*(Ileak-Iup)+Ixfer)...
            /(1+Bufc*Kbufc/(Kbufc+Cai)^2);
        CaSR=CaSR+dt*(Iup-Ileak-Irel)/(1+Bufsr*Kbufsr/(CaSR+Kbufsr)^2);
        CaSS=CaSS+dt*(-(ICaL/(2*Vss*F)*Cm)+(Vsr/Vss)*Irel-(Vc/Vss)*Ixfer)/...
            (1+BufSS*KbufSS/(CaSS+KbufSS)^2);
                
        %ALTERNATIVELY quadratic solving can be used for calcium dynamics
        % Caibufc=(Cai*Bufc)/(Cai+Kbufc);
        % Casrbufsr=(CaSR*Bufsr)/(CaSR+Kbufsr);
        % Cassbufss=(CaSS*BufSS)/(CaSS+KbufSS);
        % 
        % CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
        % dCasr=dt*(Iup-Ileak-Irel);
        % bjsr=Bufsr-CaCSQN-dCasr-CaSR+Kbufsr;
        % cjsr=Kbufsr*(CaCSQN+dCasr+CaSR);
        % CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
        % 
        % CaSSBuf=BufSS*CaSS/(CaSS+KbufSS);
        % dCaSS=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL/(2*Vss*F)*Cm));
        % bcss=BufSS-CaSSBuf-dCaSS-CaSS+KbufSS;
        % ccss=KbufSS*(CaSSBuf+dCaSS+CaSS);
        % CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;
        % 
        % CaBuf=Bufc*Cai/(Cai+Kbufc);
        % dCai=dt*((-(IbCa+IpCa-2*INaCa)/(2*Vc*F)*Cm)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
        % bc=Bufc-CaBuf-dCai-Cai+Kbufc;
        % cc=Kbufc*(CaBuf+dCai+Cai);
        % Cai=(sqrt(bc*bc+4*cc)-bc)/2;

        %sodium and potassium dynamics (no changes) 
        Nai=Nai+dt*(-(INa+IbNa+3*INaK+3*INaCa)/(Vc*F)*Cm);
        Ki=Ki+dt*(-(IK1+Ito+IKr+IKs-2*INaK+IpK+Ie+dv2dt)/(Vc*F)*Cm);

        % Vm update (no changes) 
        Iion=INa+IK1+Ito+IKr+IKs+ICaL+INaCa+INaK+IpCa+IpK+IbCa+IbNa;
        Vm=dt*(dv2dt+Ie-Iion)+Vm;

        % Rush-Larsen gate variable updates
        m=m_inf+(m-m_inf)*exp(-dt/tau_m);
        h=h_inf+(h-h_inf)*exp(-dt/tau_h);
        j=j_inf+(j-j_inf)*exp(-dt/tau_j);
        d=d_inf+(d-d_inf)*exp(-dt/tau_d);
        f=f_inf+(f-f_inf)*exp(-dt/tau_f);
        f2=f2_inf+(f2-f2_inf)*exp(-dt/tau_f2); 
        fcass=fcass_inf+(fcass-fcass_inf)*exp(-dt/tau_fcass); 
        r=r_inf+(r-r_inf)*exp(-dt/tau_r);
        s=s_inf+(s-s_inf)*exp(-dt/tau_s);
        Xs=Xs_inf+(Xs-Xs_inf)*exp(-dt/tau_Xs);
        Xr1=Xr1_inf+(Xr1-Xr1_inf)*exp(-dt/tau_Xr1);
        Xr2=Xr2_inf+(Xr2-Xr2_inf)*exp(-dt/tau_Xr2);    			     
    end
end