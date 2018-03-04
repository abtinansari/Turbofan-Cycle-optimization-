function F_mo=fnc4(X)
%X=[pi_cH,alpha,pi_f]
%%%%% constrained Analysis%%%%%
%This function gets the 
%HP compressor pressure ratio, By-pass ratio & Fan pressure ratio
% as an input and computes the Specific Thrust for a 
% UnMixed flow dual spool Turbofan Engine 
% working fluid is modeled as Modifie Specific Heat(CSH)
%%
%Inputs:
P0_P9=1;
P0_P19=1;
pi_dmax=0.99;
pi_b= 0.9557;
pi_itb=1;
pi_n= 0.99;
pi_fn=0.99;
e_cH=0.9;
e_cL=0.91;
e_tH=0.92;
e_tL=0.94;
e_f=0.91;
eta_b=0.99;
eta_itb=1;
effm_tH=0.99;
effm_tL=1;
Tt4=2429.67;
Tt4p5=1869.67;
hPRb =18600;
hPRitb =18600;
h=0;
T_o=Tamb(h);% computes Tempreture as a function of height
gc=32.174;
conversion=778.16;
hPRb=hPRb*conversion;
hPRitb=hPRitb*conversion;
M0=0.01;
pi_cL=2.1837;
%%
%Check Inputs:
if pi_dmax>1 || pi_b>1 || pi_itb>1 || pi_n>1 || pi_fn>1 || e_cH>1 || e_cL>1 || e_tH>1 || e_tL>1 || e_f>1 || eta_b>1 || eta_itb>1 || effm_tH>1 || effm_tL>1
       F_mo=Inf;
     return
else
                    
%%
fbi=0;
fitbi=0;
[hc,Cpc,gamac]=enthalpy(T_o,0);% computes enthalpy of the mixture as a function of Tempreture & fuel-air ratio
[ht,Cpt,gamat]=enthalpy(Tt4,0);
rc = (gamac - 1) * Cpc / gamac;   %'(1)
rt = (gamat - 1) * Cpt / gamat;  % '(2)
   
if gamac < 0 || gamac==1 
     F_mo=Inf;
     return
else
  a0 = (gamac * rc * T_o * gc) ^ 0.5;   % '(4)
  v0 = a0 * M0;                          %'(5)
  tau_r = (gamac - 1) * (0.5 * M0 ^ 2) + 1; %   '(6)
  pi_r = tau_r ^ (gamac / (gamac - 1)); %         '(7)
if M0 <= 1                               %'(8)
    eta_r = 1;
  end
  if M0 > 1 && M0 <= 5 
    eta_r = 1 - 0.075 * (M0 - 1) ^ 1.35;
  end
  if M0 > 5 
    eta_r = 800/ (935 + M0 ^ 4);
 end
  pi_d = pi_dmax * eta_r;                      % '(9)
  if pi_cL < X(3)
     pi_cL = X(3);
 end
pi_c=pi_cL.*X(1);
 tau_alp_b = Cpt * Tt4 / (Cpc * T_o); %                   '(10)
  tau_cH = X(1).^ ((gamac - 1) / (gamac * e_cH)); %        '(12)
  tau_cL = pi_cL ^ ((gamac - 1) / (gamac * e_cL)) ; %       '(13)
  tau_f = X(3) ^ ((gamac - 1) / (gamac * e_f)); %           '(14)
  tau_c = tau_cH.* tau_cL; %                                 '(15)
  
 if tau_c==0 || tau_f==0 || tau_cH==1 || tau_cL==1 || tau_f==0
       F_mo=Inf;
     return
 else
     
    Tt3 = tau_r * tau_c*T_o;%                        'eqn(20b)
    
% Constraint No.1: HPC Exit Temperature
    if Tt3>(840*9/5)
          F_mo=Inf;
     return
 else
 
   tau_b = Tt4./ (T_o * tau_r * tau_c);                 % '(16)
 
   eta_cH = ((tau_cH.^ e_cH) - 1)./ (tau_cH - 1);    % '(17)
  eta_cL = ((tau_cL.^ e_cL) - 1)./ (tau_cL - 1);      %'(18)
  eta_f = ((tau_f ^ e_f) - 1) / (tau_f - 1);         %'(19)
  
   h_0 = Cpc * T_o;
  temp_fb = eta_b.* hPRb./ h_0 - tau_alp_b;
  
 if temp_fb ==0
        F_mo=Inf;
     return
else
 
  [h_t3,Cp3,gama3]=enthalpy(Tt3,0);
  f_4i = (tau_alp_b - tau_r * tau_c)./ temp_fb; %       'eqn(20c),initial value for f_4i
  %3
  [h_t4,Cpt,gamat]=enthalpy(Tt4,f_4i);
  temp_fb = eta_b * hPRb - h_t4;
  
if temp_fb ==0
        F_mo=Inf;
     return
else
    
     fb = (h_t4 - h_t3)./ temp_fb     ; %                   'eqn(20d)
%%% First Iteration 
while abs(fb - f_4i) > 0.0001
    f_4i = fb;
    [h_t4,Cpt,gamat]=enthalpy(Tt4,f_4i);
    temp_fb = eta_b * hPRb - h_t4;
    
    if temp_fb ==0
        F_mo=Inf;
     break
    end
 
    fb = (h_t4 - h_t3)./ temp_fb     ; %                   'eqn(20d)
    
end
  %%
  %Iteration for cp=cp(T,f) and gama=gama(T,f)
  Err_fb = abs(fbi - fb);
  while Err_fb > 0.0001
      fbi = fb;
      [ht,Cpt,gamat]=enthalpy(Tt4,fb);
      %24
      tau_alp_b = Cpt.* Tt4./ (Cpc * T_o); %                   '(10)
      tau_cH = X(1).^ ((gamac - 1) / (gamac * e_cH)); %        '(12
      tau_cL = pi_cL.^ ((gamac - 1) / (gamac * e_cL)) ;
      tau_f = X(3) ^ ((gamac - 1) / (gamac * e_f)); %           '(14)
      tau_c = tau_cH.* tau_cL;
      
if tau_c==0 || tau_f==0 || tau_cH==1 || tau_cL==1 || tau_f==0
       F_mo=Inf;
       break
end
      tau_b = Tt4./ (T_o * tau_r * tau_c);                 % '(16)
      eta_cH = ((tau_cH.^ e_cH) - 1) / (tau_cH - 1);    % '(17)
      eta_cL = ((tau_cL.^ e_cL) - 1) / (tau_cL - 1);      %'(18)
      eta_f = ((tau_f ^ e_f) - 1) / (tau_f - 1);         %'(19)
      h_0 = Cpc * T_o;
      temp_fb = eta_b * hPRb./ h_0 - tau_alp_b;
    if temp_fb ==0
        F_mo=Inf;
     break
    end   
  Tt3 = T_o * tau_r.* tau_c ; %                        'eqn(20b)
  
% Constraint No.1: HPC Exit Temperature  
   if Tt3>(840*9/5)
          F_mo=Inf;
   break
   end
    
  [h_t3,Cp3,gama3]=enthalpy(Tt3,0);
  f_4i = (tau_alp_b - tau_r * tau_c)./ temp_fb; %       'eqn(20c),initial value for f_4i
  %3
  [h_t4,Cpt,gamat]=enthalpy(Tt4,f_4i);
  temp_fb = eta_b * hPRb - h_t4;
  
if temp_fb ==0
        F_mo=Inf;
     break
end
    
     fb = (h_t4 - h_t3)./ temp_fb     ; %                   'eqn(20d)
     
     while abs(fb - f_4i) > 0.0001
         f_4i = fb;
         [h_t4,Cpt,gamat]=enthalpy(Tt4,f_4i);
         temp_fb = eta_b * hPRb - h_t4;
         
if temp_fb ==0
        F_mo=Inf;
     break
    end
        fb = (h_t4 - h_t3)./ temp_fb     ; %                   'eqn(20d)
     end
     Err_fb = abs(fbi - fb);
  end    
  %%
  tau_tH = 1 - (tau_cH - 1).* tau_cL * tau_r./ ((1 + fb).* tau_alp_b * effm_tH)  ;%'(21)
if tau_tH < 0
     F_mo=Inf;
     return
else  
Tt4p5 = Tt4.* tau_tH;
 
% Constraint No.2: HPT Exit Temperature
if Tt4p5>(1250*9/5)
     F_mo=Inf;
     return
else
pi_tH = tau_tH.^ (gamat./ ((gamat - 1) * e_tH)); %          '(27)
 
% Constraint No.3: HPT Pressure Ratio
if pi_tH<0.0625
     F_mo=Inf;
     return
else
[hitb,Cpitb,gamaitb]=enthalpy(Tt4p5,fb);
tau_alp_itb = Cpitb.* Tt4p5./ (Cpc * T_o);%         '(11)
fitb = 0;
tau_itb = 1;
      %5
  ritb = (gamaitb - 1).* Cpitb./ gamaitb;%                   '(3)
  f = fb + fitb;%                                             '(23)
  tau_tL = 1 - (tau_cL - 1 + X(2).* (tau_f - 1)).* tau_r./ ((1 + fb + fitb).* tau_alp_itb * effm_tL); %'(26)
  Tt5 = Tt4.* tau_tH.*tau_tL;
  
  % Constraint No.4: LPT Exit Temperature
  if Tt5>(1250*9/5)
     F_mo=Inf;
     return
  else
    
  
  pi_tL = tau_tL.^ (gamaitb./ ((gamaitb - 1) * e_tL)) ; %     '(28)
  
  % Constraint No.5: LPT Pressure Ratio 
if pi_tL<0(1/256)
     F_mo=Inf;
     return
else
    
  eta_tH = (1 - tau_tH)./ (1 - tau_tH.^ (1/ e_tH))  ; %    '(29)
  eta_tL = (1 - tau_tL)./ (1 - tau_tL.^ (1 / e_tL)); %      '(30)
  Pt9_P9 = P0_P9 * pi_r * pi_d * pi_cL.* X(1).* pi_b.* pi_tH.* pi_itb.* pi_tL.* pi_n; %'(31)
  
  if Pt9_P9<1 
      F_mo=Inf;
      return
  else
  
  tempM9 = Pt9_P9.^ ((gamaitb - 1)./ gamaitb);
    Tt9_T0 = tau_r * tau_c.* tau_b.* tau_tH.* tau_itb.* tau_tL; %  '(32)
  T9_T0 = Tt9_T0./ Pt9_P9.^ ((gamaitb - 1)./ gamaitb) ; %       '(33)
 M9 = ((tempM9 - 1) * 2./ (gamaitb - 1)).^ 0.5  ; %          '(34
 M8 = M9;
 if M9 >= 1
     M8 = 1;
 end
  v9_a0 = M9.* (gamaitb.* ritb.* T9_T0./ (gamac * rc)).^ 0.5 ; % '(35)
  v9 = v9_a0.* a0;
  v9_v0 = (v9./ v0);
  Pt19_P19 = P0_P19.* pi_r * pi_d * X(3) * pi_fn    ; %          '(36)
  
   if Pt19_P19<1 
      F_mo=Inf;
      return
  else
  tempM19 = Pt19_P19.^ ((gamac - 1) / gamac);
 
M19 = ((tempM19 - 1) * 2./ (gamac - 1)).^ 0.5     ; %       '(37)
  M18 = M19;
  
if M19 >= 1
    M18 = 1;
end
Tt19_T0 = tau_r * tau_f              ; %                        '(38)
  T19_T0 = Tt19_T0./ Pt19_P19.^ ((gamac - 1) / gamac)  ; %       '(39)
  v19_a0 = M19.* T19_T0.^ 0.5    ; %                              '(40)
  v19 = v19_a0.* a0;
  v19_v0 = (v19./ v0);
  fc_mc = (1 + f).* v9_a0 - M0 + (1 + f).* (ritb./ rc).* (T9_T0./ v9_a0).* (1 - P0_P9).* (1 / gamac); %'(41)
  ff_mf = v19_a0 - M0 + (T19_T0./ v19_a0).* (1 - P0_P19)./ gamac; %'(42)
  
  %%
  %Outputs:
  F_mo = (a0 / gc) * (fc_mc + X(2).* ff_mf)./ (1 + X(2)); %   '(43)
  
  if F_mo < 0
      F_mo=Inf;
      return
  else
      
  % Constraint No.6: Specific Thrust   
  if F_mo < 30000/(4.44823*1495)
      F_mo=Inf;
      return
  else
      
end
end
 end
 
end
     end
  end
 end
 end
     end
      end
  end
 end
     end
end
end
  
