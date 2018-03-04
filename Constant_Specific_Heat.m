function SFC=fnc(X)
%X=[pi_cH]
%This function gets the HP compressor pressure ratio
% as an input and computes the SFC for a 
% UnMixed flow dual spool Turbofan Engine 
% working fluid is modeled as Constant Specific Heat(CSH)
%Inputs:
P0_P9=1;
P0_P19=1;
pi_dmax= 1;
pi_b= 0.9557;
pi_itb=1;
pi_n= 1;
pi_fn=0.99;
e_cH=0.947;
e_cL=0.87;
e_tH=0.94;
e_tL=1;
e_f=0.96;
eta_b=1;
eta_itb=1;
effm_tH=0.85;
effm_tL=0.87;
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
M0=0;
alpha=5.05;
pi_c=21.4966;
pi_f=1.5374;
%%
fbi=0;
fitbi=0;
[hc,Cpc,gamac]=enthalpy(T_o,0);% computes enthalpy of the mixture as a function of Tempreture & fuel-air ratio
[ht,Cpt,gamat]=enthalpy(Tt4,0);
rc = (gamac - 1) * Cpc / gamac;   %'(1)
rt = (gamat - 1) * Cpt / gamat;  % '(2)
 
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
pi_cL = pi_c / X;
 if pi_cL < pi_f
     pi_cL = pi_f;
 end
 
 tau_alp_b = Cpt * Tt4 / (Cpc * T_o); %                   '(10)
  tau_cH = X ^ ((gamac - 1) / (gamac * e_cH)); %        '(12)
  tau_cL = pi_cL ^ ((gamac - 1) / (gamac * e_cL)) ; %       '(13)
  tau_f = pi_f ^ ((gamac - 1) / (gamac * e_f)); %           '(14)
  tau_c = tau_cH * tau_cL; %                                 '(15)
   tau_b = Tt4 / (T_o * tau_r * tau_c);                 % '(16)
 
   eta_cH = ((tau_cH ^ e_cH) - 1) / (tau_cH - 1);    % '(17)
  eta_cL = ((tau_cL ^ e_cL) - 1) / (tau_cL - 1);      %'(18)
  eta_f = ((tau_f ^ e_f) - 1) / (tau_f - 1);         %'(19)
  
   h_0 = Cpc * T_o;
  temp_fb = eta_b * hPRb / h_0 - tau_alp_b;
  
  fb = (tau_alp_b - tau_r * tau_c) / temp_fb;%       '(20a)
  
  %%
  %Iteratio for cp=cp(T,f) and gama=gama(T,f)
  Err_fb = abs(fbi - fb);
  i=0;
  while Err_fb > 0.0001
      i=i+1;
      fbi = fb;
      [ht,Cpt,gamat]=enthalpy(Tt4,fb);
      %24
      tau_alp_b = Cpt * Tt4 / (Cpc * T_o); %                   '(10)
      tau_cH = X ^ ((gamac - 1) / (gamac * e_cH)); %        '(12
      tau_cL = pi_cL ^ ((gamac - 1) / (gamac * e_cL)) ;
      tau_f = pi_f ^ ((gamac - 1) / (gamac * e_f)); %           '(14)
      tau_c = tau_cH * tau_cL;
      tau_b = Tt4 / (T_o * tau_r * tau_c);                 % '(16)
      eta_cH = ((tau_cH ^ e_cH) - 1) / (tau_cH - 1);    % '(17)
      eta_cL = ((tau_cL ^ e_cL) - 1) / (tau_cL - 1);      %'(18)
      eta_f = ((tau_f ^ e_f) - 1) / (tau_f - 1);         %'(19)
      h_0 = Cpc * T_o;
      temp_fb = eta_b * hPRb / h_0 - tau_alp_b;
      fb = (tau_alp_b - tau_r * tau_c) / temp_fb;%       '(20a)
      Err_fb = abs(fbi - fb);
  end
%%
tau_tH = 1 - (tau_cH - 1) * tau_cL * tau_r / ((1 + fb) * tau_alp_b * effm_tH)  ;%'(21)
Tt4p5 = Tt4 * tau_tH;
[hitb,Cpitb,gamaitb]=enthalpy(Tt4p5,fb);
tau_alp_itb = Cpitb * Tt4p5 / (Cpc * T_o);%         '(11)
fitb = 0;
tau_itb = 1;
  ritb = (gamaitb - 1) * Cpitb / gamaitb;%                   '(3)
  f = fb + fitb;%                                             '(23)
  tau_tL = 1 - (tau_cL - 1 + alpha * (tau_f - 1)) * tau_r / ((1 + fb + fitb) * tau_alp_itb * effm_tL); %'(26)
  pi_tH = tau_tH ^ (gamat / ((gamat - 1) * e_tH)); %          '(27)
  pi_tL = tau_tL ^ (gamaitb / ((gamaitb - 1) * e_tL)) ; %     '(28)
  eta_tH = (1 - tau_tH) / (1 - tau_tH ^ (1/ e_tH))  ; %    '(29)
  eta_tL = (1 - tau_tL) / (1 - tau_tL ^ (1 / e_tL)); %      '(30)
  Pt9_P9 = P0_P9 * pi_r * pi_d * pi_cL * X * pi_b * pi_tH * pi_itb * pi_tL * pi_n; %'(31)
  tempM9 = Pt9_P9 ^ ((gamaitb - 1) / gamaitb);
    Tt9_T0 = tau_r * tau_c * tau_b * tau_tH * tau_itb * tau_tL; %  '(32)
  T9_T0 = Tt9_T0 / Pt9_P9 ^ ((gamaitb - 1) / gamaitb) ; %       '(33)
 M9 = ((tempM9 - 1) * 2 / (gamaitb - 1)) ^ 0.5  ; %          '(34
 M8 = M9;
 if M9 >= 1
     M8 = 1;
 end
  v9_a0 = M9 * (gamaitb * ritb * T9_T0 / (gamac * rc)) ^ 0.5 ; % '(35)
  v9 = v9_a0 * a0;
  v9_v0 = v9 / v0;
  Pt19_P19 = P0_P19 * pi_r * pi_d * pi_f * pi_fn    ; %          '(36)
  tempM19 = Pt19_P19 ^ ((gamac - 1) / gamac);
M19 = ((tempM19 - 1) * 2 / (gamac - 1)) ^ 0.5     ; %       '(37)
  M18 = M19;
  
if M19 >= 1
    M18 = 1;
end
Tt19_T0 = tau_r * tau_f              ; %                        '(38)
  T19_T0 = Tt19_T0 / Pt19_P19 ^ ((gamac - 1) / gamac)  ; %       '(39)
  v19_a0 = M19 * T19_T0 ^ 0.5    ; %                              '(40)
  v19 = v19_a0 * a0;
  v19_v0 = v19 / v0;
  fc_mc = (1 + f) * v9_a0 - M0 + (1 + f) * (ritb / rc) * (T9_T0 / v9_a0) * (1 - P0_P9) * (1 / gamac); %'(41)
  ff_mf = v19_a0 - M0 + (T19_T0 / v19_a0) * (1 - P0_P19) / gamac; %'(42)
  
  %%
  %Outputs:
  F_mo = (a0 / gc) * (fc_mc + alpha * ff_mf) / (1 + alpha); %   '(43)
  unitchange = 3600;
  SFC= f / ((1+ alpha) * F_mo) * unitchange; %  '(44)   
 
end
