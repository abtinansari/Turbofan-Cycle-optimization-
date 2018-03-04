function [h,cpf,gamaf]=enthalpy(T,ff)
% computes enthalpy of the mixture as a function of 
%Tempreture & fuel-air ratio
conversion=778.16;
Aa=zeros(1,8);
Ap=zeros(1,8);
Rgas = 53.3366; %'ft-lbf/(lbm*R)
%air
Aa(1) = 0.25020021;
Aa(2) = -0.000051536879;
Aa(3) = 0.000000065519486;
Aa(4) = -6.7178376E-12;
Aa(5) = -1.5128259E-14;
Aa(6) = 7.6215767E-18;
Aa(7) = -1.452677E-21;
Aa(8) = 1.011554E-25;
%'product
Ap(1) = 0.073816638;
Ap(2) = 0.001225863;
Ap(3) = -0.0000013771901;
Ap(4) = 9.9686793E-10;
Ap(5) = -4.2051104E-13;
Ap(6) = 1.0212913E-16;
Ap(7) = -1.3335668E-20;
Ap(8) = 7.267871E-25;
  href_a = -1.7558886; % 'Btu/lbm
  href_p = 30.58153; % 'Btu/lbm
  %'h_air
  h_a = (href_a + Aa(1).* T + Aa(2).*T.* T / 2 + Aa(3) * T.^ 3 / 3 + Aa(4) * T.^ 4 / 4 + Aa(5) * T.^ 5 / 5 + Aa(6) * T.^ 6 / 6 + Aa(7) * T.^ 7 / 7 + Aa(8) * T.^ 8 / 8);
  %h_product
  h_p = (href_p + Ap(1) * T + Ap(2) * T / 2.* T + Ap(3) * T.^ 3 / 3 + Ap(4) * T.^ 4 / 4 + Ap(5) * T.^ 5 / 5 + Ap(6) * T.^ 6 / 6 + Ap(7) * T.^ 7 / 7 + Ap(8) * T.^ 8 / 8);
  %mixture
  h = (h_a + ff.* h_p)./ (1 + ff) ; %'Btu/lbm
  h = h * conversion;
  cp_a = Aa(1) + Aa(2) * T + Aa(3) * T.* T + Aa(4) * T.^ 3 + Aa(5) * T.^ 4 + Aa(6) * T.^ 5 + Aa(7) * T.^ 6 + Aa(8) * T.^ 7;
  cp_p = Ap(1) + Ap(2) * T + Ap(3) * T.* T + Ap(4) * T.^ 3 + Ap(5) * T.^ 4 + Ap(6) * T.^ 5 + Ap(7) * T.^ 6 + Ap(8) * T.^ 7;
  cpf = (cp_a + ff.* cp_p)./ (1 + ff); % 'Btu/lbmR
  cpf = cpf * conversion;
  gamaf = cpf./ (cpf - Rgas);
end
