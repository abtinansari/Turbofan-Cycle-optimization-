function Tamb=Tamb(h)
% computes Standard Tempreture as a function of height
%Calculation is for English unit
T0ref=288.15; %K
r0 = 6356.66; % 'km
h_eng = h / 3.2808; % 'kft to km
z = r0 * h_eng / (r0 + h_eng); % 'geopotential height
if z >= 0 && z <= 11
    Temp = T0ref - 6.5 * (z); % 'Kelvin
    Tamb = Temp * 1.8 ; %'convert K to R
elseif z > 11 && z <= 20
    Temp = 216.65;
    Tamb = Temp * 1.8;
    elseif z > 20 && z <= 32
    Temp = 216.65 + 1 * (z - 20);
    Tamb = Temp * 1.8;
    elseif z > 32 &&z <= 47
    Temp = 228.65 + 2.8 * (z - 32);
    Tamb = Temp * 1.8;
    elseif z > 47 && z <= 51
    Temp = 270.65;
    Tamb = Temp * 1.8;
    elseif z > 51 z <= 71
    Temp = 270.65 - 2.8 * (z - 51);
    Tamb = Temp * 1.8;
    elseif z > 71 && z <= 84.852 
    Temp = 214.65 - 2 * (z - 71);
    Tamb = Temp * 1.8;
end
end
