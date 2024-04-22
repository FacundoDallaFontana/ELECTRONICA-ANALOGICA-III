function [W, Zo]= w_microstrip(e_r, H, t, Z)
%% Calculo de microtiras (Hammerstad)

A = (Z/60)*sqrt((e_r+1)/2)+((e_r-1)/(e_r+1))*(0.23+(0.11/e_r));
B = (377*pi)/(2*Z*sqrt(e_r));

W_H = (8*exp(A))/(exp(2*A)-2);

if W_H > 2
    W_H = (2/pi)*(B-1-log(2*B-1))+((e_r-1)/(2*e_r))*(log(B-1)+0.39-(0.61/e_r))
end

W = W_H*H

if W_H <= (1/(2*pi))
    We = W + (t/pi)*(1+log((4*pi*W)/t))
else
    We = W + (t/pi)*(1+log((2*H)/t))
end

if W_H >= 1
    e_rp = ((e_r+1)/2) + ((e_r-1)/2)*(1/(sqrt(1+(12*H)/W)))
    Zo = ((120*pi)/sqrt(e_rp))/(W_H+1.393+0.667*log(W_H+1.444))
else
    e_rp = ((e_r+1)/2) + ((e_r-1)/2)*((1/(sqrt(1+(12*H)/W)))+0.004*(1-W_H)^2)
    Zo = (60/sqrt(e_rp))*log(((8*H)/W)+(W/(4*H)))
end
