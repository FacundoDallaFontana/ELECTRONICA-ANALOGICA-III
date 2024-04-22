clear;
clc;
%esto es solo para cambiar la current folder
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end
%% Parametros
f=2e9;
Zo=50;
Zin=50;
Zout=50;
VCE=1;%V
IC=35;%mA
e_r = 4.3;%depende del dielectrico
H=1.2e-3%m
t=50e-6%m
%%
file = ['BFP640_w_noise_VCE_',num2str(VCE,1),'.0V_IC_',num2str(IC,2),'mA.s2p']
list = sparameters(file);
i = find(list.Frequencies==f,1);
S = list.Parameters(:,:,i);

%%  Prueba de Estabilidad
Delta = S(1,1)*S(2,2)-S(1,2)*S(2,1);
abs_Delta = abs(Delta) %Debe ser <1
k = (1 - abs(S(1,1))^2 - abs(S(2,2))^2 + abs_Delta^2)/(2*abs(S(1,2)*S(2,1))) %Debe ser >1 para ser incondicionalmente estable

%% Calculos
B1 = 1 + abs(S(1,1))^2 - abs(S(2,2))^2 - abs_Delta^2;
B2 = 1 + abs(S(2,2))^2 - abs(S(1,1))^2 - abs_Delta^2;
C1 = S(1,1) - (Delta*conj(S(2,2)));
C2 = S(2,2) - (Delta*conj(S(1,1)));
if B1>0
    r_Ms = (B1 - sqrt(B1^2-4*abs(C1)^2))/(2*C1);
else
    r_Ms = (B1 + sqrt(B1^2-4*abs(C1)^2))/(2*C1);
end

if B2>0
    r_ML = (B2 - sqrt(B2^2-4*abs(C2)^2))/(2*C2);
else
    r_ML = (B2 + sqrt(B2^2-4*abs(C2)^2))/(2*C2);
end
%% Coef. de Ref.
r_s = r_Ms;
r_in = conj(r_Ms);% r_in = rs*
r_L = r_ML;
r_out = conj(r_ML);% r_out = rL*

%% GTMax
G_Tmax=((1-abs(r_Ms)^2)/abs(1-S(1,1)*r_Ms)^2)*(abs(S(2,1))^2)*(1/(1-abs(r_ML)))

%% Impedancias
Z_in = Zo*(1+r_in)/(1-r_in);
Z_out = Zo*(1+r_out)/(1-r_out);
Z_s = conj(Z_in);% Zs = Zin* = Zo x (1+rs)/(1-rs)
Z_L = conj(Z_out);% ZL = Zout* = Zo x (1+rL)/(1-rL)

%% Calculo de microtiras (Hammerstad)

A = (Zo/60)*sqrt((e_r+1)/2)+((e_r-1)/(e_r+1))*(0.23+(0.11/e_r));
B = (377*pi)/(2*Zo*sqrt(e_r));

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
    Z0 = ((120*pi)/sqrt(e_rp))/(W_H+1.393+0.667*log(W_H+1.444))
else
    e_rp = ((e_r+1)/2) + ((e_r-1)/2)*((1/(sqrt(1+(12*H)/W)))+0.004*(1-W_H)^2)
    Z0 = (60/sqrt(e_rp))*log(((8*H)/W)+(W/(4*H)))
end

%% Otros calculos
Lambda_0 = (3e8)/f
Lambda_p = Lambda_0/sqrt(e_rp)

%% Calculo Acoplador de entrada
R_in = real(Z_in)
X_in = imag(Z_in)
R_inp = R_in * (1+(X_in/R_in)^2)
X_inp = R_in * (R_inp/X_in)
Z=sqrt(R_inp*Zo)
[w_in, imp_in, Lambda_p_in] = w_microstrip(e_r, H, t, Z, f)
Largo_ac_in = Lambda_p_in/4

%% Calculo Acoplador de Salida
R_out = real(Z_out)
X_out = imag(Z_out)
%%aca no pasamos a paralelo(revisar notas del profe)
Z=sqrt(R_out*Zo)
[w_out, imp_out, Lambda_p_out] = w_microstrip(e_r, H, t, Z, f)
Largo_ac_out = Lambda_p_out/4
