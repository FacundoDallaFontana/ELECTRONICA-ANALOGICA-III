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
H=1.66e-3%m
t=400e-6%m
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

[W_50, Z_50, Lambda_p_50] = w_microstrip(e_r, H, t, Zo, f)

%% Calculo Acoplador de entrada
R_in = real(Z_in)
X_in = imag(Z_in)
R_inp = R_in * (1+(X_in/R_in)^2)
X_inp = R_in * (R_inp/X_in)
Z1=sqrt(R_inp*Zo)
[w_in, imp_in, Lambda_p_in] = w_microstrip(e_r, H, t, Z1, f)
Largo_ac_in = Lambda_p_in/4%Usamos lambda eficaz

cap_in = 1/(2*pi*f*X_inp)
beta = (2*pi)/Lambda_p_in%usamos lambda eficaz
d_cap_in = (1/beta)*acot(X_inp/Zo)%Con impendacia caracteristica 50 ohm

%% Calculo Acoplador de Salida
R_out = real(Z_out)
X_out = imag(Z_out)
%%aca no pasamos a paralelo(revisar notas del profe)

Z2=sqrt(R_out*Zo)
[w_out, imp_out, Lambda_p_out] = w_microstrip(e_r, H, t, Z2, f)
Largo_ac_out = Lambda_p_out/4

X_out_2 = -(Zo^2)/X_out; %pasamos de inductiva en serie a la entrada del acopllador a capacitiva en paralelo a la salida del acoplador
cap_out = 1/(2*pi*f*X_out_2)
beta = (2*pi)/Lambda_p_out%usamos lambda eficaz
d_cap_out = (1/beta)*acot(X_out_2/Zo)%Con impendacia caracteristica 50 ohm

%% Microstrip de 80 ohm

[W_80, Z_80, Lambda_p_80] = w_microstrip(e_r, H, t, 80, f)
QWT_80 = Lambda_p_80/4

%% Microstrip de 25 ohm

[W_25, Z_25, Lambda_p_25] = w_microstrip(e_r, H, t, 25, f)
QWT_25 = Lambda_p_25/4
