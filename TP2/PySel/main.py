import os
import skrf
import math
import numpy as np

def calcular_parametros_s(frecuencia, carpetas):
    parametros_s = {}
    
    # Convertir frecuencia de GHz a Hz
    frecuencia *= 1e9
    
    for carpeta in carpetas:
        archivos_s2p = [f for f in os.listdir(carpeta) if f.endswith('.s2p')]
        for archivo_s2p in archivos_s2p:
            ruta_archivo = os.path.join(carpeta, archivo_s2p)
            red = skrf.network.Network(ruta_archivo)
            
            print(f"Archivo: {archivo_s2p}")
            
            # Encontrar la frecuencia más cercana en el archivo .s2p
            indice_frecuencia = np.argmin(np.abs(red.f - frecuencia))
            
            print(f"Frecuencia más cercana: {red.f[indice_frecuencia]}")
            
            # Acceder a los parámetros S en la frecuencia encontrada
            if archivo_s2p not in parametros_s:
                parametros_s[archivo_s2p] = red.s[indice_frecuencia]
            else:
                print(f"Advertencia: Frecuencia {frecuencia} duplicada en {archivo_s2p}")
    
    return parametros_s



def calcular_delta_y_guardar(parametros_s):
    deltas = {}
    polarizaciones_buenas_delta = {}
    for archivo, parametros in parametros_s.items():
        delta = (parametros[0, 0] * parametros[1, 1]) - (parametros[0, 1] * parametros[1, 0])
        deltas[archivo] = delta
        if abs(delta) < 1:

            polarizaciones_buenas_delta[archivo] = {'delta': delta}

    return deltas, polarizaciones_buenas_delta



def calcular_k_y_guardar(parametros_s, deltas):
    ks = {}
    polarizaciones_buenas_k = {}
    for archivo, parametros in parametros_s.items():
        S11 = parametros[0, 0]
        S12 = parametros[0, 1]
        S21 = parametros[1, 0]
        S22 = parametros[1, 1]
        abs_delta = abs(deltas[archivo])
        
        # Calcular el parámetro k
        numerador = (1 - abs(S11)**2 - abs(S22)**2 + abs_delta**2)
        denominador = 2 * abs(S12 * S21)
        k = numerador / denominador
        
        ks[archivo] = k
        if k > 1:
            # Guardar la polarización si k es mayor que 1
            polarizaciones_buenas_k[archivo] = {'k': k}
    
    return ks, polarizaciones_buenas_k



def guardar_polarizaciones_estables(polarizaciones_buenas_delta, polarizaciones_buenas_k):
    polarizaciones_estables = {}
    for archivo, info in polarizaciones_buenas_delta.items():
        if archivo in polarizaciones_buenas_k:
            polarizaciones_estables[archivo] = info
    return polarizaciones_estables



def seleccionar_polarizacion_estable(polarizaciones_estables):
    print("Polarizaciones estables disponibles:")
    for i, archivo in enumerate(polarizaciones_estables.keys(), start=1):
        print(f"{i}. {archivo}")
    
    seleccion = int(input("Seleccione una polarización estable: ")) - 1
    archivos = list(polarizaciones_estables.keys())
    archivo_seleccionado = archivos[seleccion]

    return archivo_seleccionado


def calcular_coeficientes_de_reflexion(polarizaciones_estables, parametros_s, deltas):
    coeficientes = {}
    
    for archivo_seleccionado in polarizaciones_estables:
        parametros_seleccionados = parametros_s[archivo_seleccionado]
        delta = deltas[archivo_seleccionado]

        S11 = parametros_seleccionados[0, 0]
        S22 = parametros_seleccionados[1, 1]

        B1 = 1 + abs(S11)**2 - abs(S22)**2 - abs(delta)**2
        B2 = 1 + abs(S22)**2 - abs(S11)**2 - abs(delta)**2
        C1 = S11 - (delta * np.conj(S22))
        C2 = S22 - (delta * np.conj(S11))

        if B1 > 0:
            r_Ms = (B1 - np.sqrt(B1**2 - 4 * abs(C1)**2)) / (2 * C1)
        else:
            r_Ms = (B1 + np.sqrt(B1**2 - 4 * abs(C1)**2)) / (2 * C1)

        if B2 > 0:
            r_ML = (B2 - np.sqrt(B2**2 - 4 * abs(C2)**2)) / (2 * C2)
        else:
            r_ML = (B2 + np.sqrt(B2**2 - 4 * abs(C2)**2)) / (2 * C2)

        r_in = np.conj(r_Ms)
        r_out = np.conj(r_ML)

        coeficientes[archivo_seleccionado] = {'r_in': r_in, 'r_out': r_out}

    return coeficientes


def calcular_impedancias(parametros_s, deltas, coeficientes_reflexion, Zo):
    impedancias = {}
    for archivo_seleccionado, coeficientes in coeficientes_reflexion.items():
        r_in = coeficientes['r_in']
        r_out = coeficientes['r_out']
        
        # Calcular Z_in y Z_out
        Z_in = Zo * (1 + r_in) / (1 - r_in)
        Z_out = Zo * (1 + r_out) / (1 - r_out)
        
        # Calcular Z_s y Z_L (conjugados de Z_in y Z_out respectivamente)
        Z_s = np.conj(Z_in)
        Z_L = np.conj(Z_out)
        
        # Guardar las impedancias en el diccionario
        impedancias[archivo_seleccionado] = {'Z_in': Z_in, 'Z_out': Z_out, 'Z_s': Z_s, 'Z_L': Z_L}
    
    return impedancias

def calcular_impedancias_paralelo(impedancias):
    impedancias_paralelo = {}
    for archivo, imp in impedancias.items():
        Z_in = imp['Z_in']
        R_in = Z_in.real
        X_in = Z_in.imag
        R_inp = R_in * (1 + (X_in/R_in)**2)
        X_inp = X_in * (R_inp/X_in)
        impedancias_paralelo[archivo] = {'R': R_inp, 'X': X_inp}
    return impedancias_paralelo

def calcular_Z0_microstrip_in(impedancias_paralelo, Zo):
    Z0_microstrip_in = {}
    for archivo, imp in impedancias_paralelo.items():
        R_inp = imp['R']
        Z0_microstrip_in[archivo] = math.sqrt(R_inp * Zo)
    return Z0_microstrip_in

def calcular_microstrip(e_r, H, Z0_microstrip_in):
    t = 0.05  # grosor de la placa en mm

    W_list = []
    We_list = []
    e_rp_list = []
    Z0_list = []

    for Zo in Z0_microstrip_in.values():
        A = (Zo/60)*math.sqrt((e_r+1)/2)+((e_r-1)/(e_r+1))*(0.23+(0.11/e_r))
        B = (377*math.pi)/(2*Zo*math.sqrt(e_r))

        W_H = (8*math.exp(A))/(math.exp(2*A)-2)

        if W_H > 2:
            W_H = (2/math.pi)*(B-1-math.log(2*B-1))+((e_r-1)/(2*e_r))*(math.log(B-1)+0.39-(0.61/e_r))

        W = W_H*H
        W_list.append(W)

        if W_H <= (1/(2*math.pi)):
            We = W + (t/math.pi)*(1+math.log((4*math.pi*W)/t))
        else:
            We = W + (t/math.pi)*(1+math.log((2*H)/t))
        We_list.append(We)

        if W_H >= 1:
            e_rp = ((e_r+1)/2) + ((e_r-1)/2)*(1/(math.sqrt(1+(12*H)/W)))
            Z0 = ((120*math.pi)/math.sqrt(e_rp))/(W_H+1.393+0.667*math.log(W_H+1.444))
        else:
            e_rp = ((e_r+1)/2) + ((e_r-1)/2)*((1/(math.sqrt(1+(12*H)/W)))+0.004*(1-W_H)**2)
            Z0 = (60/math.sqrt(e_rp))*math.log(((8*H)/W)+(W/(4*H)))
        e_rp_list.append(e_rp)
        Z0_list.append(Z0)

    return W_list, We_list, e_rp_list, Z0_list
'''

----------------------------------------MAIN-----------------------------------------------

'''

# Directorio que contiene los archivos .s2p
directorio_transistor = input("Ingrese modelo del transitor: ")
while directorio_transistor not in ["BFP450", "BFP640"]:
    directorio_transistor = input("Ingrese el directorio del transistor (BFP450 o BFP640): ")

# Frecuencia seleccionada por el usuario en GHz
frecuencia_usuario = float(input("Ingrese Frecuencia en GHz: "))  # Por ejemplo, 0.03 GHz

# Calcular parámetros S para todas las polarizaciones
parametros_s = calcular_parametros_s(frecuencia_usuario, [directorio_transistor])
# Mostrar los parámetros S en forma polar
for archivo, parametros in parametros_s.items():
    print(f"Parámetros S para {archivo} en {frecuencia_usuario} GHz:")
    print("S11:", abs(parametros[0, 0]), "∠", np.angle(parametros[0, 0], deg=True), "grados")
    print("S12:", abs(parametros[0, 1]), "∠", np.angle(parametros[0, 1], deg=True), "grados")
    print("S21:", abs(parametros[1, 0]), "∠", np.angle(parametros[1, 0], deg=True), "grados")
    print("S22:", abs(parametros[1, 1]), "∠", np.angle(parametros[1, 1], deg=True), "grados")
    print()

#Cálculo del delta para todas las polarizaciones y guardar las polarizaciones compatibles
deltas, polarizaciones_buenas_delta = calcular_delta_y_guardar(parametros_s)
# Mostrar el parámetro delta para cada polarización y verificar si es buena
for archivo, delta in deltas.items():
    print(f"Parámetro delta para {archivo}: {delta}")
    if abs(delta) < 1:
        print("¡La polarización es buena!")
    else:
        print("La polarización no es buena.")
    print()

# Calcular el parámetro k y guardar las polarizaciones compatibles
ks, polarizaciones_buenas_k = calcular_k_y_guardar(parametros_s, deltas)
# Mostrar el parámetro k para cada polarización y verificar si es buena
for archivo, k in ks.items():
    print(f"Parámetro k para {archivo}: {k}")
    if k > 1:
        print("¡La polarización es buena!")
        print("k:", polarizaciones_buenas_k[archivo]['k'])
    else:
        print("La polarización no es buena.")
    print()


# Guardar las polarizaciones estables
polarizaciones_estables = guardar_polarizaciones_estables(polarizaciones_buenas_delta, polarizaciones_buenas_k)
# Mostrar las polarizaciones estables
print(f'Polarizaciones estables para {frecuencia_usuario} GHz:')
i=0
for archivo, info in polarizaciones_estables.items():
    i+=1
    print(f"{i}- {archivo}")

# Seleccionar una polarización estable
#archivo_seleccionado = seleccionar_polarizacion_estable(polarizaciones_estables)
#print("Polarización elegida:", archivo_seleccionado)

# Calcular r_in y r_out para las polarizaciones incondicionalmente estables
coeficientes_reflexion = calcular_coeficientes_de_reflexion(polarizaciones_estables, parametros_s, deltas)
# Mostrar los valores de r_in y r_out
for archivo_seleccionado, coeficientes in coeficientes_reflexion.items():
    print(f"Polarización: {archivo_seleccionado}")
    print(f"r_in: {coeficientes['r_in']}")
    print(f"r_out: {coeficientes['r_out']}")
    print()

# Calcular las impedancias para cada polarización estable
Zo = 50
impedancias = calcular_impedancias(parametros_s, deltas, coeficientes_reflexion, Zo)
# Mostrar las impedancias para cada polarización estable
print("Impedancias para cada polarización estable:")
for archivo_seleccionado, valores_impedancia in impedancias.items():
    print(f"Polarización: {archivo_seleccionado}")
    print(f"Z_in: {valores_impedancia['Z_in']}")
    print(f"Z_out: {valores_impedancia['Z_out']}")
    print(f"Z_s: {valores_impedancia['Z_s']}")
    print(f"Z_L: {valores_impedancia['Z_L']}")
    print()

#Calcular las impedancias de la función de arriba pero en formato paralelo
impedancias_paralelo = calcular_impedancias_paralelo(impedancias)


#Cálculo de MicroStip para acoplamiento de entrada:
Z0_microstrip_in = calcular_Z0_microstrip_in(impedancias_paralelo, Zo)
#Cálculo de MicroStip para acoplamientos de entrada y salida:
print("Ingrese permitividad relativa del sustrato y altura de la placa: \n")
e_r = float(input("Ingrese permitividad relativa del sustrato: "))
H = float(input("Ingrese altura de la placa en mm: "))
W_list, We_list, e_rp_list, Z0_list = calcular_microstrip(e_r, H, Z0_microstrip_in)
# Imprimir los valores de W_list y los nombres de los archivos correspondientes
for We, archivo in zip(We_list, polarizaciones_estables.keys()):
    print(f"Polarización: {archivo}")
    print(f"W: {We}")