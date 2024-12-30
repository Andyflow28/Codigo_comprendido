import numpy as np
import matplotlib.pyplot as plt

# Definir los archivos de datos
archivos = [
    "data_sin_desorde/c0g001_output_sin_desorde.dat",
    "data_sin_desorde/c0g005_output_sin_desorde.dat",
    "data_sin_desorde/c0g010_output_sin_desorde.dat",
    "data_sin_desorde/c0g015_output_sin_desorde.dat",
    "data_sin_desorde/c0g020_output_sin_desorde.dat"
]

# Crear una figura para las gráficas
plt.figure(figsize=(10, 6))

# Colores y etiquetas para cada archivo
colores = ['b', 'g', 'r', 'c', 'm']
etiquetas = ['c0g001', 'c0g005', 'c0g010', 'c0g015', 'c0g020']

# Cargar y graficar los datos de cada archivo
for archivo, color, etiqueta in zip(archivos, colores, etiquetas):
    # Cargar los datos desde el archivo
    datos = np.loadtxt(archivo)
    eje_x = datos[:, 0]  # Suponiendo que la primera columna es el eje X
    eje_y = datos[:, 1]  # Suponiendo que la segunda columna es el eje Y
    
    # Graficar los datos
    plt.plot(eje_x, eje_y, label=etiqueta, color=color)

# Configurar el gráfico
plt.xlabel('Energía (MeV)')
plt.ylabel('Valor')
plt.title('Gráfico de los archivos de salida sin desorden')
plt.legend()
plt.grid(True)

# Mostrar el gráfico
plt.show()
