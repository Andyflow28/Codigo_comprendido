import numpy as np
import matplotlib.pyplot as plt

# Definir los archivos de datos
archivos = [
    "output_sin_desorden/output_sin_desorden_000E04t02.dat",
    "output_sin_desorden/output_sin_desorden_001E04t02.dat",
    "output_sin_desorden/output_sin_desorden_005E04t02.dat",
    "output_sin_desorden/output_sin_desorden_010E04t02.dat",
    "output_sin_desorden/output_sin_desorden_015E04t02.dat"
]

# Crear una figura para las gráficas
plt.figure(figsize=(10, 6))

# Colores y etiquetas para cada archivo
colores = ['b', 'g', 'r', 'c', 'm']
etiquetas = ['000E04t02', '001E04t02', '005E04t02', '010E04t02', '015E04t02']

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
