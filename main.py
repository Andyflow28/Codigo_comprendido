import matplotlib.pyplot as plt
import numpy as np

# Lista de archivos
file_paths = [
    "c0g_output/c0g001_output.dat",
    "c0g_output/c0g005_output.dat",
    "c0g_output/c0g010_output.dat",
    "c0g_output/c0g015_output.dat",
    "c0g_output/c0g020_output.dat"
]

# Iterar y graficar los datos
plt.figure(figsize=(10, 6))

for file_path in file_paths:
    # Cargar los datos
    data = np.loadtxt(file_path)
    x, y = data[:, 0], data[:, 1]
    
    # Graficar
    label = file_path.split("/")[-1]  # Usar el nombre del archivo como etiqueta
    plt.plot(x, y, label=label)

# Configurar la gráfica
plt.title("Gráfica de datos de múltiples archivos")
plt.xlabel("Eje X")
plt.ylabel("Eje Y")
plt.legend()
plt.grid()

# Guardar y mostrar
plt.savefig("multi_file_plot.png")
plt.show()
