import pandas as pd
import matplotlib.pyplot as plt

file_path = './output/pi_plus.dat'

# Read the data
data = pd.read_csv(file_path, delimiter='\t', comment='#')

# Plot dNd3p against pT
plt.figure(figsize=(10, 6))
plt.plot(data['pT'], data['dNd3p'], label='dNd3p vs pT')
plt.xlabel('pT')
plt.ylabel('dNd3p')
plt.title('dNd3p vs pT')
plt.legend()
plt.grid(True)
plt.show()

# Plot dNd3p against phi_p
plt.figure(figsize=(10, 6))
plt.plot(data['phi_p'], data['dNd3p'], label='dNd3p vs phi_p')
plt.xlabel('phi_p')
plt.ylabel('dNd3p')
plt.title('dNd3p vs phi_p')
plt.legend()
plt.grid(True)
plt.show()

# Plot dNd3p against y_p
plt.figure(figsize=(10, 6))
plt.plot(data['y_p'], data['dNd3p'], label='dNd3p vs y_p')
plt.xlabel('y_p')
plt.ylabel('dNd3p')
plt.title('dNd3p vs y_p')
plt.legend()
plt.grid(True)
plt.show()
