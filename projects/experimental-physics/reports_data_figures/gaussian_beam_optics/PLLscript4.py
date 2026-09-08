import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1. Gaussian Spot Size Evolution Equation (Eq. 2 from your manual)
# w(z) = w0 * sqrt(1 + (z/zr)^2)
def w_model(z, w0, zr):
    return w0 * np.sqrt(1 + (z / zr)**2)


z_data = np.array([0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0]) 
w_data = np.array([0.41, 0.38, 0.49, 0.65, 1.06, 2.61, 5.16, 7.65]) 

# 3. Perform the Non-linear Fit
# Initial guesses: w0 = first measured W, zr = 0.7m
popt, pcov = curve_fit(w_model, z_data, w_data, p0=[w_data[0], 0.7])
w0_fit, zr_fit = popt

# 4. Professional Visualization for your deck 
z_plot = np.linspace(0, max(z_data), 500)
plt.figure(figsize=(10, 6))
plt.scatter(z_data, w_data, color='black', label='Determined Spot Sizes ($W$)')
plt.plot(z_plot, w_model(z_plot, *popt), color='red', linewidth=2, label='Fit according to spot size evolution ')

plt.xlabel('Propagation Distance $z$ (m)', fontsize=14)
plt.ylabel('Spot Size $W(z)$ (mm)', fontsize=14)
plt.title('Spot Size vs. Distance', fontsize=16)
plt.legend()
plt.grid(alpha=0.3)
plt.show()

print(f"Final Beam Waist (W0): {w0_fit:.4f} mm")
print(f"Final Rayleigh Range (Zr): {zr_fit:.4f} m")