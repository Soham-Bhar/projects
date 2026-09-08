import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1. Define the Gaussian Model
# I(r) = I0 * exp(-2 * (r - rc)^2 / w^2)
def gaussian_beam(r, I0, rc, w):
    return I0 * np.exp(-2 * (r - rc)**2 / w**2)

# 2. Input your experimental data here
# Replace these with your actual measured values for a specific distance (z)
r_data = np.array([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]) # Horizontal pos (mm)
V_data = np.array([0.12, 0.35, 0.72, 0.95, 1.02, 0.98, 0.75, 0.38, 0.15]) # Voltage (V)

# 3. Initial Guesses (Crucial for non-linear fitting)
I0_guess = np.max(V_data)
rc_guess = r_data[np.argmax(V_data)]
w_guess = 1.0 
initial_guesses = [I0_guess, rc_guess, w_guess]

# 4. Perform the Fit
popt, pcov = curve_fit(gaussian_beam, r_data, V_data, p0=initial_guesses)
I0_fit, rc_fit, w_fit = popt
errors = np.sqrt(np.diag(pcov)) # Standard deviation of parameters

# 5. Output Results
print(f"--- Fit Results ---")
print(f"Calculated Peak (I0): {I0_fit:.4f} ± {errors[0]:.4f} V")
print(f"Center Offset (rc):   {rc_fit:.4f} ± {errors[1]:.4f} mm")
print(f"Spot Size (W):        {w_fit:.4f} ± {errors[2]:.4f} mm")

# 6. Visualization (Professional format for your slides)
r_fine = np.linspace(min(r_data), max(r_data), 100)
V_fit = gaussian_beam(r_fine, *popt)

plt.figure(figsize=(8, 5))
plt.scatter(r_data, V_data, color='red', label='Experimental Data')
plt.plot(r_fine, V_fit, color='navy', label=f'Gaussian Fit (W={w_fit:.2f}mm)')
plt.xlabel('Horizontal Position r (mm)')
plt.ylabel('Intensity / Voltage (V)')
plt.title('Transverse Beam Profile Characterization')
plt.legend()
plt.grid(alpha=0.3)
plt.show()

