import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian_beam(r, I0, rc, w):
    return I0 * np.exp(-2 * (r - rc)**2 / w**2)

experimental_data = {
    '10 cm': (np.array([0,0.07,0.14,0.21,0.28,0.35,0.42,0.49,0.56,0.63,0.7,0.77,0.84,0.91,0.98,1.05,1.12
]), np.array([0.0105,0.027,0.0665,0.1765,0.27,0.41,0.6165,0.825,0.947,1.118,1.0955,0.8295,0.5765,0.366,0.156,0.0595,0.019
])),
    '20 cm': (np.array([-1.05,-0.95,-0.85,-0.75,-0.65,-0.55,-0.45,-0.35,-0.25,-0.15,-0.05,0.05,0.15,0.25,0.35]), np.array([0.006,0.00925,0.03375,0.1295,0.2955,0.62025,0.9325,1.23325,1.15025,0.81725,0.383,0.17525,0.05325,0.012,0.0075])),
    '50 cm': (np.array([0.75,0.63,0.51,0.39,0.27,0.15,0.03,-0.09,-0.21,-0.33,-0.45,-0.57,-0.69,-0.81,-0.93]), np.array([0.004,0.01575,0.0605,0.201,0.4725,0.8695,1.15125,1.17575,0.96275,0.6675,0.39125,0.137,0.03225,0.01425,0.0035])),
    '1 m': (np.array([-0.1,0.04,0.18,0.32,0.46,0.6,0.74,0.88,1.02,1.16,1.3,1.44,1.58,1.72,1.86]), np.array([0.01975,0.068,0.16875,0.3605,0.6235,0.92,1.12325,1.10225,0.9035,0.63425,0.38375,0.186,0.0815,0.032,0.012])),
}

v_error = 0.01
colors = ['#1A5276', '#A93226', '#28B463', '#F5B041']
plt.figure(figsize=(10, 6))

for i, (label, (r_data, V_data)) in enumerate(experimental_data.items()):
    # IMPROVEMENT: Use the peak of the data to find a better initial guess for the center (rc) 
    rc_initial = r_data[np.argmax(V_data)]
    
    # Perform the fit 
    popt, _ = curve_fit(gaussian_beam, r_data, V_data, p0=[np.max(V_data), rc_initial, 1.0])
    I0_fit, rc_fit, w_fit = popt
    
    # IMPROVEMENT: Use abs() to ensure the spot size W is always displayed as a positive radius
    w_fit = abs(w_fit)
    
    V_norm_err = v_error / I0_fit
    V_norm = V_data / I0_fit
    r_fine = np.linspace(min(r_data), max(r_data), 200)
    V_fit_norm = gaussian_beam(r_fine, I0_fit, rc_fit, w_fit) / I0_fit
    
    # Plotting: Raw data [cite: 200, 202]
    plt.errorbar(r_data - rc_fit, V_norm, yerr=V_norm_err, fmt='o', 
                 color=colors[i], alpha=0.6, capsize=3) 
    
    # Plotting: Fit [cite: 201]
    plt.plot(r_fine - rc_fit, V_fit_norm, color=colors[i], 
             label=f'Fit {label} ($W={w_fit:.2f}$ mm)')

# Formatting [cite: 181-183, 201]
plt.xlabel('Radial Position $r - r_c$ (mm)', fontsize=14)
plt.ylabel('Normalized Voltage', fontsize=14)
plt.title('Evolution of Transverse Beam Profile with Distance', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.ylim(0, 1.1)
plt.show()