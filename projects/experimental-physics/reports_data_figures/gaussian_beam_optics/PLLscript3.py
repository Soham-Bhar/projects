import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gaussian_beam(r, I0, rc, w):
    return I0 * np.exp(-2 * (r - rc)**2 / w**2)

# 1. Organize your data: { 'Distance Label': (r_data, V_data) }
experimental_data = {
    '2 m': (np.array([-1.6,-1.4,-1.2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6]), np.array([0.0035, 0.01075, 0.0395, 0.10875, 0.23825, 0.44725, 0.69425, 0.915, 1.075, 1.08025, 0.955, 0.721, 0.4775, 0.283, 0.14675, 0.069, 0.03275])),
    '5 m':   (np.array([0,0.6,1.2,1.8,2.4,3,3.6,4.2,4.8,5.4,6,6.6,7.2,7.8,8.4,9]), np.array([0.00525, 0.019, 0.0555, 0.14425, 0.31425, 0.5455, 0.80125, 0.90325, 0.84175, 0.6445, 0.404, 0.19375, 0.0655, 0.01625, 0.0035, 0.001])),
    '10 m':   (np.array([0,0.9,1.8,2.7,3.6,4.5,5.4,6.3,7.2,8.1,9,9.9,10.8,11.7,12.6,13.5]), np.array([0.0155, 0.0415, 0.091, 0.19425, 0.3555, 0.57175, 0.78525, 0.96725, 1.0825, 1.05575, 0.91925, 0.69975, 0.47825, 0.28425, 0.14575, 0.059])),
    '15 m':   (np.array([0,1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,11,12.1,13.2,14.3,15.4,16.3]), np.array([0.07725, 0.15325, 0.27975, 0.4385, 0.63625, 0.84025, 1.01625, 1.13575, 1.16525, 1.11025, 0.97425, 0.7845, 0.587, 0.4, 0.25125, 0.158])),
}

v_error = 0.01


colors = ['#1A5276', '#A93226', '#28B463', '#F5B041'] # Professional palette (Navy, Red, Green, Orange)
plt.figure(figsize=(10, 6))

# 2. Loop through each dataset to fit and plot
for i, (label, (r_data, V_data)) in enumerate(experimental_data.items()):
    popt, _ = curve_fit(gaussian_beam, r_data, V_data, p0=[np.max(V_data), 0, 1.0])
    I0_fit, rc_fit, w_fit = popt
    V_norm_err = v_error / I0_fit
    
    # 3. Normalization [cite: 811]
    V_norm = V_data / I0_fit
    r_fine = np.linspace(min(r_data), max(r_data), 200)
    V_fit_norm = gaussian_beam(r_fine, I0_fit, rc_fit, w_fit) / I0_fit
    
    # Plotting: Raw data as points with error bars [cite: 200, 202]
    # (Label removed to hide from legend for simplicity )
    plt.errorbar(r_data - rc_fit, V_norm, yerr=V_norm_err, fmt='o', 
                 color=colors[i], alpha=0.6, capsize=3) 
    
    # Plotting: Continuous line for the fit 
    # This remains as the single legend entry per distance 
    plt.plot(r_fine - rc_fit, V_fit_norm, color=colors[i], 
             label=f'Fit {label} ($W={w_fit:.2f}$ mm)')

# 4. Professional Formatting [cite: 181-183, 201]
plt.xlabel('Radial Position $r - r_c$ (mm)', fontsize=14)
plt.ylabel('Normalized Voltage ', fontsize=14)
plt.title('Evolution of Transverse Beam Profile with Distance', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.ylim(0, 1.1)
plt.show()