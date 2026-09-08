import numpy as np
import matplotlib.pyplot as plt
from scipy.special import airy, ai_zeros

# Define Airy parameters
n_states = 5  # number of eigenstates to plot
z = np.linspace(0, 20, 500)

# Get first few zeros of the Airy Ai function
# ai_zeros returns the negative roots of Ai(x) [2, 3, 4, 5, 6, 7]
zeros = ai_zeros(n_states)  
# The scaled energies E_n' are positive 
E_n_prime = -zeros  

# Compute eigenfunctions ψ_n(z')
psi =[]
for n in range(n_states):
    # airy(x) returns the Ai(x) function [8]
    psi_n = airy(z - E_n_prime[n])  # Ai(z - E'_n)
    # Normalize the eigenfunction using trapezoidal integration
    psi.append(psi_n / np.sqrt(np.trapz(psi_n**2, z)))  

# Plot results
plt.figure(figsize=(8, 6))

# Add a horizontal line at y=0 to highlight the zeros
plt.axhline(0, color='black', linestyle='--', lw=1)

for n in range(n_states):
    # Plot without the y-axis shift
    plt.plot(z, psi[n], label=f"$n={n+1}$, $E'_n={E_n_prime[n]:.3f}$")

plt.title("Solutions to the Airy Equation $\\psi'' - (z' - E'_n)\\psi = 0$")
plt.xlabel("$z'$ (dimensionless height)")
plt.ylabel("$\\psi_n(z')$") # Updated label
plt.legend()
plt.grid(True)
plt.show()