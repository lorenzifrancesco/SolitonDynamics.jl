import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import brentq

def my_function(eta, g, V0, kL):
    if eta <= 0:
        return np.nan
    
    exp_term = np.exp(kL**2 * eta**2)
    numerator = (-8 * exp_term * kL**2 * V0 * eta**3 * (g * np.sqrt(2 * np.pi) + 2 * np.pi * eta)
                 + 16 * kL**4 * V0**2 * eta**6 * (g * np.sqrt(2 * np.pi) + 2 * np.pi * eta)
                 + np.exp(2 * kL**2 * eta**2) * (2 * np.pi * eta + g * (np.sqrt(2 * np.pi) - g * eta**5)))
    
    denominator = g * np.sqrt((-1 + 4 * np.exp(-kL**2 * eta**2) * kL**2 * V0 * eta**3) / g)
    
    if denominator == 0:
        return np.nan
    
    return numerator / denominator

def count_zeros(g, V0, kL, eta_min=1e-4, eta_max=15, num_points=10000):
    eta_vals = np.linspace(eta_min, eta_max, num_points)
    y_vals = [my_function(eta, g, V0, kL) for eta in eta_vals]
    
    zero_crossings = []
    for i in range(len(eta_vals) - 1):
        if np.isnan(y_vals[i]) or np.isnan(y_vals[i+1]):
            continue  # Skip invalid values
        if y_vals[i] * y_vals[i+1] < 0:
            try:
                root = brentq(my_function, eta_vals[i], eta_vals[i+1], args=(g, V0, kL))
                if not any(np.isclose(root, z, atol=1e-3) for z in zero_crossings):
                    zero_crossings.append(root)
            except ValueError:
                pass
    if len(zero_crossings) == 2 and zero_crossings[1]<1:
      return 5
    else:
      return len(zero_crossings)

# Define parameter ranges
n_points = 50
g_values = np.linspace(-1.6, -0.2, n_points)
V0_values = np.linspace(0, 3, n_points)
V0_values = V0_values * 1.386 # Er for dL/l_perp = 1.887
kL = 2
heatmap_data = np.zeros((len(V0_values), len(g_values)))

if False:
  # Compute number of zeros for each (g, V0) pair
  for i, V0 in enumerate(V0_values):
    for j, g in enumerate(g_values):
      heatmap_data[i, j] = count_zeros(g, V0, kL)
  np.save("results/phase_diagram/variational_v0.npy", heatmap_data)

heatmap_data = np.load("results/phase_diagram/variational_v0.npy")

V0_values /= 1.386 # go back to units of Er
print(f"Plotting the {heatmap_data.shape} phase diagram...")
plt.figure(figsize=(3.7, 3))
ax = sns.heatmap(heatmap_data[::-1, :], 
            # yticks = [0, n_points-1],
            yticklabels=[V0_values[0], V0_values[-1]],
            # xticks=[0, n_points-1],
            xticklabels=[g_values[0], g_values[-1]], cmap='viridis')
ytick_positions = [0, len(V0_values) - 1]
xtick_positions = [0, len(g_values) - 1]
ytick_labels = [round(V0_values[0], 2), round(V0_values[-1], 2)][::-1]
xtick_labels = [round(g_values[0], 2), round(g_values[-1], 2)]
ax.set_xticks(xtick_positions)
ax.set_xticklabels(xtick_labels)
ax.set_yticks(ytick_positions)
ax.set_yticklabels(ytick_labels)

plt.xlabel(r'$g$')
plt.ylabel(r'$V_0$ [$E_r$]')
plt.tight_layout()
plt.savefig('media/variational_v0_phase_diagram.png', dpi=300)
