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
    
    return len(zero_crossings)

l_perp = 1.59e-6
# Define parameter ranges
n_points = 50
g_values = np.linspace(-1.6, 0.0, n_points)
kL_values = np.linspace(0, 10/(l_perp*1e6), n_points)
V0 = 1.3 * 1.386 # Er for dL/l_perp = 1.887
heatmap_data = np.zeros((len(kL_values), len(g_values)))

if True:
  # Compute number of zeros for each (g, V0) pair
  for i, kL in enumerate(kL_values):
    for j, g in enumerate(g_values):
      heatmap_data[i, j] = count_zeros(g, V0, kL)
  np.save("results/phase_diagram/variational.npy", heatmap_data)

heatmap_data = np.load("results/phase_diagram/variational_kl.npy")

print(f"Plotting the {heatmap_data.shape} phase diagram...")
plt.figure(figsize=(3.7, 3))
ax = sns.heatmap(heatmap_data[::-1, :], 
            # yticks = [0, n_points-1],
            yticklabels=[kL_values[0], kL_values[-1]],
            # xticks=[0, n_points-1],
            xticklabels=[g_values[0], g_values[-1]], cmap='viridis')
ytick_positions = [0, len(kL_values) - 1][::-1]
xtick_positions = [0, len(g_values) - 1]
ytick_labels = [round(kL_values[0], 2), round(kL_values[-1], 2)]
xtick_labels = [round(g_values[0], 2), round(g_values[-1], 2)]
ax.set_xticks(xtick_positions)
ax.set_xticklabels(xtick_labels)
ax.set_yticks(ytick_positions)
ax.set_yticklabels(ytick_labels)

plt.xlabel(r'$g$')
plt.ylabel(r'$k_L$')
plt.tight_layout()
plt.savefig('media/variational_kl_phase_diagram.png')
