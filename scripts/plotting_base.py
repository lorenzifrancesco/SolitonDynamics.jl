import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.special import genlaguerre

# Function to calculate the radial part of the Laguerre-Gaussian mode
def laguerre_gauss_radial(r, p, l, w0):
    normalization = np.sqrt(2 * np.math.factorial(p) / (np.pi * np.math.factorial(p + abs(l))))
    radial_term = (r * np.sqrt(2) / w0) ** abs(l)
    gaussian_envelope = np.exp(-r**2 / w0**2)
    laguerre_poly = genlaguerre(p, abs(l))(2 * r**2 / w0**2)
    return normalization * radial_term * gaussian_envelope * laguerre_poly

def plot_axial():
  data = pd.read_csv("results/data.csv")
  x = data["x"]
  y = data["y"]
  plt.figure(figsize=(3, 2.4), dpi=300)
  plt.plot(x, y, linestyle="-", color="blue", linewidth=1.1)
  plt.xlabel(r'$t$')
  plt.ylabel(r'$|A|^2$')
  plt.grid(False)
  plt.tight_layout()
  plt.savefig("media/axial.pdf", format="pdf", dpi=300)
  
def plot_radial():
    # Parameters
    l = 1        # Azimuthal index (topological charge)
    r = np.linspace(0, 3, 500)  # Radial coordinate

    p_color = ["r", "b"]
    labels = [r'$p=0$, $m=1$', r'$p=1$, $m=1$']
    plt.figure(figsize=(3, 2.4))
    for p in [0, 1]:
      # Different rescaling factors for w0
      w0_values = [0.8, 1]
      w0_styles = ["-", "--"]
      for iw, w0 in enumerate(w0_values):
          radial_profile = laguerre_gauss_radial(r, p, l, w0)
          if iw == 0:
            plt.plot(r, np.abs(radial_profile)**2,
                    ls=w0_styles[iw], 
                    color=p_color[p], 
                    label=labels[p])
          else:
            plt.plot(r, np.abs(radial_profile)**2,
                    ls=w0_styles[iw], 
                    color=p_color[p])
    plt.xlabel(r'$r$')
    plt.ylabel(r'$|T|^2$')
    plt.grid(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig("media/radial.pdf", format="pdf", dpi=300)