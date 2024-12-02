import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.special import genlaguerre
import seaborn as sns
import matplotlib.animation as animation

def plot_heatmap():
  # Set up matplotlib for LaTeX-compatible fonts
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')

  # Load the CSV file
  filename = "results/experiment1.csv"  # Replace with your CSV file path
  data = pd.read_csv(filename)
  time = data.iloc[0, :].values  # First column is time
  print("\033[91mWarn:\033[0m playing sketchy stuff with data extraction")
  psi2_values = data.iloc[1:, :].values  # Remaining columns are ψ² values

  # Heatmap plot
  fig, ax = plt.subplots(figsize=(4, 3))
  sns.heatmap(
      psi2_values,  # Transpose so time is along y-axis
      xticklabels=False,  # Suppress x-ticks
      yticklabels=False,  # Suppress y-ticks
      cmap="viridis",  # Choose a colormap
      cbar_kws={'label': r'$|\psi|^2$'},  # LaTeX for colorbar label
      ax=ax
  )
  ax.set_xlabel(r'$t$')
  ax.set_ylabel(r'$x$')
  plt.tight_layout()

  # Save the heatmap as a PDF
  heatmap_filename = "media/td_heatmap.png"
  plt.savefig(heatmap_filename, dpi=300)
  plt.close()
  print(f"Heatmap saved as {heatmap_filename}")



def plot_animation():
  # Set up matplotlib for LaTeX-compatible fonts
  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')

  # Load the CSV file
  filename = "results/experiment1.csv"  # Replace with your CSV file path
  data = pd.read_csv(filename)
  time = data.iloc[0, :].values  # First row is time
  print("\033[91mWarn:\033[0m playing sketchy stuff with data extraction")
  psi2_values = data.iloc[1:, :].values  # Remaining rows are ψ² values
  print(np.shape(psi2_values))

  # Set up the plot for animation
  fig, ax = plt.subplots(figsize=(4, 3))
  ax.set_xlabel(r'$x$')
  ax.set_ylabel(r'$|\psi|^2$')
  # ax.set_xlim(min(time), max(time))
  ax.set_ylim(0, psi2_values.max())

  # Initialize an empty line plot
  space = range(len(psi2_values[:, 1]))
  line, = ax.plot(space, psi2_values[:, 0], lw=1.5)
  
  # Function to update the line plot for each time frame
  def update_frame(i):
      x = space
      y = psi2_values[:, i]
      line.set_data(x, y)  # Update the line with new data
      ax.set_title(f"Time: {time[i]:.2f} $2\pi\omega_\perp^{-1}$")
      return [line]
    
  # Create the animation
  ani = animation.FuncAnimation(fig, update_frame, frames=len(time), interval=10, blit=True)
  fig.tight_layout()
  # Save the animation as a GIF
  gif_filename = "media/td_line_animation.gif"
  ani.save(gif_filename, writer='pillow', dpi=100)
  print(f"Animation saved as {gif_filename}")
  plt.close()
