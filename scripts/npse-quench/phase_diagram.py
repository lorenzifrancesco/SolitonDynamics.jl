import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import re
import toml

def parse_file(filepath):
    """
    Reads a file with a specific format and parses the header (i, j indices) and associated data.

    Args:
        filepath (str): Path to the file.

    Returns:
        tuple: A dictionary of tuples (indices, data), where indices are (i, j) tuples and data is a numpy array,
               and a tuple of phase diagram bounds (v_0_range, a_s_range).
    """
    results = {}
    v_0_range, a_s_range = None, None

    with open(filepath, 'r') as file:
        for line in file:
            # Check for the header line containing the ranges
            if line.startswith("==="):
                match = re.search(r"v_0 in \[(.+?), (.+?)\], a_s in\[(.+?), (.+?)\]", line)
                if match:
                    v_0_range = (float(match.group(1)), float(match.group(2)))
                    a_s_range = (float(match.group(3)), float(match.group(4)))
            # Parse the indices and wavefunction data
            elif "|" in line:
                header, data = line.split("|")

                # Parse the indices (i, j)
                i, j = map(int, header.strip().split(","))

                # Parse the data values
                data_values = np.array(list(map(float, data.strip().split())))

                results[(i, j)] = data_values

    return results, (v_0_range, a_s_range)


def read_phase_diagram_config(config_path):
    """
    Reads the phase diagram configuration from a TOML file.

    Args:
        config_path (str): Path to the TOML configuration file.

    Returns:
        tuple: Grids for a_s and v_0 values.
    """
    config = toml.load(config_path)

    a_s_min = config["a_s_min"]
    a_s_max = config["a_s_max"]
    n_a_s = config["n_a_s"]
    v_0_min = config["v_0_min"]
    v_0_max = config["v_0_max"]
    n_v_0 = config["n_v_0"]

    a_s_values = np.linspace(a_s_min, a_s_max, n_a_s)
    v_0_values = np.linspace(v_0_min, v_0_max, n_v_0)

    return a_s_values, v_0_values

def save_wavefunction_plot(data, indices, output_folder):
    """
    Saves the plot of the wavefunction to a file.

    Args:
        data (np.array): The wavefunction data.
        indices (tuple): The (i, j) indices.
        output_folder (str): Path to the output folder.
    """
    ix, iy = indices
    filename = os.path.join(output_folder, f"name_{ix}_{iy}.png")

    plt.figure(figsize=(4, 3))
    plt.plot(data)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$|\psi|^2$")
    plt.title(f"GS of ({ix}, {iy})")
    plt.tight_layout() 
    plt.savefig(filename, dpi=300)
    plt.close()

def compute_peaks(data):
    """
    Computes the number of peaks in the wavefunction data.

    Args:
        data (np.array): The wavefunction data.

    Returns:
        int: The number of peaks.
    """
    peaks, _ = find_peaks(data)
    return len(peaks)

def generate_phase_diagram(peaks_data, a_s_values, v_0_values, output_file):
    """
    Generates and saves a phase diagram based on the number of peaks.

    Args:
        peaks_data (dict): A dictionary with keys as (i, j) and values as number of peaks.
        output_file (str): Path to save the phase diagram.
    """
    # Extract unique i and j values
    indices = np.array(list(peaks_data.keys()))
    values = np.array(list(peaks_data.values()))

    i_vals = np.unique(indices[:, 0])
    j_vals = np.unique(indices[:, 1])

    # Create a grid for the phase diagram
    phase_grid = np.zeros((len(i_vals), len(j_vals)))

    for (ix, iy), peaks in peaks_data.items():
        i_idx = np.where(i_vals == ix)[0][0]
        j_idx = np.where(j_vals == iy)[0][0]
        phase_grid[i_idx, j_idx] = peaks
        
    # Plot the phase diagram
    plt.figure(figsize=(4, 3.5))
    plt.imshow(phase_grid, origin="lower", aspect="auto", cmap="viridis", 
               extent=[a_s_values[0], a_s_values[-1], v_0_values[0], v_0_values[-1]])
    plt.colorbar(label=r"$\mathrm{Number \; of \; peaks}$")
    plt.xlabel(r"$a_s \quad [a_0]$")
    plt.ylabel(r"$V_0 \quad [E_r]$")
    # plt.title("Phase Diagram of Peaks")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

def fill_phase_diagram():
    # Filepath and output configuration
    filepath = "results/phase_diagram.txt"
    config_file = "input/phase_diagram.toml"
    output_folder = "media/gs"
    phase_diagram_file = os.path.join(output_folder, "phase_diagram.png")

    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Parse the file
    parsed_data, _ = parse_file(filepath)
    a_s_values, v_0_values = read_phase_diagram_config(config_file)

    # print(parsed_data[0])
    # Process wavefunctions and compute peaks
    print(parsed_data)
    peaks_data = {}
    for indices, data in parsed_data.items():
        save_wavefunction_plot(data, indices, output_folder)
        peaks_count = compute_peaks(data)
        print(f"Found {peaks_count:5d} peaks")
        peaks_data[indices] = peaks_count

    # Generate and save the phase diagram
    generate_phase_diagram(peaks_data, a_s_values, v_0_values, phase_diagram_file)

if __name__ == "__main__":
    fill_phase_diagram()
