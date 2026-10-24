#!/usr/bin/env python3
"""
Plotting tool for rotational evolution data
Plots angular velocity from wconvevol.dat and wradevol.dat on log-log scale
"""

import os

import matplotlib.pyplot as plt
import numpy as np


def plot_rotation_evolution(output_dir="output", save=True):
    """
    Plot rotation evolution: wconv and wrad vs time on log-log scale

    Parameters:
    -----------
    output_dir : str
        Directory containing output data files (default: 'output')
    save : bool
        Whether to save the figure (default: True)
    """

    print("\n" + "=" * 80)
    print("ROTATIONAL EVOLUTION PLOT")
    print("=" * 80)

    # File paths
    wconv_file = os.path.join(output_dir, "wconvevol.dat")
    wrad_file = os.path.join(output_dir, "wradevol.dat")

    # Check if files exist
    if not os.path.exists(wconv_file):
        print(f"ERROR: {wconv_file} not found!")
        return
    if not os.path.exists(wrad_file):
        print(f"ERROR: {wrad_file} not found!")
        return

    print(f"Reading data from: {output_dir}/")

    # Load data
    try:
        data_conv = np.loadtxt(wconv_file)
        data_rad = np.loadtxt(wrad_file)

        # Handle single-line data
        if data_conv.ndim == 1:
            data_conv = data_conv.reshape(1, -1)
        if data_rad.ndim == 1:
            data_rad = data_rad.reshape(1, -1)

        # Extract columns
        time_conv = data_conv[:, 0]  # Time in Myr
        w_conv = data_conv[:, 1]  # Angular velocity in Wsol

        time_rad = data_rad[:, 0]  # Time in Myr
        w_rad = data_rad[:, 1]  # Angular velocity in Wsol

        print(f"✓ Loaded wconvevol.dat: {len(time_conv)} points")
        print(f"  Time: {time_conv[0]:.3e} - {time_conv[-1]:.3e} Myr")
        print(f"  Wconv: {np.min(w_conv):.3e} - {np.max(w_conv):.3e} Wsol")
        print(f"✓ Loaded wradevol.dat: {len(time_rad)} points")
        print(f"  Time: {time_rad[0]:.3e} - {time_rad[-1]:.3e} Myr")
        print(f"  Wrad: {np.min(w_rad):.3e} - {np.max(w_rad):.3e} Wsol")

    except Exception as e:
        print(f"ERROR reading data files: {e}")
        return

    # Create figure
    print("\nCreating figure...")
    fig, ax = plt.subplots(figsize=(12, 8), dpi=150)

    # Plot envelope (plain line)
    ax.loglog(
        time_conv,
        w_conv,
        "b-",
        linewidth=2.5,
        label="Envelope (Wconv)",
        marker="o",
        markersize=2.5,
        markevery=15,
        alpha=0.85,
    )

    # Plot core (dashed line)
    ax.loglog(
        time_rad, w_rad, "b--", linewidth=2.5, label="Core (Wrad)", marker="s", markersize=2.5, markevery=15, alpha=0.85
    )

    # Format axes
    ax.set_xlabel("Time (Myr)", fontsize=13, fontweight="bold")
    ax.set_ylabel("Angular Velocity (Wsol)", fontsize=13, fontweight="bold")
    ax.set_title("Rotational Evolution: Envelope vs Core (Log-Log Scale)", fontsize=14, fontweight="bold", pad=15)

    # Grid
    ax.grid(True, which="both", linestyle="--", alpha=0.3, linewidth=0.8)
    ax.grid(True, which="major", linestyle="-", alpha=0.5, linewidth=1.0)

    # Legend
    ax.legend(fontsize=12, loc="best", framealpha=0.95, shadow=True, edgecolor="black")

    # Tight layout
    ax.set_xlim(1e0, 5e3)
    ax.set_ylim(5e-1, 1e2)
    plt.tight_layout()

    print("✓ Figure created")

    # Save figure
    if save:
        output_file = os.path.join(output_dir, "evolution_plot.png")
        try:
            plt.savefig(output_file, dpi=150, bbox_inches="tight")
            print(f"✓ Saved to: {output_file}")
        except Exception as e:
            print(f"WARNING: Could not save figure: {e}")

    # Show plot
    print("✓ Displaying plot...")
    plt.show()

    print("=" * 80 + "\n")


if __name__ == "__main__":
    import sys

    # Get output directory from command line argument if provided
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "output"

    plot_rotation_evolution(output_dir=output_dir, save=True)
