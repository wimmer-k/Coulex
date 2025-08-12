#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
import sys

def read_Pm():
    if not sys.stdin.isatty():
        data = sys.stdin.read().strip()
        P_m = np.array([float(x) for x in data.split()])
    else:
        print("Enter five P_m values separated by spaces (for m = -2 to +2):")
        P_m = np.array([float(x) for x in input().strip().split()])
    if len(P_m) != 5:
        raise ValueError("Expected 5 P_m values for m = -2, -1, 0, +1, +2.")
    if not np.isclose(np.sum(P_m), 1.0, atol=1e-5):
        raise ValueError("P_m values must sum to 1.")
    return P_m

def compute_alignment(P_m):
    rho2_0 = (1 / np.sqrt(14)) * (2 * (P_m[0] + P_m[4]) - (P_m[1] + P_m[3]))
    rho4_0 = (1 / (2 * np.sqrt(70))) * (6 * (P_m[0] + P_m[4]) - 4 * (P_m[1] + P_m[3]) + 2 * P_m[2])
    return rho2_0, rho4_0

def angular_coefficients(rho2_0, rho4_0, delta=np.inf):
    if np.isinf(delta):  # pure E2
        a2 = np.sqrt(5 / 14) * rho2_0
        a4 = np.sqrt(9 / 70) * rho4_0
    else:
        delta2 = delta ** 2
        norm = 1 / (1 + delta2)
        # Approximate: E2 + M1 contribution with no interference terms for simplicity
        a2 = norm * (np.sqrt(5 / 14) * rho2_0)
        a4 = norm * (np.sqrt(9 / 70) * rho4_0)
    return a2, a4

def plot_W_theta(P_m, a2, a4):
    J = 2
    m_values = np.arange(-J, J + 1)
    theta = np.linspace(0, np.pi, 180)
    cos_theta = np.cos(theta)
    P2 = legendre(2)(cos_theta)
    P4 = legendre(4)(cos_theta)

    W_theta = 1 + a2 * P2 + a4 * P4
    W_theta /= np.max(W_theta)

    fig, ax = plt.subplots()
    ax.plot(np.degrees(theta), W_theta, label="W(theta)")
    ax.set_xlabel("theta (degrees)")
    ax.set_ylabel("W(theta)")
    ax.set_title("Gamma Angular Distribution")

    inset_ax = fig.add_axes([0.6, 0.6, 0.3, 0.3])
    inset_ax.bar(m_values, P_m, color="gray", edgecolor="black")
    inset_ax.set_title("P(m)")
    inset_ax.set_xticks(m_values)
    #inset_ax.set_ylim(0, 0.35)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate and plot angular distribution W(theta) from P_m.")
    parser.add_argument("--Pm", type=float, nargs=5, metavar=('Pm_-2', 'Pm_-1', 'Pm_0', 'Pm_+1', 'Pm_+2'),
                        help="Five P_m values for m = -2 to +2, must sum to 1.")
    parser.add_argument("--delta", type=float, default=float("inf"),
                        help="Multipole mixing ratio delta = E2/M1 (default: pure E2)")
    args = parser.parse_args()

    if args.Pm is None:
        print("Error: You must provide the -Pm argument with 5 values.")
        sys.exit(1)

    P_m = np.array(args.Pm)
    if not np.isclose(np.sum(P_m), 1.0, atol=1e-4):
        raise ValueError("P_m values must sum to 1. Sum(P_m) = ", np.sum(P_m))

    rho2_0, rho4_0 = compute_alignment(P_m)
    a2, a4 = angular_coefficients(rho2_0, rho4_0, delta=args.delta)
    print(f"rho2_0 = {rho2_0:.5f}")
    print(f"rho4_0 = {rho4_0:.5f}")
    print(f"a2 = {a2:.5f}")
    print(f"a4 = {a4:.5f}")
    plot_W_theta(P_m, a2, a4)
