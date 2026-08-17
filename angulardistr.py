#!/usr/bin/env python3

import sys
sys.path.insert(0, "/home/wimmer/programs/nilsson")

import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import eval_legendre
from wigner import CG, wigner3J, wigner6J

def m_values(J):
    return np.arange(-J, J + 1, 1.0)

def rhok(J, k, Pm):
    # Yamazaki Eq1
    s = 0.0
    for m, p in zip(m_values(J), Pm):
        s += (-1)**int(round(J + m)) * wigner3J(J, J, k, -m, m, 0) * p        
    return np.sqrt((2*k + 1)*(2*J + 1)) * s

def Ak(Ji, Jf, L1, L2, delta, k):
    # Yamazaki Eq 3, excl rho!!
    if np.isinf(delta):
        return Fk(Ji, Jf, L2, L2, k)
    return 1/ (1 + delta**2) * (Fk(Ji, Jf, L1, L1, k) + 2*delta*Fk(Ji, Jf, L1, L2, k) + delta**2*Fk(Ji, Jf, L2, L2, k))

def Fk(Ji, Jf, L1, L2, k):
    # Yamazaki Eq 4
    return (-1)**int(round(Jf + Ji + 1)) * np.sqrt((2*L1 + 1)*(2*L2 + 1)*(2*Ji + 1)) * CG(L1, 1, L2, -1, k, 0) * wigner6J(L1, L2, k, Ji, Ji, Jf)

def coefficients(Ji, Jf, Pm, L1, L2, delta):
    kmax = int(min(2*Ji, 2*max(L1, L2)))
    ks = range(0, kmax + 1, 2)

    rho = {k: rhok(Ji, k, Pm) for k in ks}
    A = {k: Ak(Ji, Jf, L1, L2, delta, k) for k in ks}

    norm = A[0] * rho[0]
    a = {k: A[k]*rho[k]/norm for k in ks if k > 0}
    return A, rho, a

def Wtheta(theta_deg, a):
    theta = np.deg2rad(theta_deg)
    W = np.ones_like(theta)
    for k, ak in a.items():
        W += ak * eval_legendre(k, np.cos(theta))
    return W / (4*np.pi)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--Ji", type=float, required=True)
    parser.add_argument("--Jf", type=float, required=True)
    parser.add_argument("-Pm", "--Pm", type=float, nargs="+", required=True)
    parser.add_argument("--L1", type=int, default=None)
    parser.add_argument("--L2", type=int, default=None)
    parser.add_argument("--delta", type=float, default=0.0)
    args = parser.parse_args()

    Ji, Jf = args.Ji, args.Jf
    Pm = np.array(args.Pm, dtype=float)

    if len(Pm) != int(2*Ji + 1):
        raise ValueError(f"Ji={Ji} requires {int(2*Ji+1)} Pm values.")
    if not np.isclose(Pm.sum(), 1.0):
        raise ValueError(f"Pm must sum to 1, currently {Pm.sum()}")

    L1 = args.L1 if args.L1 is not None else max(1, int(abs(Ji - Jf)))
    L2 = args.L2 if args.L2 is not None else L1 + 1

    A, rho, a = coefficients(Ji, Jf, Pm, L1, L2, args.delta)

    print(f"Ji={Ji:g} -> Jf={Jf:g}, L1={L1}, L2={L2}, delta={args.delta:g}")
    print(f"Pm = {Pm.tolist()}")
    for k in sorted(rho):
        print(f"rho{k} = {rho[k]: .8f}   A{k} = {A[k]: .8f}", end="")
        if k in a:
            print(f"   a{k} = {a[k]: .8f}")
        else:
            print()

    theta = np.linspace(0, 180, 181)
    W = Wtheta(theta, a)

    fig, ax = plt.subplots()
    ax.plot(theta, W, lw=2)
    ax.set_xlabel(r"$\theta$ (degrees)")
    ax.set_ylabel(r"$W(\theta)$")
    ax.set_xlim(0, 180)
    ax.set_title(rf"$J_i={Ji:g}\rightarrow J_f={Jf:g}$")

    inset = ax.inset_axes([0.6, 0.60, 0.30, 0.30])
    inset.bar(m_values(Ji), Pm, color="gray", edgecolor="black")
    inset.set_xlabel("m")
    inset.set_ylabel("P(m)")
    inset.set_xticks(m_values(Ji))

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
