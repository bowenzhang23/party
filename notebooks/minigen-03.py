import numpy as np
import matplotlib.pyplot as plt

from party.histogram import *
from party.vector import *


def boost_z_direction(beta, p):
    """Simply boost to Z-axis
    """
    gamma = np.sqrt(1.0 / (1.0 - np.square(beta)))
    p_boost = p.copy()
    p_boost[:, 0] = gamma * p[:, 0] - gamma * beta * p[:, 3]
    p_boost[:, 3] = gamma * p[:, 3] - gamma * beta * p[:, 0]
    return p_boost


def main():
    """An example of examing generated pp->Zy->mumu events
    """
    cost, beta, ecm_hat = (
        np.loadtxt("minigen/outputs/pp_zy_mumu_cost.txt"),
        np.loadtxt("minigen/outputs/pp_zy_mumu_beta.txt"),
        np.loadtxt("minigen/outputs/pp_zy_mumu_ecm_hat.txt"),
    )

    theta = np.arccos(cost)
    phi = np.random.uniform(0.0, 1.0, len(theta)) * 2 * np.pi

    # mu-
    p_cm_mu = np.ones((len(theta), 4))
    p_cm_mu[:, 1] = np.sin(theta) * np.cos(phi)
    p_cm_mu[:, 2] = np.sin(theta) * np.sin(phi)
    p_cm_mu[:, 3] = cost
    # mu+
    p_cm_mu_ = -p_cm_mu
    p_cm_mu_[:, 0] = 1
    p_cm_mu *= 0.5 * ecm_hat.reshape(-1, 1)
    p_cm_mu_ *= 0.5 * ecm_hat.reshape(-1, 1)
    # boost into lab frame
    p_mu = boost_z_direction(beta, p_cm_mu)
    p_mu_ = boost_z_direction(beta, p_cm_mu_)

    p = p_mu + p_mu_
    m = np.sqrt(p[:, 0] ** 2 - (p[:, 1] ** 2 + p[:, 2] ** 2 + p[:, 3] ** 2))
    pt = np.sqrt(p_mu[:, 1] ** 2 + p_mu[:, 2] ** 2)
    np.allclose(m, ecm_hat)

    h_m = make_hist1d(m, 100, 60, 200)
    h_pt = make_hist1d(pt, 100, 0, 100)
    h_cost = make_hist1d(cost, 100, -1, 1)

    plt.ion()

    plt.figure()
    h_m.plot(label="", color="black", density=False)
    plt.xlabel(r"$M_{\mu\mu}$ [GeV]")
    plt.ylabel("Events")
    plt.yscale("log")
    plt.ylim(1, 1e4)
    plt.show()

    plt.figure()
    h_pt.plot(label="", color="black", density=False)
    plt.xlabel(r"$\mu^{-}$ $p_{T}$ [GeV]")
    plt.ylabel("Events")
    plt.show()

    plt.figure()
    h_cost.plot(label="", color="black", density=False)
    plt.xlabel(r"$\mu^{-}$ $\cos(\theta)$")
    plt.ylabel("Events")
    plt.ylim(0, 200)
    plt.show()

    plt.waitforbuttonpress()


if __name__ == "__main__":
    main()
