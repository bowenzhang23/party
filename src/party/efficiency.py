import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.stats import norm, beta, binom
from party.histogram import Hist1D
from typing import Any, Tuple


class EfficiencyInconsistencyError(RuntimeError): ...


class Efficiency1D(object):
    def __init__(self, h_pass: Hist1D, h_total: Hist1D, interval=0.682689492) -> None:
        """1-D Efficiency

        Args:
            h_pass (Hist1D): Histogram of passed events
            h_total (Hist1D): Histogram of total events
            interval (float, optional): confidence interval. Defaults to 0.683 (1-sigma).

        Raises:
            EfficiencyInconsistencyError: h_pass and h_total are not consistent in nbins
        """
        self._h_pass = h_pass
        self._h_total = h_total
        self._weighted = h_pass.is_weighted() or h_total.is_weighted()
        if not self._check_consistency():
            raise EfficiencyInconsistencyError()
        self._n_bins = self._h_pass.get_n_bins()
        self._bin_edge = self._h_pass.data("edge")
        self._bin_centre = self._h_pass.data("centre")
        self._bin_eff = np.ones((self._n_bins,))
        self._bin_error_up = np.zeros((self._n_bins,))
        self._bin_error_dn = np.zeros((self._n_bins,))
        self._interval = interval
        self._calculate_efficiency()

    def _debug_data(self) -> Any:
        """Return detailed data hold by the object for debug purpose"""
        data_dict = {
            "h_pass": self._h_pass,
            "h_total": self._h_total,
            "edge": self._bin_edge,
            "centre": self._bin_centre,
            "eff": self._bin_eff,
            "error_up": self._bin_error_up,
            "error_dn": self._bin_error_dn,
        }
        return data_dict

    def data(self, spec: str = "eff") -> np.ndarray:
        """Return data specifying a name

        Args:
            spec (str, optional): name. Defaults to "eff".

        Returns:
            np.ndarray: the underlying data
        """
        return self._debug_data()[spec]

    def get_bin_eff(self, i: int) -> Any:
        return self._bin_eff[i]

    def get_bin_error_up(self, i) -> Any:
        return self._bin_error_up[i]

    def get_bin_error_dn(self, i) -> Any:
        return self._bin_error_dn[i]

    def _check_consistency(self) -> bool:
        if self._h_pass.get_n_bins() != self._h_total.get_n_bins():
            return False
        return True
    
    def _calculate_clopper_pearson_interval(self, total: float, passed: float, eff: float) -> Tuple:
        """Calculate binomial error based on number of total and passed events,
        using Clopper Pearson interval

        Args:
            total (float): number of total events
            passed (float): number of passed events
            eff (float): efficiency

        Returns:
            Tuple: error (up) and error (down)
        """
        # round the total and passed events for weights
        nt = round(total)
        np = round(passed)
        alpha = 0.5 - self._interval / 2
        # Clopper-Pearson interval
        error_dn = 0.0 if np == 0 else eff - beta.ppf(alpha, np, nt - np + 1)
        error_up = 0.0 if np == nt else beta.ppf(1 - alpha, np + 1, nt - np) - eff
        return error_up, error_dn

    def _calculate_normal_approximated_interval(self, tw: float, tw2: float, pw: float, pw2: float, eff: float) -> Tuple:
        """Calculate binomial error based on number of total and passed events,
        using normal approximation

        Args:
            tw (float): number of total events
            tw2 (float): squared sum of weight of total events
            pw (float): number of passed events
            pw2 (float): squared sum of weight of passed events
            eff (float): efficiency

        Returns:
            Tuple: error (up) and error (down)
        """
        if tw == 0:
            return 0., 0.
        variance = (pw2 * (1. - 2*eff) + tw2 * eff * eff) / (tw * tw)
        sigma = np.sqrt(variance)
        prob = 0.5 + self._interval / 2
        delta = norm.ppf(prob, 0, sigma)
        # normal interval ..
        error_dn = 0.0 if pw == 0 else eff if eff < delta else delta
        error_up = 0.0 if pw == tw else 1 - eff if 1 - eff < delta else delta 
        return error_up, error_dn

    def _calculate_interval(self, i: int) -> Tuple:
        """Calculate binomial error based on number of total and passed events,
        if weighted, then force to use normal approximation,
        else if not, then use Clopper-Pearson

        Returns:
            Tuple: error (up) and error (down)
        """
        tw = self._h_total.get_bin_sumw(i)
        tw2 = self._h_total.get_bin_sumw2(i)
        pw = self._h_pass.get_bin_sumw(i)
        pw2 = self._h_pass.get_bin_sumw2(i)
        eff = self._bin_eff[i]
        
        if self._weighted:
            return self._calculate_normal_approximated_interval(tw, tw2, pw, pw2, eff)
        else:
            return self._calculate_clopper_pearson_interval(tw, pw, eff)

    def _calculate_efficiency(self) -> None:
        """Calculate the efficiency as well as its up and down errors"""
        h_ratio = self._h_pass.divide(self._h_total)
        self._bin_eff = h_ratio.data("sumw")
        for i in range(self._n_bins):
            self._bin_error_up[i], self._bin_error_dn[i] = self._calculate_interval(i)

    def plot(self, label: str, color: str, underflow=False, overflow=False) -> Any:
        """plotting function

        TODO: plotting with underflow and overflow are not implemented yet

        Args:
            label (str): label displayed in legend.
            color (str): the color.

        Returns:
            Any: legend handle
        """
        bin_content = self._bin_eff[1:-1].copy()
        bin_content_dn = bin_content - self._bin_error_dn[1:-1].copy()
        bin_content_up = bin_content + self._bin_error_up[1:-1].copy()
        hline_data = bin_content, self._bin_edge[:-1], self._bin_edge[1:]
        vline_data = self._bin_centre, bin_content_dn, bin_content_up
        plt.hlines(*hline_data, colors=color, label=label)
        plt.vlines(*vline_data, colors=color)
        return mlines.Line2D([], [], color=color, marker='|', linestyle='-', markersize=10, label=label)
