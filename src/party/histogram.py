import numpy as np
from scipy.stats import gamma
import matplotlib.pyplot as plt
from typing import Any, Tuple


class BinErrorType:
    """
    | Name     | Description                                                        |
    | -------- | ------------------------------------------------------------------ |
    | Normal   | Errors with Normal (Wald) approximation: errorUp=errorLow= sqrt(N) |
    | Poisson  | Errors from Poisson interval at 68.3% (1 sigma)                    |
    | Poisson2 | Errors from Poisson interval at 95% CL (~ 2 sigma)                 |
    """

    Normal = 0
    Poisson = 1
    Poisson2 = 2


class Hist1D(object):
    def __init__(self, bins: Any, start: float = 0.0, end: float = 1.0) -> None:
        """1-D Histogram

        Args:
            bins (Any): if it is array-like, then it is interpreted as bin edges,
            else if it is a number, then it is interpreted as the number of bins.
            start (float, optional): lower edge of the first bin. Defaults to 0.0.
            end (float, optional): upper edge of the last bin. Defaults to 1.0.

        Raises:
            RuntimeError:
        """
        if isinstance(bins, np.ndarray) or isinstance(bins, list):
            self._bin_edges = np.array(bins)
        elif isinstance(bins, int):
            self._bin_edges = np.linspace(start, end, bins + 1)
        else:
            raise RuntimeError(f"bins of type {type(bins)} is not supported!")
        self._xmin = self._bin_edges[0]
        self._xmax = self._bin_edges[-1]
        self._max = 0.0
        self._min = 0.0
        self._bin_centre = 0.5 * (self._bin_edges[:-1] + self._bin_edges[1:])
        self._bin_content = np.zeros((len(self._bin_edges) + 1,))
        self._bin_error_up = np.zeros((len(self._bin_edges) + 1,))
        self._bin_error_dn = np.zeros((len(self._bin_edges) + 1,))
        self._bin_error_type = BinErrorType.Normal

    def _update_bin_content(self, idx: int, weight: float) -> None:
        self._bin_content[idx] += weight
        c = self._bin_content[idx]
        self._max = max(self._max, c)
        self._min = min(self._min, c)
        self._update_bin_error(idx)

    def _update_bin_error(self, idx: int) -> None:
        c = self._bin_content[idx]
        if self._bin_error_type == BinErrorType.Normal:
            self._bin_error_up[idx] = self._bin_error_dn[idx] = np.sqrt(np.abs(c))
            return
        if self._bin_error_type == BinErrorType.Poisson:
            interval = 0.682689492
        elif self._bin_error_type == BinErrorType.Poisson2:
            interval = 0.95
        n = int(c)
        if n > 0:
            self._bin_error_dn[idx] = c - gamma.ppf(0.5 - interval / 2, a=n)
            self._bin_error_up[idx] = gamma.ppf(0.5 + interval / 2, a=n + 1) - c
        else:
            self._bin_error_up[idx] = self._bin_error_dn[idx] = np.sqrt(np.abs(c))

    def fill(self, value: float, weight: float = None) -> None:
        """Fill a value with weight to the histogram

        Args:
            value (float): the value
            weight (float, optional): the weight. Defaults to None.
        """
        weight = 1.0 if weight is None else weight
        if value < self._xmin:
            idx = 0
        elif value > self._xmax:
            idx = len(self._bin_edges)
        else:
            idx = np.searchsorted(self._bin_edges, value, side="right")
        self._update_bin_content(idx, weight)

    def fill_array(self, values: np.ndarray, weights: np.ndarray = None) -> None:
        """Fill an array of value with an array of weights of the same size

        Args:
            values (np.ndarray): the values
            weights (np.ndarray, optional): the weights. Defaults to None.
        """
        if weights is not None:
            assert values.shape == weights.shape, (
                f"shapes of value and weight array must be identical, "
                f"however, found {values.shape} and {weights.shape}"
            )
        else:
            weights = np.ones_like(values)
        for value, weight in zip(values, weights):
            self.fill(value, weight)

    def set_error_type(self, error_type: BinErrorType) -> None:
        self._bin_error_type = error_type

    # region Getters
    def data(self) -> Tuple:
        return self._bin_content, self._bin_error_up, self._bin_error_dn

    def get_max(self) -> Any:
        return self._max

    def get_min(self) -> Any:
        return self._min

    def get_bin_content(self, i) -> Any:
        return self._bin_content[i]

    def get_bin_error_up(self, i) -> Any:
        return self._bin_error_up[i]

    def get_bin_error_dn(self, i) -> Any:
        return self._bin_error_dn[i]

    def get_bin_centre(self, i) -> Any:
        return self._bin_centre[i - 1]

    # endregion

    def plot(self, label="", color=None, underflow=False, overflow=False) -> None:
        """plotting function

        Args:
            label (str, optional): label displayed in legend. Defaults to "".
            color (_type_, optional): the color. Defaults to None.
            underflow (bool, optional): include underflow in the first bin. Defaults to False.
            overflow (bool, optional): include overflow in the last bin. Defaults to False.
        """
        bin_content = self._bin_content[1:-1]
        bin_content_dn = bin_content - self._bin_error_dn[1:-1]
        bin_content_up = bin_content + self._bin_error_up[1:-1]
        hline_data = bin_content, self._bin_edges[:-1], self._bin_edges[1:]
        vline_data = self._bin_centre, bin_content_dn, bin_content_up
        plt.hlines(*hline_data, colors=color, label=label)
        plt.vlines(*vline_data, colors=color)
