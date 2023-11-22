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
        self._n_bins = len(self._bin_edges)
        self._xmin = self._bin_edges[0]
        self._xmax = self._bin_edges[-1]
        self._max = 0.0
        self._min = 0.0
        self._bin_centre = 0.5 * (self._bin_edges[:-1] + self._bin_edges[1:])
        self._bin_content = np.zeros((len(self._bin_edges) + 1,))
        self._bin_error_up = np.zeros((len(self._bin_edges) + 1,))
        self._bin_error_dn = np.zeros((len(self._bin_edges) + 1,))
        self._bin_error_type = BinErrorType.Normal

    @staticmethod
    def _get_bin_content_s(idx: int, weight: float, content: np.ndarray) -> Any:
        content[idx] += weight
        return content[idx]

    @staticmethod
    def _get_error_s(c: float, error_type: BinErrorType) -> Tuple:
        error_up = error_dn = 0
        if error_type == BinErrorType.Normal:
            error_up = error_dn = np.sqrt(np.abs(c))
            return error_up, error_dn
        if error_type == BinErrorType.Poisson:
            interval = 0.682689492
        elif error_type == BinErrorType.Poisson2:
            interval = 0.95
        n = int(c)
        if n > 0:
            error_dn = c - gamma.ppf(0.5 - interval / 2, a=n)
            error_up = gamma.ppf(0.5 + interval / 2, a=n + 1) - c
        else:
            error_up = error_dn = np.sqrt(np.abs(c))
        return error_up, error_dn

    @staticmethod
    def _get_bin_error_s(
        idx: int,
        c: float,
        error_up: np.ndarray,
        error_dn: np.ndarray,
        error_type: BinErrorType,
    ) -> Tuple:
        error_up[idx], error_dn[idx] = Hist1D._get_error_s(c, error_type)
        return error_up[idx], error_dn[idx]

    @staticmethod
    def _get_integral_s(content: np.ndarray) -> Any:
        return np.sum(content)

    @staticmethod
    def _integral_s(h: Any, error_type: BinErrorType) -> Tuple:
        integral = Hist1D._get_integral_s(h.data("content"))
        error_up, error_dn = Hist1D._get_error_s(integral, error_type)
        return integral, error_up, error_dn

    def _update_bin_content(self, idx: int, weight: float) -> None:
        c = Hist1D._get_bin_content_s(idx, weight, self._bin_content)
        self._max = max(self._max, c)
        self._min = min(self._min, c)
        self._update_bin_error(idx)

    def _update_bin_error(self, idx: int) -> None:
        Hist1D._get_bin_error_s(
            idx,
            self._bin_content[idx],
            self._bin_error_up,
            self._bin_error_dn,
            self._bin_error_type,
        )

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

    def data(self, spec: str = "content") -> np.ndarray:
        data_dict = {
            "content": self._bin_content,
            "error_up": self._bin_error_up,
            "error_dn": self._bin_error_dn,
        }
        return data_dict[spec]

    def get_n_bins(self) -> Any:
        return self._n_bins

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

    # region arithmetic

    def copy(self) -> Any:
        """
        Returns:
            Any: the copied Hist1D object
        """
        _copy = Hist1D(self._bin_edges)
        _copy._bin_content = self._bin_content.copy()
        _copy._bin_error_up = self._bin_error_up.copy()
        _copy._bin_error_dn = self._bin_error_dn.copy()

        return _copy

    def add(self, other: Any, scale: float = 1.0) -> Any:
        copy = self.copy()
        n = copy.get_n_bins()
        assert n == other.get_n_bins(), (
            f"Must have same number of bins, "
            f"however get {n} and {other.get_n_bins()}"
        )
        for i in range(n):
            copy._update_bin_content(i, scale * other.get_bin_content(i))
        return copy

    def __add__(self, other: Any) -> Any:
        return self.add(other, 1.0)

    def __sub__(self, other: Any) -> Any:
        return self.add(other, -1.0)

    def scale(self, k: float = 1.0) -> Any:
        copy = self.copy()
        n = copy.get_n_bins()
        for i in range(n):
            copy._update_bin_content(i, k * copy.get_bin_content(i))
        return copy

    def __mul__(self, k: float) -> Any:
        return self.scale(k)

    def __rmul__(self, k: float) -> Any:
        return self.scale(k)

    # endregion

    def plot(
        self, label="", color=None, density=False, underflow=False, overflow=False
    ) -> None:
        """plotting function

        Args:
            label (str, optional): label displayed in legend. Defaults to "".
            color (_type_, optional): the color. Defaults to None.
            underflow (bool, optional): include underflow in the first bin. Defaults to False.
            overflow (bool, optional): include overflow in the last bin. Defaults to False.
        """
        bin_content = self._bin_content[1:-1].copy()
        bin_content_dn = bin_content - self._bin_error_dn[1:-1].copy()
        bin_content_up = bin_content + self._bin_error_up[1:-1].copy()
        if density:
            integral = Hist1D._get_integral_s(bin_content)
            bin_content /= integral
            bin_content_dn /= integral
            bin_content_up /= integral
            # print(bin_content, bin_content_dn, bin_content_up)
        hline_data = bin_content, self._bin_edges[:-1], self._bin_edges[1:]
        vline_data = self._bin_centre, bin_content_dn, bin_content_up
        plt.hlines(*hline_data, colors=color, label=label)
        plt.vlines(*vline_data, colors=color)
