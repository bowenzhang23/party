import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma
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
            self._bin_edge = np.array(bins)
        elif isinstance(bins, int):
            self._bin_edge = np.linspace(start, end, bins + 1)
        else:
            raise RuntimeError(f"bins of type {type(bins)} is not supported!")
        self._n_bins = len(self._bin_edge) + 1
        self._xmin = self._bin_edge[0]
        self._xmax = self._bin_edge[-1]
        self._max = 0.0
        self._min = 0.0
        self._bin_centre = 0.5 * (self._bin_edge[:-1] + self._bin_edge[1:])
        self._bin_sumw = np.zeros((self._n_bins,))
        self._bin_sumw2 = np.zeros((self._n_bins,))
        self._bin_error_up = np.zeros((self._n_bins,))
        self._bin_error_dn = np.zeros((self._n_bins,))
        self._bin_error_type = BinErrorType.Normal

    def _use_normal_error(self) -> bool:
        """Check if normal error should be forced to use"""
        return self._bin_error_type == BinErrorType.Normal or not np.allclose(
            self._bin_sumw, self._bin_sumw2
        )

    def _calculate_poisson_error(self, sumw: float) -> Tuple:
        """Calculate poisson error given the sum of weights

        Args:
            sumw (float): sum of weights

        Returns:
            Tuple: error (up) and error (down)
        """
        if self._use_normal_error():
            return None
        if self._bin_error_type == BinErrorType.Poisson:
            interval = 0.682689492
        elif self._bin_error_type == BinErrorType.Poisson2:
            interval = 0.95
        n = int(sumw)
        if n > 0:
            error_dn = sumw - gamma.ppf(0.5 - interval / 2, a=n)
            error_up = gamma.ppf(0.5 + interval / 2, a=n + 1) - sumw
        else:
            # here use sumw! should warning!
            error_up = error_dn = np.sqrt(np.abs(sumw))
        return error_up, error_dn

    def _calculate_error(self, sumw: float, sumw2: float) -> Tuple:
        """Calculate error given sum of weights and sum of weight squares
        final error type is determined automatically

        Args:
            sumw (float): sum of weights
            sumw2 (float): sum of weight squares

        Returns:
            Tuple: error (up) and error (down)
        """
        if self._use_normal_error():
            return (np.sqrt(np.abs(sumw2)),) * 2
        return self._calculate_poisson_error(sumw)

    def _reset_bin_content(self, idx: int, sum: float, sum2: float) -> None:
        """Reset the content of a bin, will call updating of the errors

        Args:
            idx (int): index
            sum (float): sum of weights to accumulate
            sum2 (float): sum of weight squares to accumulate
        """
        self._bin_sumw[idx] = sum
        self._bin_sumw2[idx] = sum2
        c = self._bin_sumw[idx]
        self._max = max(self._max, c)
        self._min = min(self._min, c)
        self._update_bin_error(idx)

    def _update_bin_content(self, idx: int, sum: float, sum2: float) -> None:
        """Update the content of a bin, will call updating of the errors

        Args:
            idx (int): index
            sum (float): sum of weights to accumulate
            sum2 (float): sum of weight squares to accumulate
        """
        self._reset_bin_content(
            idx, self._bin_sumw[idx] + sum, self._bin_sumw2[idx] + sum2
        )

    def _update_bin_error(self, idx: int) -> None:
        """Update the errors of a bin

        Args:
            idx (int): index
        """
        self._bin_error_up[idx], self._bin_error_dn[idx] = self._calculate_error(
            self._bin_sumw[idx], self._bin_sumw2[idx]
        )

    def fill(self, value: float, weight=None) -> None:
        """Fill a value with weight to the histogram

        Args:
            value (float): the value
            weight (optional): the weight. Defaults to None.
        """
        weight = 1.0 if weight is None else weight
        if value < self._xmin:
            idx = 0
        elif value > self._xmax:
            idx = self._n_bins - 1
        else:
            idx = np.searchsorted(self._bin_edge, value, side="right")
        self._update_bin_content(idx, weight, weight * weight)

    def fill_array(self, values: np.ndarray, weights=None) -> None:
        """Fill an array of value with an array of weights of the same size

        Args:
            values (np.ndarray): the values
            weights (optional): the weights. Defaults to None.
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
        """Set error type manually"""
        self._bin_error_type = error_type

    # region Getters
    def _debug_data(self) -> Any:
        """Return detailed data hold by the object for debug purpose"""
        data_dict = {
            "edge": self._bin_edge,
            "centre": self._bin_centre,
            "sumw": self._bin_sumw,
            "sumw2": self._bin_sumw2,
            "error_up": self._bin_error_up,
            "error_dn": self._bin_error_dn,
        }
        return data_dict

    def data(self, spec: str = "sumw") -> np.ndarray:
        """Return data specifying a name

        Args:
            spec (str, optional): name. Defaults to "sumw".

        Returns:
            np.ndarray: the underlying data
        """
        return self._debug_data()[spec]

    def get_n_bins(self) -> Any:
        return self._n_bins

    def get_max(self) -> Any:
        return self._max

    def get_min(self) -> Any:
        return self._min

    def get_bin_sumw(self, i) -> Any:
        return self._bin_sumw[i]

    def get_bin_sumw2(self, i) -> Any:
        return self._bin_sumw2[i]

    def get_bin_error_up(self, i) -> Any:
        return self._bin_error_up[i]

    def get_bin_error_dn(self, i) -> Any:
        return self._bin_error_dn[i]

    # endregion

    # region arithmetic

    def copy(self) -> Any:
        """
        Returns:
            Any: the copied Hist1D object
        """
        _copy = Hist1D(self._bin_edge)
        _copy._bin_sumw = self._bin_sumw.copy()
        _copy._bin_sumw2 = self._bin_sumw2.copy()
        _copy._bin_error_up = self._bin_error_up.copy()
        _copy._bin_error_dn = self._bin_error_dn.copy()
        _copy._bin_error_type = self._bin_error_type

        return _copy

    def add(self, other: Any, scale: float = 1.0) -> Any:
        """Add scale x other (histogram)

        Args:
            other (Any): Another Hist1D object
            scale (float, optional): scaling of other. Defaults to 1.0.

        Returns:
            Any: new Hist1D that equals this + scale x other,
            errors are recalculated
        """
        copy = self.copy()
        n = copy.get_n_bins()
        assert n == other.get_n_bins(), (
            f"Must have same number of bins, "
            f"however get {n} and {other.get_n_bins()}"
        )
        for i in range(n):
            copy._update_bin_content(
                i, scale * other.get_bin_sumw(i), scale * scale * other.get_bin_sumw2(i)
            )
        return copy

    def __add__(self, other: Any) -> Any:
        return self.add(other, 1.0)

    def __sub__(self, other: Any) -> Any:
        return self.add(other, -1.0)

    def scale(self, k: float = 1.0) -> Any:
        """Scalar scaling of this Hist1D

        Args:
            k (float, optional): scaling factor. Defaults to 1.0.

        Returns:
            Any: new Hist1D that equals k x this
            errors are recalculated
        """
        copy = self.copy()
        n = copy.get_n_bins()
        for i in range(n):
            copy._reset_bin_content(
                i, k * copy.get_bin_sumw(i), k * k * copy.get_bin_sumw2(i)
            )
        return copy

    def multiply(self, other: Any) -> Any:
        """Multiplication of two Hist1D

        Args:
            other (Any): Another Hist1D object

        Returns:
            Any: new Hist1D that equals this x other,
            errors are recalculated
        """
        copy = self.copy()
        n = copy.get_n_bins()
        assert n == other.get_n_bins(), (
            f"Must have same number of bins, "
            f"however get {n} and {other.get_n_bins()}"
        )
        for i in range(n):
            sumw_this, sumw2_this = copy.get_bin_sumw(i), copy.get_bin_sumw2(i)
            sumw_that, sumw2_that = other.get_bin_sumw(i), other.get_bin_sumw2(i)
            prod = sumw_this * sumw_that
            prod_sumw2_estimate = 0
            if sumw_this != 0:
                prod_sumw2_estimate += sumw2_this / (sumw_this**2)
            if sumw_that != 0:
                prod_sumw2_estimate += sumw2_that / (sumw_that**2)
            prod_sumw2_estimate *= prod**2
            copy._reset_bin_content(i, prod, prod_sumw2_estimate)
        return copy

    def divide(self, other: Any) -> Any:
        """Division of two Hist1D

        Args:
            other (Any): Another Hist1D object

        Returns:
            Any: new Hist1D that equals this / other,
            errors are recalculated, zero division are restricted to zero
        """
        copy = self.copy()
        n = copy.get_n_bins()
        assert n == other.get_n_bins(), (
            f"Must have same number of bins, "
            f"however get {n} and {other.get_n_bins()}"
        )
        for i in range(n):
            sumw_this, sumw2_this = copy.get_bin_sumw(i), copy.get_bin_sumw2(i)
            sumw_that, sumw2_that = other.get_bin_sumw(i), other.get_bin_sumw2(i)
            div = 0
            div_sumw2_estimate = 0
            if sumw_this != 0 and sumw_that != 0:
                div = sumw_this / sumw_that
                div_sumw2_estimate += sumw2_this / (sumw_this**2)
                div_sumw2_estimate += sumw2_that / (sumw_that**2)
                div_sumw2_estimate *= div**2
            copy._reset_bin_content(i, div, div_sumw2_estimate)
        return copy

    def __mul__(self, other: Any) -> Any:
        if isinstance(other, float):
            return self.scale(other)
        elif isinstance(other, Hist1D):
            return self.multiply(other)
        else:
            raise ArithmeticError(f"Can not multiply Hist1D with {type(other)}")

    def __truediv__(self, other: Any) -> Any:
        if isinstance(other, float):
            return self.scale(1.0 / other)
        elif isinstance(other, Hist1D):
            return self.divide(other)
        else:
            raise ArithmeticError(f"Can not divide Hist1D with {type(other)}")

    def __rmul__(self, other: float) -> Any:
        return self.scale(other)

    def integral(self) -> Tuple:
        """Get integral and up/down error of the histogram

        Returns:
            Tuple: integral, up error, down error
        """
        integral = np.sum(self._bin_sumw)
        error_up, error_dn = self._calculate_error(integral, np.sum(self._bin_sumw2))
        return integral, error_up, error_dn

    # endregion

    def plot(
        self, label: str, color: str, density=False, underflow=False, overflow=False
    ) -> None:
        """plotting function

        TODO: plotting with underflow and overflow are not implemented yet

        Args:
            label (str): label displayed in legend.
            color (str): the color.
            density (bool, optional): plot density or not. Defaults to False.
            underflow (bool, optional): include underflow in the first bin. Defaults to False.
            overflow (bool, optional): include overflow in the last bin. Defaults to False.
        """
        bin_content = self._bin_sumw[1:-1].copy()
        bin_content_dn = bin_content - self._bin_error_dn[1:-1].copy()
        bin_content_up = bin_content + self._bin_error_up[1:-1].copy()
        if density:
            integral = np.sum(bin_content)
            bin_content /= integral
            bin_content_dn /= integral
            bin_content_up /= integral
            # print(bin_content, bin_content_dn, bin_content_up)
        hline_data = bin_content, self._bin_edge[:-1], self._bin_edge[1:]
        vline_data = self._bin_centre, bin_content_dn, bin_content_up
        plt.hlines(*hline_data, colors=color, label=label)
        plt.vlines(*vline_data, colors=color)


def make_Hist1D(
    arr: np.ndarray,
    nbin: int,
    start: float,
    end: float,
    error_type: BinErrorType = BinErrorType.Normal,
) -> Hist1D:
    """Helper function to create a hist1d object with unity weight

    Args:
        arr (np.ndarray): values to fill
        nbin (int): number of bins
        start (float): lower bin edge
        end (float): upper bin edge
        error_type (BinErrorType, optional): error type. Defaults to BinErrorType.Normal.

    Returns:
        Hist1D: Hist1D object
    """
    h = Hist1D(nbin, start, end)
    h.set_error_type(error_type)
    h.fill_array(arr)
    return h
