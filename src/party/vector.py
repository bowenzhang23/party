import numpy as np
from typing import Any, List, Set, Dict, Callable


# constants
dim_spacetime: int = 4
g_uv: np.ndarray = np.diag([1.0, -1.0, -1.0, -1.0]).astype(np.float64)
varname_set: Set[str] = (
    "E",
    "m",
    "px",
    "py",
    "pz",
    "pt",
    "phi",
    "theta",
    "eta",
    "rapidity",
)


# helper functions
def gamma(beta: Any) -> np.float64:
    """$$\gamma = \frac{1}{\sqrt{1 - {\beta}^{2}}}$$

    Args:
        beta (Any): velocity

    Returns:
        np.float64: gamma, aka Lorentz factor
    """
    return np.reciprocal(np.sqrt(1 - np.square(beta)))


def rapidity(beta: Any) -> np.float64:
    """$$\tanh{y} = \beta$$

    Args:
        beta (Any): velocity

    Returns:
        np.float64: rapidity
    """
    return np.arctanh(beta)


# duality
class Duality(object):
    """Vector dual space type

    Args:
        object (_type_):
    """

    Contravariant: int = 0
    Covariant: int = 1


# four vector class
class FourVector(object):
    def __init__(
        self,
        t: Any,
        x: Any,
        y: Any,
        z: Any,
        label: str = "",
        duality: Duality = Duality.Contravariant,
    ) -> None:
        """Four Vector (Lorentz Vector)

        Args:
            t (Any): time
            x (Any): cartesian x
            y (Any): cartesian y
            z (Any): cartesian z
            label (str, optional): label of the variable for display. Defaults to "".
            duality (Duality, optional): contra- or covariant. Defaults to Duality.Contravariant.
        """
        contra: np.ndarray = np.array([t, x, y, z], dtype=np.float64)
        if duality == Duality.Contravariant:
            self.txyz: np.ndarray = contra
        elif duality == Duality.Covariant:
            self.txyz: np.ndarray = g_uv @ contra
        self.label: str = label
        self.duality: Duality = duality

    @staticmethod
    def zero() -> Any:
        """Zero vector

        Returns:
            Any: zero vector
        """
        return FourVector(*(0 for _ in range(dim_spacetime)))

    def __repr__(self) -> str:
        greek: str = r"u"
        label: str = self.label
        list_view: List[float] = [round(i, 5) for i in self.txyz.tolist()]
        if self.duality == Duality.Contravariant:
            return f"{label}^{greek} = {list_view}"
        elif self.duality == Duality.Covariant:
            return f"{label}_{greek} = {list_view}"
        return ""

    def __getitem__(self, i) -> np.float64:
        assert i >= 0 and i < dim_spacetime, f"must satisfy 0 <= {i} < dim_spacetime"
        return self.txyz[i]

    def _copy(self) -> Any:
        copy = FourVector.zero()
        copy.txyz: np.ndarray = self.txyz
        copy.label: str = self.label
        copy.duality: Duality = self.duality
        return copy

    def _binary_op(self, other: Any, op: str, fun_op: Callable) -> Any:
        res: FourVector = FourVector.zero()
        res.label = (
            f"({self.label}"
            f"{op}"
            f"{other.label if hasattr(other, 'label') else other})"
        )
        res.duality: Duality = self.duality
        other_value: Any = other.txyz if hasattr(other, "txyz") else other
        res.txyz: np.ndarray = fun_op(other_value)
        return res

    def __add__(self, other: Any) -> Any:
        assert isinstance(other, FourVector), "must add with FourVector"
        assert self.duality == other.duality, "must have same duality"
        return self._binary_op(other, "+", self.txyz.__add__)

    def __sub__(self, other: Any) -> Any:
        assert isinstance(other, FourVector), "must subtract with FourVector"
        assert self.duality == other.duality, "must have same duality"
        return self._binary_op(other, "-", self.txyz.__sub__)

    def __mul__(self, other: Any) -> Any:
        if isinstance(other, FourVector):
            # dot product
            res: np.float64 = self.txyz.dot(other.txyz)
        else:
            # scaler multiplication
            assert isinstance(
                other, (int, float, np.int64, np.float64)
            ), "must be int or float"
            res: FourVector = self._binary_op(other, "x", self.txyz.__mul__)
        return res

    def __rmul__(self, other: Any) -> Any:
        # scaler multiplication
        assert isinstance(
            other, (int, float, np.int64, np.float64)
        ), "must be int or float"
        res: FourVector = FourVector.zero()
        return self._binary_op(other, "x", self.txyz.__rmul__)

    def dual_transform(self) -> Any:
        """x_u = g_uv * x^v, x^u = g^uv * x_v

        Returns:
            Any: dual vector
        """
        dual: FourVector = FourVector.zero()
        dual.label: str = self.label
        dual.duality: Duality = self.duality ^ 1
        dual.txyz: np.ndarray = g_uv @ self.txyz
        return dual

    def square(self) -> np.float64:
        """
        returns ||vector||^2
        """
        return self * self

    def rotate(self, angle: Any, axis: int) -> Any:
        """spatial rotation

        Args:
            angle (Any): unit is radian
            axis (int): which axis to rotate on

        Returns:
            Any: transformed four vector
        """
        # passive transformation
        assert (
            axis > 0 and axis < dim_spacetime
        ), "must rotate along x (1), y (2), z (3)"
        a_uv: np.ndarray = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
        cos: np.float64 = np.cos(angle, dtype=np.float64)
        sin: np.float64 = np.sin(angle, dtype=np.float64)
        i, j = (i if i != 0 else 3 for i in [(axis + 1) % 3, (axis + 2) % 3])
        a_uv[i][i] = a_uv[j][j] = cos
        a_uv[i][j] = a_uv[j][i] = sin
        a_uv[j][i] *= -1
        out: FourVector = self._copy()
        out.txyz = a_uv @ out.txyz
        out.label = (
            f"rotate({out.label}, axis={axis}, degree={angle * 180. / np.pi:.3f})"
        )
        return out, a_uv

    def boost(self, velocity: Any, axis: int) -> Any:
        """Lorentz boost

        Args:
            velocity (Any): velocity
            axis (int): which axis to boost on

        Returns:
            np.float64: transformed four vector
        """
        assert velocity > -1 and velocity < 1, "must be lower than the speed of light"
        assert (
            axis > 0 and axis < dim_spacetime
        ), "must rotate along x (1), y (2), z (3)"
        a_uv: np.ndarray = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
        b, y = velocity, gamma(velocity)
        a_uv[0][0] = a_uv[axis][axis] = y
        a_uv[0][3] = a_uv[3][0] = -y * b
        out: FourVector = self._copy()
        out.txyz = a_uv @ out.txyz
        out.label = f"boost({out.label}, axis={axis}, velocity={velocity:.3f})"
        return out, a_uv


# class FourMomentum(FourVector):
#     def __init__(self, vars: Dict[str, Any]):
