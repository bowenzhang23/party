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


def eta2rapidity(eta: Any, m: Any, p: Any) -> np.float64:
    """y = arctanh(1/sqrt(1+m^2/p^2)tanh(n))

    Args:
        eta (Any): pseudo-rapidity

    Returns:
        np.float64: actual rapidity
    """
    return np.arctanh(1 / np.sqrt(1 + m ^ 2 / p ^ 2) * np.tanh(eta))


def theta2eta(theta: Any) -> np.float64:
    """azimuth angle to pseudo-rapidity, i.e. definition of eta:
    n = -ln(tan(t/2))

    Args:
        theta (Any): azimuth angle

    Returns:
        _type_: pseudo-rapidity
    """
    return -np.log(np.tan(theta * 0.5))


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

    def rotate(self, angle: Any, param: Any) -> Any:
        if isinstance(param, tuple):
            return self._rotate_vector(angle, param)
        else:
            return self._rotate_axis(angle, param)

    def boost(self, velocity: Any, param: Any) -> Any:
        if isinstance(param, tuple):
            return self._boost_vector(velocity, param)
        else:
            return self._boost_axis(velocity, param)

    def boost_to_zmf(self) -> Any:
        E, p = self.txyz[0], self.txyz[1:]
        b: np.ndarray = p / E
        v: np.ndarray = np.linalg.norm(b)
        return self.boost(v, tuple(b.tolist()))

    def _rotate_axis(self, angle: Any, axis: int) -> Any:
        """spatial rotation along a given axis
        x, y, z = 1, 2, 3

        Args:
            angle (Any): unit is radian
            axis (int): which axis to rotate on

        Returns:
            Any: transformed four vector and rotation matrix
        """
        # passive transformation
        assert (
            axis > 0 and axis < dim_spacetime
        ), "must rotate along x (1), y (2), z (3)"
        a_uv: np.ndarray = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
        cos: np.float64 = np.cos(angle, dtype=np.float64)
        sin: np.float64 = np.sin(angle, dtype=np.float64)
        i, j = (i if i != 0 else 3 for i in [(axis + 1) % 3, (axis + 2) % 3])
        a_uv[i, i] = a_uv[j, j] = cos
        a_uv[i, j] = a_uv[j, i] = sin
        a_uv[j, i] *= -1
        out: FourVector = self._copy()
        out.txyz = a_uv @ out.txyz
        out.label = (
            f"rotate({out.label}, axis={axis}, degree={angle * 180. / np.pi:.3f})"
        )
        return out, a_uv

    def _rotate_vector(self, angle: Any, u: tuple) -> Any:
        """spatial rotation along a given direction

        Args:
            angle (Any): counter-clock angle, axis point to observer
            u (tuple): direction

        Returns:
            Any: transformed four vector and rotation matrix

        Reference:
            https://en.wikipedia.org/wiki/Rotation_matrix#cite_ref-5
        """
        assert len(u) == 3, "must provide 3-D vector!"
        un: np.ndarray = np.array(u, dtype=np.float64)
        un = un / np.linalg.norm(un, keepdims=True)
        ux, uy, uz = un[0], un[1], un[2]
        cos: np.float64 = np.cos(angle, dtype=np.float64)
        sin: np.float64 = np.sin(angle, dtype=np.float64)
        a_uv: np.ndarray = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
        a_uv[1:, 1:] = np.array(
            [
                [
                    cos + ux * ux * (1 - cos),
                    ux * uy * (1 - cos) + uz * sin,
                    ux * uz * (1 - cos) - uy * sin,
                ],
                [
                    uy * ux * (1 - cos) - uz * sin,
                    cos + uy * uy * (1 - cos),
                    uy * uz * (1 - cos) + uy * sin,
                ],
                [
                    uz * ux * (1 - cos) + uy * sin,
                    ux * uy * (1 - cos) - ux * sin,
                    cos + uz * uz * (1 - cos),
                ],
            ],
            dtype=np.float64,
        )
        out: FourVector = self._copy()
        out.txyz = a_uv @ out.txyz
        out.label = f"rotate({out.label}, axis={u}, degree={angle * 180. / np.pi:.3f})"
        return out, a_uv

    def _boost_axis(self, velocity: Any, axis: int) -> Any:
        """Lorentz boost along a given axis
        x, y, z = 1, 2, 3

        Args:
            velocity (Any): velocity
            axis (int): which axis to boost on

        Returns:
            Any: transformed four vector and Lorentz transformation matrix
        """
        assert velocity > -1 and velocity < 1, "must be lower than the speed of light"
        assert (
            axis > 0 and axis < dim_spacetime
        ), "must rotate along x (1), y (2), z (3)"
        a_uv: np.ndarray = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
        b, y = velocity, gamma(velocity)
        a_uv[0, 0] = a_uv[axis, axis] = y
        a_uv[0, 3] = a_uv[3, 0] = -y * b
        out: FourVector = self._copy()
        out.txyz = a_uv @ out.txyz
        out.label = f"boost({out.label}, axis={axis}, velocity={velocity:.3f})"
        return out, a_uv

    def _boost_vector(self, velocity: Any, u: tuple):
        """Lorentz boost along a given direction

        Args:
            velocity (Any): velocity
            u (tuple): direction

        Returns:
            Any: transformed four vector and Lorentz transformation matrix

        Reference:
            https://en.wikipedia.org/wiki/Lorentz_transformation
        """
        assert len(u) == 3, "must provide 3-D vector!"
        un: np.ndarray = np.array(u, dtype=np.float64)
        v = un / np.linalg.norm(un, keepdims=True) * velocity
        v2 = velocity**2
        vx, vy, vz = v[0], v[1], v[2]
        a_uv: np.ndarray = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
        y = gamma(velocity)
        a_uv[0][0] = y
        a_uv[0, 1:] = -y * v
        a_uv[1:, 0] = -y * v
        a_uv[1:, 1:] += (y - 1) * np.array(
            [
                [vx * vx / v2, vx * vy / v2, vx * vz / v2],
                [vy * vx / v2, vy * vy / v2, vy * vz / v2],
                [vz * vx / v2, vz * vy / v2, vz * vz / v2],
            ],
            dtype=np.float64,
        )
        out: FourVector = self._copy()
        out.txyz = a_uv @ out.txyz
        out.label = f"boost({out.label}, axis={u}, velocity={velocity:.3f})"
        return out, a_uv


class FourMomentum(FourVector):
    def __init__(
        self,
        vars: Dict[str, Any],
        label: str = "",
        duality: Duality = Duality.Contravariant,
    ):
        """Four momentum can be
        - without rapidity
            - (E, px, py, pz) -> the standard p4
            - (m, px, py, pz) -> p^2 = m^2, i.e. E^2 = p^2 + m^2, take the positive E
            - (E, pt, phi, pz) -> px = pt*cos(phi), py = pt*sin(phi)
            - (m, pt, phi, pz) -> px = pt*cos(phi), py = pt*sin(phi)
        - with rapidity (by convention, always boost along z-direction).
            - (E, px, py, rapidity) -> pz from tanhy = b = pz / E
            - (m, px, py, rapidity) -> pz from sinhy = yb = pz / m
            - (E, pt, phi, rapidity) -> pz from tanhy = b = pz / E
            - (m, pt, phi, rapidity) -> pz from sinhy = yb = pz / m

        Notes:
            rapidity (y) <-> eta (n) <-> theta (t)
            y = arctanh(1/sqrt(1+m^2/p^2)tanh(n))
            n = -ln(tan(t/2))

        Args:
            vars (Dict[str, Any]): _description_
        """
        for var in vars.keys():
            if var not in varname_set:
                raise RuntimeError(
                    f"{var} is not supported, cannot construct FourMomentum!"
                )

        if {"E", "m"}.intersection(vars.keys()) == set():
            raise RuntimeError(f"must have at least one of E and m!")
