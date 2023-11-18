import numpy as np
from typing import Any, Dict, Callable


__author__ = "Bowen Zhang"
__copyright__ = "Copyright (C) 2023 Bowen Zhang"
__license__ = "MIT License"
__version__ = "0.1"


# constants
dim_spacetime = 4
g_uv = np.diag([1.0, -1.0, -1.0, -1.0]).astype(np.float64)
varname_set = (
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
    """y = 1/sqrt(1 - beta^2)

    Args:
        beta (Any): velocity

    Returns:
        np.float64: gamma, aka Lorentz factor
    """
    return np.reciprocal(np.sqrt(1 - np.square(beta)))


def rapidity(beta: Any) -> np.float64:
    """tanhy = beta

    Args:
        beta (Any): velocity

    Returns:
        np.float64: rapidity
    """
    return np.arctanh(beta)


def pz_from_pt_eta(pt: Any, eta: Any) -> np.float64:
    """See README.md: Calculating pz with $eta$ ($theta$)

    Args:
        pt (Any): transverse momentum
        eta (Any): pseudorapidity from zenith angle theta (w.r.t. z-axis)

    Returns:
        np.float64: momentum in z-direction
    """
    return pt * np.sinh(eta)


def theta2eta(theta: Any) -> np.float64:
    """azimuth angle to pseudo-rapidity, i.e. definition of eta:
    n = -ln(tan(t/2))

    Args:
        theta (Any): azimuth angle

    Returns:
        _type_: pseudo-rapidity
    """
    return -np.log(np.tan(theta * 0.5))


def cos_theta(a: Any, b: Any) -> np.float64:
    """Cosine angle between the two spatial vectors

    Returns:
        np.float64: cos_theta
    """
    assert a.shape == (3,) and b.shape == (3,)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if np.allclose(norm_a, 0, atol=1e-8) or np.allclose(norm_b, 0, atol=1e-8):
        rnd = np.random.rand(1) * 2 - 1
        return rnd.item()
    u_a = a / norm_a
    u_b = b / norm_b
    return np.sum(u_a * u_b)


# duality
class Duality(object):
    """Vector dual space type

    Args:
        object (_type_):
    """

    Contravariant: int = 0
    Covariant: int = 1


# four vector class
class FourVectorBase(object):
    def __init__(
        self,
        vec: np.ndarray,
        label: str = "",
        duality: Duality = Duality.Contravariant,
    ) -> None:
        """Four Vector (Lorentz Vector)

        Args:
            vec (np.ndarray): four vector
            label (str, optional): label of the variable for display. Defaults to "".
            duality (Duality, optional): contra- or covariant. Defaults to Duality.Contravariant.
        """
        contra = vec.copy().astype(np.float64).reshape(-1)
        assert len(contra) == dim_spacetime, f"must satisfy len == {dim_spacetime}"
        if duality == Duality.Contravariant:
            self.txyz = contra
        elif duality == Duality.Covariant:
            self.txyz = g_uv @ contra
        self.label = label
        self.duality = duality

    @staticmethod
    def zero() -> Any:
        """Zero vector

        Returns:
            Any: zero vector
        """
        return FourVectorBase(np.zeros((dim_spacetime,)))

    def __repr__(self) -> str:
        greek = r"u"
        label = self.label
        list_view = [round(i, 5) for i in self.txyz.tolist()]
        if self.duality == Duality.Contravariant:
            return f"{label}^{greek} = {list_view}"
        elif self.duality == Duality.Covariant:
            return f"{label}_{greek} = {list_view}"
        return ""

    def __getitem__(self, i) -> np.float64:
        assert i >= 0 and i < dim_spacetime, f"must satisfy 0 <= {i} < dim_spacetime"
        return self.txyz[i]

    def _copy(self) -> Any:
        copy = FourVectorBase.zero()
        copy.txyz = self.txyz
        copy.label = self.label
        copy.duality = self.duality
        return copy

    def _binary_op(self, other: Any, op: str, fun_op: Callable) -> Any:
        res: FourVectorBase = FourVectorBase.zero()
        res.label = (
            f"({self.label}"
            f"{op}"
            f"{other.label if hasattr(other, 'label') else other})"
        )
        res.duality = self.duality
        other_value = other.txyz if hasattr(other, "txyz") else other
        res.txyz = fun_op(other_value)
        return res

    def __neg__(self) -> Any:
        """negative = -1 * self

        Returns:
            Any: p -> -p
        """
        return -1 * self

    def __add__(self, other: Any) -> Any:
        """add = self + other

        Args:
            other (Any): another four vector

        Returns:
            Any: new four vector
        """
        assert isinstance(other, FourVectorBase), "must add with FourVector"
        assert self.duality == other.duality, "must have same duality"
        return self._binary_op(other, "+", self.txyz.__add__)

    def __sub__(self, other: Any) -> Any:
        """sub = self - other

        Args:
            other (Any): another four vector

        Returns:
            Any: new four vector
        """
        assert isinstance(other, FourVectorBase), "must subtract with FourVector"
        assert self.duality == other.duality, "must have same duality"
        return self._binary_op(other, "-", self.txyz.__sub__)

    def __mul__(self, other: Any) -> Any:
        """sub = self * other, other can be vector (dot product) or scalar

        Args:
            other (Any): four vector or scalar

        Returns:
            Any: scalar or four vector
        """
        if isinstance(other, FourVectorBase):
            # dot product
            if self.duality == other.duality:
                other = other.dual_transform()
            res = self.txyz.dot(other.txyz)
        else:
            # scaler multiplication
            assert isinstance(
                other, (int, float, np.int64, np.float64)
            ), "must be int or float"
            res = self._binary_op(other, "x", self.txyz.__mul__)
        return res

    def __rmul__(self, other: Any) -> Any:
        """sub = other * self

        Args:
            other (Any): scalar

        Returns:
            Any: new four vector
        """
        # scaler multiplication
        assert isinstance(
            other, (int, float, np.int64, np.float64)
        ), "must be int or float"
        return self._binary_op(other, "x", self.txyz.__rmul__)

    def copy(self) -> Any:
        """
        Returns:
            Any: a copy of FourVector
        """
        return self._copy()

    def dual_transform(self) -> Any:
        """x_u = g_uv * x^v, x^u = g^uv * x_v

        Returns:
            Any: dual vector
        """
        dual = FourVectorBase.zero()
        dual.label = self.label
        dual.duality = self.duality ^ 1
        dual.txyz = g_uv @ self.txyz
        return dual

    def opposite(self) -> Any:
        """
        Returns:
            Any: opposite direction in space
        """
        opposite = self._copy()
        opposite.txyz[1:] *= -1
        return opposite

    def vector3(self) -> np.ndarray:
        """
        Returns:
            np.ndarray: cloned 3-D space vector
        """
        return self.txyz[1:].copy()

    def square(self) -> np.float64:
        """
        Returns:
            np.float64: ||vector||^2
        """
        return self * self._copy()

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

    def boost_to_restframe_of(self, other: Any) -> Any:
        """Boost to the rest frame of a specific FourVector

        Returns:
            Any: transformed vector
        """
        E, p = other.txyz[0], other.txyz[1:]
        b = p / E
        v = np.linalg.norm(b)
        return self.boost(v, tuple(b.tolist()))

    def boost_to_restframe(self) -> Any:
        """Boost to the rest frame

        Returns:
            Any: transformed vector
        """
        return self.boost_to_restframe_of(self)

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
        a_uv = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
        cos = np.cos(angle, dtype=np.float64)
        sin = np.sin(angle, dtype=np.float64)
        i, j = (i if i != 0 else 3 for i in [(axis + 1) % 3, (axis + 2) % 3])
        a_uv[i, i] = a_uv[j, j] = cos
        a_uv[i, j] = a_uv[j, i] = sin
        a_uv[j, i] *= -1
        out = self._copy()
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
        un = np.array(u, dtype=np.float64)
        un = un / np.linalg.norm(un, keepdims=True)
        ux, uy, uz = un[0], un[1], un[2]
        cos = np.cos(angle, dtype=np.float64)
        sin = np.sin(angle, dtype=np.float64)
        a_uv = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
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
        out = self._copy()
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
        a_uv = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
        b, y = velocity, gamma(velocity)
        a_uv[0, 0] = a_uv[axis, axis] = y
        a_uv[0, 3] = a_uv[3, 0] = -y * b
        out = self._copy()
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
        un = np.array(u, dtype=np.float64)
        v = un / np.linalg.norm(un, keepdims=True) * velocity
        v2 = velocity**2
        vx, vy, vz = v[0], v[1], v[2]
        a_uv = np.diag([1.0, 1.0, 1.0, 1.0]).astype(np.float64)
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
        out = self._copy()
        out.txyz = a_uv @ out.txyz
        out.label = f"boost({out.label}, axis={u}, velocity={velocity:.3f})"
        return out, a_uv

    def cos_theta(self, other: Any) -> np.float64:
        """Cosine angle between the spatial vector of two FourVectors

        Args:
            other (Any): the other four vector

        Returns:
            np.float64: cos_theta
        """
        return cos_theta(self.vector3(), other.vector3())


class FourVector(FourVectorBase):
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
        super().__init__(np.array([t, x, y, z], dtype=np.float64), label, duality)


class FourMomentum(FourVector):
    def __init__(
        self,
        vars: Dict[str, Any],
        label: str = "",
        duality: Duality = Duality.Contravariant,
    ):
        """Flexible four-momentum construction

        Args:
            vars (Dict[str, Any]): variable dict
            label (str, optional): see FourVector. Defaults to "".
            duality (Duality, optional): see FourVector. Defaults to Duality.Contravariant.

        Raises:
            RuntimeError: failed to construct FourMomentum
        """
        for var in vars.keys():
            if var not in varname_set:
                raise RuntimeError(
                    f"{var} is not supported, cannot construct FourMomentum!"
                )

        if "pz" in vars.keys():
            op_list = [self._get_pz, self._get_pt, self._get_emass]
        elif "rapidity" in vars.keys():
            op_list = [self._get_rapidity, self._get_pz, self._get_pt, self._get_emass]
        elif "eta" in vars.keys() or "theta" in vars.keys():
            op_list = [self._get_eta, self._get_pt, self._get_pz, self._get_emass]
        else:
            raise RuntimeError(f"Can not construct FourMomemtum given {vars.keys()}")

        for op in op_list:
            op(vars)

        super().__init__(vars["E"], vars["px"], vars["py"], vars["pz"], label, duality)

    def _get_pz(self, vars: Dict[str, Any]):
        """Head, Leaf node

        Args:
            vars (Dict[str, Any]): vars
        """
        if "pz" in vars.keys():
            pass
        elif "eta" in vars.keys() and "pt" in vars.keys():
            vars["pz"] = pz_from_pt_eta(vars["pt"], vars["eta"])
        elif "rapidity" in vars.keys() and "E" in vars.keys():
            vars["pz"] = vars["E"] * np.cosh(vars["rapidity"])
        elif "rapidity" in vars.keys() and "m" in vars.keys():
            vars["pz"] = vars["m"] * np.sinh(vars["rapidity"])
        else:
            raise RuntimeError("Pz can not be determined!")

    def _get_rapidity(self, vars: Dict[str, Any]):
        """Head node

        Args:
            vars (Dict[str, Any]): _description_
        """
        if "rapidity" in vars.keys():
            pass

    def _get_eta(self, vars: Dict[str, Any]):
        """Head node

        Args:
            vars (Dict[str, Any]): _description_
        """
        if "eta" in vars.keys():
            pass
        elif "theta" in vars.keys():
            vars["eta"] = theta2eta(vars["theta"])

    def _get_pt(self, vars: Dict[str, Any]):
        """Leaf node

        Args:
            vars (Dict[str, Any]): vars
        """
        if "px" in vars.keys() and "py" in vars.keys():
            vars["pt"] = np.sqrt(vars["px"] ** 2 + vars["py"] ** 2)
            vars["phi"] = np.arctan2(vars["py"], vars["px"])
        elif "pt" in vars.keys() and "phi" in vars.keys():
            vars["px"] = vars["pt"] * np.cos(vars["phi"])
            vars["py"] = vars["pt"] * np.sin(vars["phi"])
        else:
            raise RuntimeError("Pt can not be determined!")

    def _get_emass(self, vars: Dict[str, Any]):
        """Leaf node

        Args:
            vars (Dict[str, Any]): vars
        """
        if "E" in vars.keys():
            pass
        elif "m" in vars.keys():
            px, py, pz = vars["px"], vars["py"], vars["pz"]
            p2 = px**2 + py**2 + pz**2
            m2 = vars["m"] ** 2
            vars["E"] = np.sqrt(m2 + p2)
        else:
            raise RuntimeError("E can not be determined!")
