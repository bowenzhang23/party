import numpy as np

from typing import Any, List
from pyminiroot import MinirootBase


class DataType:
    """Data Types"""

    Float = 0
    Double = 1
    Byte = 2
    Short = 3
    Integer = 4
    Long = 5


class Miniroot(MinirootBase):
    def __init__(self, file_name: str):
        """Inherit from MinirootBase, extend get functionality

        Args:
            file_name (str): Path to root file
        """
        super().__init__(file_name)
        self.func = {
            DataType.Float: self.get_float,
            DataType.Double: self.get_double,
            DataType.Byte: self.get_byte,
            DataType.Short: self.get_short,
            DataType.Integer: self.get_integer,
            DataType.Long: self.get_long,
        }
        self.dtype = {
            DataType.Float: np.float32,
            DataType.Double: np.float64,
            DataType.Byte: np.byte,
            DataType.Short: np.int16,
            DataType.Integer: np.int32,
            DataType.Long: np.int64,
        }

    def branches(self) -> List[str]:
        """
        Get names of branches

        Returns:
            List[str]: List of branch names
        """
        return super().branches()

    def get(self, branch_name: str, type: DataType) -> Any:
        """Get numpy array specifying branch name and data type

        Args:
            branch_name (str): name of the branch
            type (DataType): data type to interpret to

        Returns:
            Any: numpy array
        """
        return np.array(self.func[type](branch_name), self.dtype[type])
