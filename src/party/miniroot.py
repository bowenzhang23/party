import numpy as np

from typing import Any, Tuple
from pyminiroot import MinirootBase, DataType


class Miniroot(MinirootBase):
    def __init__(self, file_name):
        super().__init__(file_name)

    def get(self, branch_name: str, type: DataType) -> Any:
        func = {
            DataType.Float: self.get_float,
            DataType.Double: self.get_double,
            DataType.Byte: self.get_byte,
            DataType.Short: self.get_short,
            DataType.Integer: self.get_integer,
            DataType.Long: self.get_long,
        }
        dtype = {
            DataType.Float: np.float32,
            DataType.Double: np.float64,
            DataType.Byte: np.byte,
            DataType.Short: np.int16,
            DataType.Integer: np.int32,
            DataType.Long: np.int64,
        }
        return np.array(func[type](branch_name), dtype[type])
