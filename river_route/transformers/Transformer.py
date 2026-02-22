from .AbstractTransformer import AbstractBaseTransformer
from ..types import FloatArray

__all__ = ['Transformer', ]


class Transformer(AbstractBaseTransformer):
    """
    Concrete transformer initialized exclusively via AbstractBaseTransformer.from_kernel().

    This class provides no kernel-computation logic of its own.  Its only purpose is to give
    a concrete, instantiable type for pre-computed kernels so that AbstractBaseTransformer can remain
    a proper ABC.  Direct instantiation (``Transformer(dt)``) is intentionally not possible —
    load a transformer from a cached parquet file instead::

    t = (
        Transformer
        .from_kernel(dt=3600, kernel='kernel.parquet')
        .set_state('state.parquet')
    )
    """

    def _build_kernel(self) -> FloatArray:
        raise NotImplementedError(
            'Transformer should only be instantiated with Transformer.from_kernel() with your pre-computed kernel.'
        )
