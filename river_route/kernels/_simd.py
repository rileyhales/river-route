"""
Four float32 lanes in one SIMD register for the numba kernels: a numba type for them and the LLVM operations on it that
numba has no way to write, such as moving values between lanes, which lets a recurrence be solved four steps at a time.
simd.h names its operations the same.
"""

from collections.abc import Callable
from typing import Any, Concatenate, cast

from llvmlite import ir
from numba import types
from numba.core import cgutils
from numba.extending import intrinsic as numba_intrinsic
from numba.extending import models, register_model

__all__ = [
    'FourLanes',
    'add',
    'broadcast',
    'four_lanes',
    'last_lane',
    'last_lane_in_every_lane',
    'load_four',
    'multiply',
    'multiply_add',
    'shift_in_one_lane',
    'shift_in_two_lanes',
    'store_four',
]

LANE = ir.FloatType()
VECTOR = ir.VectorType(LANE, 4)
LANE_NUMBER = ir.IntType(32)


class FourLanes(types.Type):
    """The numba type of four float32 values held in one SIMD register."""

    def __init__(self) -> None:
        super().__init__(name='FourLanes')


four_lanes_type = FourLanes()


@register_model(FourLanes)
class _FourLanesModel(models.PrimitiveModel):
    def __init__(self, dmm, fe_type) -> None:
        super().__init__(dmm, fe_type, VECTOR)


def intrinsic[**P](definition: Callable[Concatenate[Any, P], Any]) -> Callable[P, Any]:
    """numba's intrinsic decorator, typed as kernels call the intrinsic: without the typing context numba passes."""
    return cast(Callable[P, Any], numba_intrinsic(definition))


def _is_float32_row(values: types.Type) -> bool:
    """Whether numba typed values as a 1D float32 array. Every row routing reads is checked to have its values
    contiguous before routing, which the lane loads and stores rely on, since an 'A' layout row does not show it."""
    return isinstance(values, types.Array) and values.ndim == 1 and values.dtype == types.float32


def _lane_mask(*lanes: int) -> ir.Constant:
    """The lanes a shuffle takes, numbered 0 to 3 in its first vector and 4 to 7 in its second."""
    return ir.Constant(ir.VectorType(LANE_NUMBER, 4), list(lanes))


@intrinsic
def load_four(typingctx, values, index):
    """values[index] to values[index + 3] in one load, from anywhere in a row whose values are contiguous."""
    if not _is_float32_row(values) or not isinstance(index, types.Integer):
        return None

    def codegen(context, builder, signature, args):
        data = context.make_array(signature.args[0])(context, builder, args[0]).data
        return builder.load(builder.bitcast(builder.gep(data, [args[1]]), VECTOR.as_pointer()), align=4)

    return four_lanes_type(values, index), codegen


@intrinsic
def store_four(typingctx, values, index, vector):
    """Write the four lanes into values[index] to values[index + 3] in one store."""
    if not _is_float32_row(values) or not isinstance(index, types.Integer) or vector != four_lanes_type:
        return None

    def codegen(context, builder, signature, args):
        data = context.make_array(signature.args[0])(context, builder, args[0]).data
        builder.store(args[2], builder.bitcast(builder.gep(data, [args[1]]), VECTOR.as_pointer()), align=4)
        return context.get_dummy_value()

    return types.none(values, index, vector), codegen


@intrinsic
def broadcast(typingctx, value):
    """The float32 value in all four lanes."""
    if value != types.float32:
        return None

    def codegen(context, builder, signature, args):
        first = builder.insert_element(ir.Constant(VECTOR, ir.Undefined), args[0], ir.Constant(LANE_NUMBER, 0))
        return builder.shuffle_vector(first, first, _lane_mask(0, 0, 0, 0))

    return four_lanes_type(value), codegen


@intrinsic
def four_lanes(typingctx, first, second, third, fourth):
    """Four float32 values, one per lane in order."""
    if any(value != types.float32 for value in (first, second, third, fourth)):
        return None

    def codegen(context, builder, signature, args):
        vector = ir.Constant(VECTOR, ir.Undefined)
        for lane_number, value in enumerate(args):
            vector = builder.insert_element(vector, value, ir.Constant(LANE_NUMBER, lane_number))
        return vector

    return four_lanes_type(first, second, third, fourth), codegen


@intrinsic
def multiply_add(typingctx, left, right, addend):
    """left * right + addend in each lane, rounded once."""
    if any(vector != four_lanes_type for vector in (left, right, addend)):
        return None

    def codegen(context, builder, signature, args):
        function = cgutils.get_or_insert_function(
            builder.module, ir.FunctionType(VECTOR, [VECTOR] * 3), 'llvm.fma.v4f32'
        )
        return builder.call(function, args)

    return four_lanes_type(left, right, addend), codegen


@intrinsic
def add(typingctx, left, right):
    """left + right in each lane."""
    if left != four_lanes_type or right != four_lanes_type:
        return None

    def codegen(context, builder, signature, args):
        return builder.fadd(*args)

    return four_lanes_type(left, right), codegen


@intrinsic
def multiply(typingctx, left, right):
    """left * right in each lane."""
    if left != four_lanes_type or right != four_lanes_type:
        return None

    def codegen(context, builder, signature, args):
        return builder.fmul(*args)

    return four_lanes_type(left, right), codegen


@intrinsic
def shift_in_one_lane(typingctx, before, vector):
    """Each lane of vector moved up one lane, and the last lane of before in the first: (b3, v0, v1, v2)."""
    if before != four_lanes_type or vector != four_lanes_type:
        return None

    def codegen(context, builder, signature, args):
        return builder.shuffle_vector(args[0], args[1], _lane_mask(3, 4, 5, 6))

    return four_lanes_type(before, vector), codegen


@intrinsic
def shift_in_two_lanes(typingctx, before, vector):
    """Each lane of vector moved up two lanes, and the last two lanes of before in the first two: (b2, b3, v0, v1)."""
    if before != four_lanes_type or vector != four_lanes_type:
        return None

    def codegen(context, builder, signature, args):
        return builder.shuffle_vector(args[0], args[1], _lane_mask(2, 3, 4, 5))

    return four_lanes_type(before, vector), codegen


@intrinsic
def last_lane_in_every_lane(typingctx, vector):
    """The value of the last lane in all four lanes."""
    if vector != four_lanes_type:
        return None

    def codegen(context, builder, signature, args):
        return builder.shuffle_vector(args[0], args[0], _lane_mask(3, 3, 3, 3))

    return four_lanes_type(vector), codegen


@intrinsic
def last_lane(typingctx, vector):
    """The float32 value of the last lane."""
    if vector != four_lanes_type:
        return None

    def codegen(context, builder, signature, args):
        return builder.extract_element(args[0], ir.Constant(LANE_NUMBER, 3))

    return types.float32(vector), codegen
