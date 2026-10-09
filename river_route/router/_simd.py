"""
Four float32 lanes in one SIMD register for the numba kernels: a numba type for them and the LLVM operations on it that
numba has no way to write, such as moving values between lanes, which lets a recurrence be solved four steps at a time.
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
    'last_lane_in_every_lane',
    'load_up_to_four',
    'multiply',
    'multiply_add',
    'shift_in_one_lane',
    'shift_in_two_lanes',
    'store_up_to_four',
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


def _pointer_to_lane(builder, data, index, lane_number: int):
    """A pointer to values[index + lane_number], given the pointer to the row's first value."""
    return builder.gep(data, [builder.add(index, ir.Constant(index.type, lane_number))])


def _lane_mask(*lanes: int) -> ir.Constant:
    """The lanes a shuffle takes, numbered 0 to 3 in its first vector and 4 to 7 in its second."""
    return ir.Constant(ir.VectorType(LANE_NUMBER, 4), list(lanes))


@intrinsic
def load_up_to_four(typingctx, values, index, count):
    """
    values[index : index + count] in the first count lanes and 0 in the rest: four lanes in one load unless a series
    ends, where each of its last values is loaded on its own. Called with a count of 4 that the compiler knows, the
    choice between them is compiled away.
    """
    if not _is_float32_row(values) or not isinstance(index, types.Integer) or not isinstance(count, types.Integer):
        return None

    def codegen(context, builder, signature, args):
        row, index, count = args
        data = context.make_array(signature.args[0])(context, builder, row).data
        lanes = cgutils.alloca_once_value(builder, ir.Constant(VECTOR, [0.0] * 4))
        with builder.if_else(builder.icmp_signed('>=', count, ir.Constant(count.type, 4))) as (four, fewer):
            with four:
                pointer = builder.bitcast(builder.gep(data, [index]), VECTOR.as_pointer())
                builder.store(builder.load(pointer, align=4), lanes)
            with fewer:
                for lane_number in range(3):
                    with builder.if_then(builder.icmp_signed('>', count, ir.Constant(count.type, lane_number))):
                        value = builder.load(_pointer_to_lane(builder, data, index, lane_number))
                        lane = ir.Constant(LANE_NUMBER, lane_number)
                        builder.store(builder.insert_element(builder.load(lanes), value, lane), lanes)
        return builder.load(lanes)

    return four_lanes_type(values, index, count), codegen


@intrinsic
def store_up_to_four(typingctx, values, index, vector, count):
    """
    Write the first count lanes into values[index : index + count]: four lanes in one store unless a series ends, where
    each of its last values is stored on its own. Called with a count of 4 that the compiler knows, the choice between
    them is compiled away.
    """
    if not _is_float32_row(values) or not isinstance(index, types.Integer) or vector != four_lanes_type:
        return None
    if not isinstance(count, types.Integer):
        return None

    def codegen(context, builder, signature, args):
        row, index, vector, count = args
        data = context.make_array(signature.args[0])(context, builder, row).data
        with builder.if_else(builder.icmp_signed('>=', count, ir.Constant(count.type, 4))) as (four, fewer):
            with four:
                builder.store(vector, builder.bitcast(builder.gep(data, [index]), VECTOR.as_pointer()), align=4)
            with fewer:
                for lane_number in range(3):
                    with builder.if_then(builder.icmp_signed('>', count, ir.Constant(count.type, lane_number))):
                        value = builder.extract_element(vector, ir.Constant(LANE_NUMBER, lane_number))
                        builder.store(value, _pointer_to_lane(builder, data, index, lane_number))
        return context.get_dummy_value()

    return types.none(values, index, vector, count), codegen


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
