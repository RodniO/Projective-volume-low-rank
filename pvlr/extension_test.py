import numpy as np

from pvlr.extension import IntVector, Matrix, maxvol


def test_maxvol(size=1_000, rank=10, max_steps=2):
    max_swaps = 4 * rank
    rng = np.random.default_rng(42)
    A = Matrix.from_buffer(rng.uniform(-1, 1, size=(size, size)))
    C = Matrix()
    CA = Matrix()
    AR = Matrix()

    index = np.arange(size, dtype=np.int32) + 1
    ix = IntVector.from_buffer(index)
    jx = IntVector.from_buffer(index)
    maxvol(A, size, size, rank, ix, jx, C, AR, max_steps, max_swaps, CA, True)

    arr_C = C.to_buffer()
    assert isinstance(arr_C, np.ndarray)
    assert arr_C.shape == (size, rank)

    arr_CA = CA.to_buffer()
    assert isinstance(arr_CA, np.ndarray)
    assert arr_CA.shape == (size, rank)

    arr_AR = AR.to_buffer()
    assert isinstance(arr_AR, np.ndarray)
    assert arr_AR.shape == (rank, size)

    arr_A = A.to_buffer()
    assert isinstance(arr_A, np.ndarray)
    assert arr_A.shape == (size, size)

    diff = arr_A - arr_C @ arr_AR
    aerr = np.linalg.norm(diff)
    print(aerr)
    breakpoint()
