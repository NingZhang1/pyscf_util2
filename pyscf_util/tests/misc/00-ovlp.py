import unittest
import numpy as np
from pyscf_util.misc.math import _orthogonalize


class TestOrthogonalize(unittest.TestCase):
    def test_orthogonalize_with_identity_overlap(self):
        # Test with identity matrix as overlap
        vec = np.array([[1, 1], [1, 2], [1, 3]])
        result = _orthogonalize(vec)

        # Check if the result is orthogonal
        self.assertTrue(np.allclose(result.T @ result, np.eye(2), atol=1e-7))

        # Check if the result is still in the span of original vectors
        self.assertTrue(
            np.allclose(np.linalg.matrix_rank(vec), np.linalg.matrix_rank(result))
        )

    def test_orthogonalize_with_custom_overlap(self):
        # Test with a custom overlap matrix
        vec = np.array([[1, 1], [1, 2], [1, 3]])
        ovlp = np.array([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
        result = _orthogonalize(vec, ovlp)

        # Check if the result is orthogonal under the custom inner product
        self.assertTrue(np.allclose(result.T @ ovlp @ result, np.eye(2), atol=1e-7))

        # Check if the result is still in the span of original vectors
        self.assertTrue(
            np.allclose(np.linalg.matrix_rank(vec), np.linalg.matrix_rank(result))
        )

    def test_orthogonalize_with_random_matrix(self):
        # Test with a random matrix
        np.random.seed(42)  # For reproducibility
        vec = np.random.rand(5, 3)
        result = _orthogonalize(vec)

        # Check if the result is orthogonal
        self.assertTrue(np.allclose(result.T @ result, np.eye(3), atol=1e-7))

        # Check if the result is still in the span of original vectors
        self.assertTrue(
            np.allclose(np.linalg.matrix_rank(vec), np.linalg.matrix_rank(result))
        )


if __name__ == "__main__":
    unittest.main()
