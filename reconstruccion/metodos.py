"""Métodos de reconstrucción.

Convención: en Tikhonov, lambda es directamente el coeficiente de la penalización,
    x_lambda = argmin ||A x - b||^2 + lambda ||x||^2,
igual que en el artículo de referencia (no se usa lambda^2).
"""

from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla


@dataclass
class SVD:
    """SVD compacta A = U diag(sigma) V^T, con sigma en orden decreciente."""
    U: np.ndarray
    sigma: np.ndarray
    V: np.ndarray

    @classmethod
    def de(cls, A: np.ndarray) -> "SVD":
        U, sigma, Vt = np.linalg.svd(A, full_matrices=False)
        return cls(U=U, sigma=sigma, V=Vt.T)

    def coeficientes(self, b: np.ndarray) -> np.ndarray:
        """Coeficientes de Picard u_i^T b."""
        return self.U.T @ b


# Pseudoinversa


def tolerancia_rango_scipy(A: np.ndarray, sigma1: float) -> float:
    """Umbral por defecto de scipy.linalg.pinv: rtol = max(m, n) * eps, aplicado como rtol * sigma_1."""
    m, n = A.shape
    return max(m, n) * np.finfo(A.dtype).eps * sigma1


def pseudoinversa_numerica(A: np.ndarray, b: np.ndarray):
    """x = pinv(A) b con la tolerancia por defecto de SciPy.

    Devuelve (x, rango_numerico, umbral). La función descarta los valores singulares
    menores que el umbral, así que NO es una 'inversa sin regularizar': el rango que
    usa debe informarse.
    """
    A_pinv, rango = sla.pinv(A, return_rank=True)
    sigma1 = np.linalg.norm(A, 2)
    return A_pinv @ b, int(rango), tolerancia_rango_scipy(A, sigma1)


def pseudoinversa_completa(svd: SVD, b: np.ndarray) -> np.ndarray:
    """x = sum_i (u_i^T b / sigma_i) v_i sobre TODOS los sigma_i > 0 (sin umbral).

    Sirve para mostrar la amplificación del ruido; no es un método recomendable.
    """
    positivos = svd.sigma > 0
    c = svd.coeficientes(b)
    return svd.V[:, positivos] @ (c[positivos] / svd.sigma[positivos])


# SVD truncada


def svd_truncada(svd: SVD, b: np.ndarray, k: int) -> np.ndarray:
    """x_k = sum_{i<=k} (u_i^T b / sigma_i) v_i."""
    if not 0 <= k <= len(svd.sigma):
        raise ValueError(f"k = {k} fuera de rango")
    c = svd.coeficientes(b)
    return svd.V[:, :k] @ (c[:k] / svd.sigma[:k])



# Tikhonov


def factores_filtro(sigma: np.ndarray, lam: float) -> np.ndarray:
    """f_i(lambda) = sigma_i^2 / (sigma_i^2 + lambda)."""
    return sigma ** 2 / (sigma ** 2 + lam)


def tikhonov_espectral(svd: SVD, b: np.ndarray, lam: float) -> np.ndarray:
    """x_lambda = sum_i sigma_i / (sigma_i^2 + lambda) (u_i^T b) v_i.

    Se usa sigma/(sigma^2+lambda) en lugar de f_i/sigma_i para no dividir entre
    valores singulares diminutos.
    """
    if lam <= 0:
        raise ValueError("lambda debe ser positivo")
    c = svd.coeficientes(b)
    return svd.V @ (svd.sigma / (svd.sigma ** 2 + lam) * c)


def tikhonov_aumentado(A: np.ndarray, b: np.ndarray, lam: float) -> np.ndarray:
    """Resuelve min || [A; sqrt(lambda) I] x - [b; 0] ||^2 con QR (sin ecuaciones normales)."""
    m, n = A.shape
    M = np.vstack([A, np.sqrt(lam) * np.eye(n)])
    rhs = np.concatenate([b, np.zeros(n)])
    Q, R = np.linalg.qr(M, mode="reduced")
    return sla.solve_triangular(R, Q.T @ rhs, lower=False)


def tikhonov_normales(A: np.ndarray, b: np.ndarray, lam: float) -> np.ndarray:
    """(A^T A + lambda I)^{-1} A^T b. Solo para comprobaciones con lambda moderado."""
    n = A.shape[1]
    return np.linalg.solve(A.T @ A + lam * np.eye(n), A.T @ b)
