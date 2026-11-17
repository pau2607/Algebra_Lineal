"""Tikhonov generalizado: penalización lambda ||L x||^2 en lugar de lambda ||x||^2.

    x_lambda = argmin ||A x - b||^2 + lambda ||L x||^2,       L de tamaño p x n.

Con L = I se recupera el Tikhonov de `metodos.py`. Con L de primeras diferencias
se penalizan los saltos entre puntos vecinos (oscilaciones), no el tamaño de x.

Unicidad: la función objetivo es estrictamente convexa si y solo si
ker A ∩ ker L = {0}, porque z^T (A^T A + lambda L^T L) z = ||A z||^2 + lambda ||L z||^2
solo se anula cuando A z = 0 y L z = 0 a la vez. Para L de diferencias
ker L = span{1} (constantes) y A 1 = 1 (filas de A suman 1), así que la condición
se cumple aunque A fuera singular.

Dos formas de calcular, que las pruebas comparan entre sí:

* `tikhonov_general_aumentado`: minimos cuadrados || [A; sqrt(lambda) L] x - [b; 0] ||
  resuelto con QR (sin formar A^T A).
* `GSVD` + `tikhonov_general_espectral`: descomposición en valores singulares
  generalizada (Paige-Saunders)

      A = U C X^{-1},   L = V S X^{-1},   C = diag(c_i), S = diag(s_i), c_i^2 + s_i^2 = 1,

  calculada con una QR de [A; L] seguida de una SVD del bloque superior de Q.
  En esa base la solución es

      x_lambda = sum_i  f_i(lambda) (u_i^T b / c_i) x_i,
      f_i(lambda) = c_i^2 / (c_i^2 + lambda s_i^2) = gamma_i^2 / (gamma_i^2 + lambda),

  con gamma_i = c_i / s_i los valores singulares generalizados. Las direcciones del
  núcleo de L (s_i = 0) tienen f_i = 1: nunca se penalizan.
"""

from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla


# --------------------------------------------------------------------------
# Matrices de penalización
# --------------------------------------------------------------------------

def primeras_diferencias(n: int) -> np.ndarray:
    """L_1 de tamaño (n-1) x n: (L x)_j = x_{j+1} - x_j. ker L_1 = constantes."""
    return np.diff(np.eye(n), axis=0)


def segundas_diferencias(n: int) -> np.ndarray:
    """L_2 de tamaño (n-2) x n: (L x)_j = x_j - 2 x_{j+1} + x_{j+2}. ker L_2 = rectas."""
    return np.diff(np.eye(n), n=2, axis=0)


def nucleos_disjuntos(A: np.ndarray, L: np.ndarray) -> bool:
    """¿ker A ∩ ker L = {0}?  Equivale a que [A; L] tenga columnas independientes."""
    M = np.vstack([A, L])
    return np.linalg.matrix_rank(M) == A.shape[1]


# --------------------------------------------------------------------------
# Sistema aumentado
# --------------------------------------------------------------------------

def tikhonov_general_aumentado(A: np.ndarray, L: np.ndarray, b: np.ndarray, lam: float) -> np.ndarray:
    """min || [A; sqrt(lambda) L] x - [b; 0] ||^2 con QR."""
    if lam <= 0:
        raise ValueError("lambda debe ser positivo")
    M = np.vstack([A, np.sqrt(lam) * L])
    rhs = np.concatenate([b, np.zeros(L.shape[0])])
    Q, R = np.linalg.qr(M, mode="reduced")
    return sla.solve_triangular(R, Q.T @ rhs, lower=False)


# --------------------------------------------------------------------------
# GSVD y fórmula espectral
# --------------------------------------------------------------------------

@dataclass
class GSVD:
    """A = U diag(c) X^{-1}, L = V diag(s) X^{-1} (V puede faltar columnas donde s_i = 0)."""
    U: np.ndarray
    c: np.ndarray
    s: np.ndarray
    X: np.ndarray

    @classmethod
    def de(cls, A: np.ndarray, L: np.ndarray) -> "GSVD":
        m, n = A.shape
        if not nucleos_disjuntos(A, L):
            raise ValueError("ker A ∩ ker L ≠ {0}: la GSVD así definida requiere [A; L] de rango completo")
        Q, R = np.linalg.qr(np.vstack([A, L]), mode="reduced")   # (m+p) x n, R invertible
        Q1, Q2 = Q[:m], Q[m:]                                     # A = Q1 R,  L = Q2 R
        U, c, Wt = np.linalg.svd(Q1, full_matrices=False)         # Q1 = U diag(c) W^T, 0 <= c <= 1
        # Q1^T Q1 + Q2^T Q2 = I  =>  las columnas de Q2 W son ortogonales con normas s_i = sqrt(1 - c_i^2).
        # Se calculan directamente (sqrt(1 - c^2) pierde precisión cuando c ≈ 1, en las direcciones de ker L).
        s = np.linalg.norm(Q2 @ Wt.T, axis=0)
        # Tolerancia de rango (como en pinv): un s_i por debajo es una dirección de ker L con redondeo.
        s[s < 100 * max(m + L.shape[0], n) * np.finfo(float).eps] = 0.0
        X = sla.solve_triangular(R, Wt.T, lower=False)            # X = R^{-1} W  =>  X^{-1} = W^T R
        return cls(U=U, c=c, s=s, X=X)

    @property
    def gamma(self) -> np.ndarray:
        """Valores singulares generalizados gamma_i = c_i / s_i (inf donde s_i = 0)."""
        with np.errstate(divide="ignore"):
            return np.where(self.s > 0, self.c / np.where(self.s > 0, self.s, 1.0), np.inf)

    def coeficientes(self, b: np.ndarray) -> np.ndarray:
        return self.U.T @ b


def factores_filtro_general(gsvd: GSVD, lam: float) -> np.ndarray:
    """f_i = c_i^2 / (c_i^2 + lambda s_i^2)."""
    return gsvd.c ** 2 / (gsvd.c ** 2 + lam * gsvd.s ** 2)


def tikhonov_general_espectral(gsvd: GSVD, b: np.ndarray, lam: float) -> np.ndarray:
    """x_lambda = X diag(c_i / (c_i^2 + lambda s_i^2)) U^T b."""
    if lam <= 0:
        raise ValueError("lambda debe ser positivo")
    coef = gsvd.c / (gsvd.c ** 2 + lam * gsvd.s ** 2) * gsvd.coeficientes(b)
    return gsvd.X @ coef


# --------------------------------------------------------------------------
# Selección del parámetro (misma regla de discrepancia; no usa x_true)
# --------------------------------------------------------------------------

def discrepancia_general(A, L, b, lambdas, delta: float, tau: float = 1.1, gsvd: GSVD | None = None):
    """Mayor lambda de la malla con ||A x_lambda - b|| <= tau * delta. Devuelve (lambda, idx) o (None, None)."""
    lambdas = np.asarray(lambdas)
    if gsvd is None:
        gsvd = GSVD.de(A, L)
    residuos = np.array([np.linalg.norm(A @ tikhonov_general_espectral(gsvd, b, lam) - b) for lam in lambdas])
    viables = np.flatnonzero(residuos <= tau * delta)
    if viables.size == 0:
        return None, None
    idx = viables[np.argmax(lambdas[viables])]
    return float(lambdas[idx]), int(idx)
