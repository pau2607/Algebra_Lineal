"""Reconstrucción con operadores separables: imágenes y audio por bloques.

Modelo de observación (X es una matriz de datos, no un vector):

    B = A_v X A_h^T + E

* Imagen: A_v desenfoca en la dirección vertical (mezcla filas) y A_h en la
  horizontal (mezcla columnas). Es un desenfoque gaussiano 2D separable.
* Audio por bloques: cada columna de X es un bloque de muestras y A_h = I,
  así que B = A_v X + E es el problema 1D aplicado a todas las columnas a la vez.

Vectorizando por columnas, vec(A_v X A_h^T) = (A_h ⊗ A_v) vec(X), y la SVD del
producto de Kronecker se obtiene de las SVD de los factores:

    A_h ⊗ A_v = (U_h ⊗ U_v) (Σ_h ⊗ Σ_v) (V_h ⊗ V_v)^T.

Por eso nunca se forma la matriz grande (para una imagen de 256 x 256 sería de
65 536 x 65 536): todas las fórmulas espectrales del caso 1D se aplican entrada
por entrada sobre la matriz de coeficientes C = U_v^T B U_h, con valores
singulares S_ij = sigma_i(A_v) * sigma_j(A_h).
"""

from dataclasses import dataclass
from functools import cached_property

import numpy as np


@dataclass
class OperadorSeparable:
    """B = A_v X A_h^T con las SVD de los dos factores precalculadas."""
    A_v: np.ndarray
    A_h: np.ndarray
    U_v: np.ndarray
    s_v: np.ndarray
    V_v: np.ndarray
    U_h: np.ndarray
    s_h: np.ndarray
    V_h: np.ndarray

    @classmethod
    def de(cls, A_v: np.ndarray, A_h: np.ndarray) -> "OperadorSeparable":
        U_v, s_v, Vt_v = np.linalg.svd(A_v, full_matrices=False)
        U_h, s_h, Vt_h = np.linalg.svd(A_h, full_matrices=False)
        return cls(A_v, A_h, U_v, s_v, Vt_v.T, U_h, s_h, Vt_h.T)

    # --- operador y espectro ------------------------------------------------

    def aplicar(self, X: np.ndarray) -> np.ndarray:
        """A_v X A_h^T."""
        return self.A_v @ X @ self.A_h.T

    @cached_property
    def S(self) -> np.ndarray:
        """Valores singulares del operador completo, S_ij = sigma_i(A_v) sigma_j(A_h)."""
        return np.outer(self.s_v, self.s_h)

    @property
    def sigma1(self) -> float:
        return float(self.s_v[0] * self.s_h[0])

    def coeficientes(self, B: np.ndarray) -> np.ndarray:
        """C = U_v^T B U_h (coeficientes de Picard en forma matricial)."""
        return self.U_v.T @ B @ self.U_h

    def _sintetizar(self, Z: np.ndarray) -> np.ndarray:
        """V_v Z V_h^T: vuelve de coordenadas espectrales al espacio de la imagen."""
        return self.V_v @ Z @ self.V_h.T

    # --- métodos de reconstrucción --------------------------------------------

    def rango_numerico(self) -> int:
        """Cuántos S_ij superan el umbral por defecto de scipy.linalg.pinv, max(m, n) eps S_max,
        con (m, n) el tamaño de la matriz de Kronecker A_h ⊗ A_v (que no se forma)."""
        m = self.A_v.shape[0] * self.A_h.shape[0]
        n = self.A_v.shape[1] * self.A_h.shape[1]
        umbral = max(m, n) * np.finfo(float).eps * self.S.max()
        return int((self.S > umbral).sum())

    def pseudoinversa_completa(self, B: np.ndarray) -> np.ndarray:
        """Divide cada coeficiente entre su S_ij (todos los S_ij > 0, sin umbral)."""
        S = self.S
        C = self.coeficientes(B)
        Z = np.where(S > 0, C / np.where(S > 0, S, 1.0), 0.0)
        return self._sintetizar(Z)

    @cached_property
    def orden_espectral(self) -> np.ndarray:
        """Índices (aplanados) de S ordenados de mayor a menor."""
        return np.argsort(self.S, axis=None, kind="stable")[::-1]

    def mascara_truncada(self, k: int) -> np.ndarray:
        """Máscara que conserva los k valores singulares S_ij más grandes."""
        S = self.S
        if not 0 <= k <= S.size:
            raise ValueError(f"k = {k} fuera de rango")
        m = np.zeros(S.size, dtype=bool)
        m[self.orden_espectral[:k]] = True
        return m.reshape(S.shape)

    def svd_truncada(self, B: np.ndarray, k: int) -> np.ndarray:
        """X_k: conserva solo las k direcciones con mayor S_ij."""
        S = self.S
        m = self.mascara_truncada(k)
        Z = np.where(m, self.coeficientes(B) / np.where(m, S, 1.0), 0.0)
        return self._sintetizar(Z)

    def tikhonov(self, B: np.ndarray, lam: float) -> np.ndarray:
        """X_lambda = V_v [ S/(S^2+lambda) ∘ C ] V_h^T  (∘: producto entrada por entrada)."""
        if lam <= 0:
            raise ValueError("lambda debe ser positivo")
        S = self.S
        return self._sintetizar(S / (S ** 2 + lam) * self.coeficientes(B))

    # --- medidas sin x_true ----------------------------------------------------

    def residuo(self, X: np.ndarray, B: np.ndarray) -> float:
        """||A_v X A_h^T - B||_F (norma euclídea del vector vec)."""
        return float(np.linalg.norm(self.aplicar(X) - B))

    def residuo_tikhonov_espectral(self, B: np.ndarray, lam: float) -> float:
        """Residuo de Tikhonov sin reconstruir: ||(lambda/(S^2+lambda)) ∘ C||_F.

        Válido cuando U_v y U_h son cuadradas (A_v, A_h cuadradas), como aquí.
        """
        S = self.S
        return float(np.linalg.norm(lam / (S ** 2 + lam) * self.coeficientes(B)))


# --- selección del parámetro (no usa la verdad) -----------------------------------

def discrepancia_tikhonov(op: OperadorSeparable, B, lambdas, delta: float, tau: float = 1.1):
    """Mayor lambda de la malla con ||A_v X_lambda A_h^T - B||_F <= tau * delta.

    Devuelve (lambda, índice) o (None, None) si ningún lambda es viable.
    """
    lambdas = np.asarray(lambdas)
    residuos = np.array([op.residuo(op.tikhonov(B, lam), B) for lam in lambdas])
    viables = np.flatnonzero(residuos <= tau * delta)
    if viables.size == 0:
        return None, None
    idx = viables[np.argmax(lambdas[viables])]
    return float(lambdas[idx]), int(idx)


def discrepancia_tsvd(op: OperadorSeparable, B, ks, delta: float, tau: float = 1.1):
    """Menor k de la malla con residuo <= tau * delta (más regularización posible)."""
    ks = np.asarray(ks)
    residuos = np.array([op.residuo(op.svd_truncada(B, int(k)), B) for k in ks])
    viables = np.flatnonzero(residuos <= tau * delta)
    if viables.size == 0:
        return None, None
    idx = viables[np.argmin(ks[viables])]
    return int(ks[idx]), int(idx)


# --- utilidades ----------------------------------------------------------------

def ruido_matricial(B_limpio: np.ndarray, eta: float, semilla: int) -> np.ndarray:
    """E con ||E||_F / ||B_limpio||_F = eta exactamente (misma receta que el caso 1D)."""
    Z = np.random.default_rng(semilla).standard_normal(B_limpio.shape)
    return eta * np.linalg.norm(B_limpio) * Z / np.linalg.norm(Z)


def psnr(X: np.ndarray, X_true: np.ndarray, rango: float = 1.0) -> float:
    """Relación señal/ruido de pico en dB, medida usual en imágenes (más alto = mejor)."""
    mse = float(np.mean((X - X_true) ** 2))
    return float("inf") if mse == 0 else 10.0 * np.log10(rango ** 2 / mse)
