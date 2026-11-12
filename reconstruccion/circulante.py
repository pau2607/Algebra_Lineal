"""Suavizado circular de una señal larga (audio) y su SVD vía la transformada de Fourier.

Para un audio de N = 16 384 muestras no es razonable formar A (N x N) ni calcular
su SVD. Si el suavizado es una convolución circular con un núcleo h,

    A x = h ⊛ x,   A = circulante(h),

la matriz es circulante y la transformada discreta de Fourier la diagonaliza:
A = F^{-1} diag(ĥ) F con ĥ = fft(h). Si h es simétrico (h_j = h_{N-j}), ĥ es real
y |ĥ_k| son los valores singulares de A (con vectores singulares de Fourier).
Las fórmulas espectrales del caso general quedan, frecuencia por frecuencia:

    pseudoinversa:  x̂_k = b̂_k / ĥ_k
    Tikhonov:       x̂_k = ĥ_k / (ĥ_k^2 + lambda) * b̂_k

y el residuo podría calcularse en frecuencia por la identidad de Parseval,
||r|| = ||r̂|| / sqrt(N); aquí se calcula en el tiempo, que da lo mismo.
`tests/test_aplicaciones.py` comprueba que esto coincide con la SVD de la
matriz circulante formada explícitamente para N pequeño.
"""

from dataclasses import dataclass

import numpy as np


def nucleo_gaussiano_circular(N: int, sigma_muestras: float) -> np.ndarray:
    """h_j = exp(-d_j^2 / (2 sigma^2)) / suma, con d_j = min(j, N-j) (distancia circular).

    Suma 1 (conserva constantes, como las filas de `datos.matriz_suavizado`) y es simétrico.
    """
    j = np.arange(N)
    d = np.minimum(j, N - j)
    h = np.exp(-(d ** 2) / (2.0 * sigma_muestras ** 2))
    return h / h.sum()


@dataclass
class OperadorCirculante:
    h: np.ndarray

    @property
    def N(self) -> int:
        return self.h.size

    @property
    def h_hat(self) -> np.ndarray:
        """Autovalores ĥ_k (reales para núcleo simétrico)."""
        return np.fft.fft(self.h).real

    @property
    def sigma(self) -> np.ndarray:
        """Valores singulares |ĥ_k| (en el orden de las frecuencias, no ordenados)."""
        return np.abs(self.h_hat)

    @property
    def sigma1(self) -> float:
        return float(self.sigma.max())

    def frecuencias(self, fs: float) -> np.ndarray:
        """Frecuencia en Hz asociada a cada k."""
        return np.fft.fftfreq(self.N, d=1.0 / fs)

    def aplicar(self, x: np.ndarray) -> np.ndarray:
        return np.fft.ifft(self.h_hat * np.fft.fft(x)).real

    def pseudoinversa_completa(self, b: np.ndarray) -> np.ndarray:
        """Divide cada coeficiente de Fourier entre ĥ_k, en todas las frecuencias con ĥ_k != 0
        (sin umbral, como `metodos.pseudoinversa_completa`). Amplifica el ruido."""
        hh = self.h_hat
        no_nulos = hh != 0
        cociente = np.where(no_nulos, np.fft.fft(b) / np.where(no_nulos, hh, 1.0), 0.0)
        return np.fft.ifft(cociente).real

    def tikhonov(self, b: np.ndarray, lam: float) -> np.ndarray:
        if lam <= 0:
            raise ValueError("lambda debe ser positivo")
        hh = self.h_hat
        return np.fft.ifft(hh / (hh ** 2 + lam) * np.fft.fft(b)).real

    def residuo(self, x: np.ndarray, b: np.ndarray) -> float:
        return float(np.linalg.norm(self.aplicar(x) - b))


def discrepancia_tikhonov(op: OperadorCirculante, b, lambdas, delta: float, tau: float = 1.1):
    """Mayor lambda de la malla con ||A x_lambda - b|| <= tau * delta (no usa la verdad)."""
    lambdas = np.asarray(lambdas)
    residuos = np.array([op.residuo(op.tikhonov(b, lam), b) for lam in lambdas])
    viables = np.flatnonzero(residuos <= tau * delta)
    if viables.size == 0:
        return None, None
    idx = viables[np.argmax(lambdas[viables])]
    return float(lambdas[idx]), int(idx)


def snr_db(x: np.ndarray, x_true: np.ndarray) -> float:
    """10 log10(||x_true||^2 / ||x - x_true||^2): más alto = más parecido a la señal original."""
    return float(10.0 * np.log10(np.sum(x_true ** 2) / np.sum((x - x_true) ** 2)))
