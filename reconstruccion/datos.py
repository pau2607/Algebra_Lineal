
import numpy as np


def malla(n: int = 128) -> np.ndarray:
    """Puntos t_j = j/(n-1), j = 0, ..., n-1."""
    return np.arange(n) / (n - 1)


def senal_verdadera(t: np.ndarray) -> np.ndarray:
    """x_j = exp[-80 (t_j - 0.3)^2] + 0.6 * 1_[0.6, 0.75](t_j).

    Una región suave (gaussiana centrada en 0.3) y un pulso rectangular en [0.6, 0.75].
    """
    suave = np.exp(-80.0 * (t - 0.3) ** 2)
    pulso = 0.6 * ((t >= 0.6) & (t <= 0.75)).astype(float)
    return suave + pulso


def matriz_suavizado(t: np.ndarray, s: float = 0.03) -> np.ndarray:
    """Suavizado gaussiano discreto con filas normalizadas.

    a~_ij = exp(-(t_i - t_j)^2 / (2 s^2)),   a_ij = a~_ij / sum_l a~_il.

    Cada fila suma 1 (se conservan las constantes), pero cerca de los bordes las
    filas tienen sumas distintas antes de normalizar, por lo que en general A != A^T.
    No se usa convolución periódica.
    """
    diferencias = t[:, None] - t[None, :]
    A_tilde = np.exp(-(diferencias ** 2) / (2.0 * s ** 2))
    return A_tilde / A_tilde.sum(axis=1, keepdims=True)


def direccion_ruido(m: int, semilla: int) -> np.ndarray:
    """Vector gaussiano z de longitud m generado con una semilla fija."""
    rng = np.random.default_rng(semilla)
    return rng.standard_normal(m)


def ruido(b_limpio: np.ndarray, eta: float, z: np.ndarray) -> np.ndarray:
    """e = eta * ||A x_true||_2 * z / ||z||_2, de modo que ||e|| / ||A x_true|| = eta exactamente."""
    return eta * np.linalg.norm(b_limpio) * z / np.linalg.norm(z)


def problema(n: int = 128, s: float = 0.03, eta: float = 1e-3, semilla: int = 2026):
    """Construye (t, A, x_true, b_limpio, e, b) para un nivel de ruido eta."""
    t = malla(n)
    A = matriz_suavizado(t, s)
    x_true = senal_verdadera(t)
    b_limpio = A @ x_true
    e = ruido(b_limpio, eta, direccion_ruido(n, semilla))
    return t, A, x_true, b_limpio, e, b_limpio + e
