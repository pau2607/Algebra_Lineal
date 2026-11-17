"""Selección del parámetro de regularización.

Se separan dos cosas que no deben mezclarse:

* Criterio de selección (discrepancia): usa A, b y una cota de ruido delta.
  Las funciones de discrepancia NO reciben x_true.
* Evaluación oráculo: usa x_true para hallar el mejor parámetro de la malla.
  Solo existe en un experimento sintético y sirve como referencia.
"""

import numpy as np

from .metodos import SVD, svd_truncada, tikhonov_espectral
from .metricas import error_relativo, residuo



# Principio de discrepancia (no usa x_true)


def discrepancia_tikhonov(A, b, svd: SVD, lambdas, delta: float, tau: float = 1.1):
    """Mayor lambda de la malla cuyo residuo no excede tau * delta.

    Devuelve (lambda_elegido, indice) o (None, None) si ningún lambda de la malla
    alcanza el umbral (el criterio no es viable en esa malla).
    """
    lambdas = np.asarray(lambdas)
    residuos = np.array([residuo(A, tikhonov_espectral(svd, b, lam), b) for lam in lambdas])
    viables = np.flatnonzero(residuos <= tau * delta)
    if viables.size == 0:
        return None, None
    idx = viables[np.argmax(lambdas[viables])]
    return float(lambdas[idx]), int(idx)


def discrepancia_tsvd(A, b, svd: SVD, ks, delta: float, tau: float = 1.1):
    """Menor k (más regularización) cuyo residuo no excede tau * delta.

    Devuelve (k, indice) o (None, None) si no es viable.
    """
    ks = np.asarray(ks)
    residuos = np.array([residuo(A, svd_truncada(svd, b, int(k)), b) for k in ks])
    viables = np.flatnonzero(residuos <= tau * delta)
    if viables.size == 0:
        return None, None
    idx = viables[np.argmin(ks[viables])]
    return int(ks[idx]), int(idx)


# Referencia oráculo (usa x_true; NO es un método aplicable a datos reales)


def oraculo(reconstrucciones, x_true):
    """Índice y error del parámetro con menor error relativo respecto de x_true."""
    errores = np.array([error_relativo(x, x_true) for x in reconstrucciones])
    idx = int(np.argmin(errores))
    return idx, float(errores[idx])


# --------------------------------------------------------------------------
# Estimación del nivel de ruido a partir de los datos (no usa e ni x_true)
# --------------------------------------------------------------------------

def estimar_delta(svd: SVD, b, fraccion_cola: float = 0.5) -> float:
    """Estima ||e|| con la cola de los coeficientes de Picard.

    Supuesto: en las últimas direcciones (i > (1 - fraccion_cola) * n) los sigma_i son tan
    pequeños que u_i^T b = u_i^T e, es decir, la cola es solo ruido. Si e es isotrópico,
    E[(u_i^T e)^2] = ||e||^2 / m, así que

        delta_hat = sqrt( m * media_{i en la cola} (u_i^T b)^2 ).

    Para ruido de norma fija y dirección uniforme en la esfera, delta_hat / ||e|| depende solo
    de la dirección de e y su desviación relativa es ~ sqrt((m - q) / (2 q (m + 2))) (≈ 6 % con
    m = 128, q = 64). El valor 1/sqrt(2q) ≈ 9 % corresponde a ruido gaussiano comparado con su
    sigma*sqrt(m), no con la norma realizada.
    """
    c = svd.coeficientes(b)
    m = svd.U.shape[0]          # número de filas de A (dimensión de e), no min(m, n)
    i0 = int(round((1.0 - fraccion_cola) * m))
    cola = c[i0:]
    return float(np.sqrt(m * np.mean(cola ** 2)))


def discrepancia_tikhonov_refinada(A, b, svd: SVD, lambdas, delta: float, tau: float = 1.1, iteraciones: int = 40):
    """Como `discrepancia_tikhonov`, pero refina por bisección entre el último lambda viable y el
    siguiente de la malla, usando que el residuo es no decreciente en lambda (desarrollo §6.2).
    Devuelve lambda con ||A x_lambda - b|| = tau * delta hasta la precisión de la bisección, o None.
    """
    lambdas = np.asarray(lambdas)
    lam, idx = discrepancia_tikhonov(A, b, svd, lambdas, delta, tau)
    if lam is None or idx + 1 >= lambdas.size:
        return lam
    lo, hi = lambdas[idx], lambdas[idx + 1]          # residuo(lo) <= tau*delta < residuo(hi)
    for _ in range(iteraciones):
        mid = np.sqrt(lo * hi)                        # bisección en escala logarítmica
        if residuo(A, tikhonov_espectral(svd, b, mid), b) <= tau * delta:
            lo = mid
        else:
            hi = mid
    return float(lo)
