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


# --------------------------------------------------------------------------
# Principio de discrepancia (no usa x_true)
# --------------------------------------------------------------------------

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


# --------------------------------------------------------------------------
# Referencia oráculo (usa x_true; NO es un método aplicable a datos reales)
# --------------------------------------------------------------------------

def oraculo(reconstrucciones, x_true):
    """Índice y error del parámetro con menor error relativo respecto de x_true."""
    errores = np.array([error_relativo(x, x_true) for x in reconstrucciones])
    idx = int(np.argmin(errores))
    return idx, float(errores[idx])
