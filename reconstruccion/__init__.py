"""Reconstrucción de señales con ruido: pseudoinversa, SVD truncada y Tikhonov.

Módulos:
    datos      -- señal verdadera, matriz de suavizado y ruido controlado
    metodos    -- SVD, pseudoinversa, SVD truncada y Tikhonov (fórmula espectral y sistema aumentado)
    seleccion  -- principio de discrepancia (sin x_true) y referencia oráculo (con x_true)
    metricas   -- error relativo, residuo y norma
"""

from . import datos, metodos, seleccion, metricas  # noqa: F401
