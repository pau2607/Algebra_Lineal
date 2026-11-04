"""Métricas de evaluación."""

import numpy as np


def error_relativo(x: np.ndarray, x_true: np.ndarray) -> float:
    """||x - x_true|| / ||x_true||. Requiere la verdad: solo para evaluación posterior."""
    return float(np.linalg.norm(x - x_true) / np.linalg.norm(x_true))


def residuo(A: np.ndarray, x: np.ndarray, b: np.ndarray) -> float:
    """||A x - b||_2. Calculable sin conocer x_true."""
    return float(np.linalg.norm(A @ x - b))


def norma(x: np.ndarray) -> float:
    """||x||_2."""
    return float(np.linalg.norm(x))
