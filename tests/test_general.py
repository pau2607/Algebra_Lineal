"""Pruebas de la ampliación: Tikhonov generalizado con penalización lambda ||L x||^2."""

import numpy as np
import pytest

from reconstruccion import datos
from reconstruccion import general as G
from reconstruccion import metodos as M
from reconstruccion.metricas import residuo

N = 128


@pytest.fixture(scope="module")
def problema():
    t, A, x_true, b_limpio, e, b = datos.problema(n=N, eta=1e-3)
    return t, A, x_true, b_limpio, e, b


def test_matrices_de_diferencias():
    L1, L2 = G.primeras_diferencias(6), G.segundas_diferencias(6)
    x = np.array([1.0, 4.0, 9.0, 16.0, 25.0, 36.0])          # cuadrados: diferencias 3,5,7,9,11; segundas 2,2,2,2
    assert L1.shape == (5, 6) and L2.shape == (4, 6)
    assert np.allclose(L1 @ x, [3, 5, 7, 9, 11])
    assert np.allclose(L2 @ x, [2, 2, 2, 2])
    assert np.allclose(L1 @ np.ones(6), 0)                  # ker L1 contiene las constantes
    assert np.allclose(L2 @ np.arange(6), 0)                # ker L2 contiene las rectas


def test_nucleos_disjuntos(problema):
    _, A, *_ = problema
    assert G.nucleos_disjuntos(A, G.primeras_diferencias(N))      # A 1 = 1 ≠ 0
    assert G.nucleos_disjuntos(A, G.segundas_diferencias(N))
    # Si A anula las constantes (A centrada), la condición falla con L1
    A_centrada = A - A.mean(axis=1, keepdims=True)
    assert np.allclose(A_centrada @ np.ones(N), 0, atol=1e-12)
    assert not G.nucleos_disjuntos(A_centrada, G.primeras_diferencias(N))


def test_gsvd_reconstruye_A_y_L(problema):
    _, A, *_ = problema
    for L in (np.eye(N), G.primeras_diferencias(N), G.segundas_diferencias(N)):
        g = G.GSVD.de(A, L)
        Xinv = np.linalg.inv(g.X)
        assert np.allclose(g.U @ np.diag(g.c) @ Xinv, A, atol=1e-10)
        # L = V diag(s) X^{-1}  =>  L X = V diag(s)  =>  las columnas de L X tienen norma s_i
        assert np.allclose(np.linalg.norm(L @ g.X, axis=0), g.s, atol=1e-10)
        assert np.allclose(g.c ** 2 + g.s ** 2, 1.0)


@pytest.mark.parametrize("lam", [1e-8, 1e-3, 1.0, 1e3])
def test_espectral_coincide_con_aumentado(problema, lam):
    _, A, _, _, _, b = problema
    for L in (G.primeras_diferencias(N), G.segundas_diferencias(N)):
        g = G.GSVD.de(A, L)
        xe = G.tikhonov_general_espectral(g, b, lam)
        xa = G.tikhonov_general_aumentado(A, L, b, lam)
        assert np.allclose(xe, xa, rtol=1e-6, atol=1e-8)


def test_con_L_identidad_es_tikhonov_clasico(problema):
    _, A, _, _, _, b = problema
    g = G.GSVD.de(A, np.eye(N))
    svd = M.SVD.de(A)
    for lam in (1e-6, 1e-3, 1e-1):
        assert np.allclose(G.tikhonov_general_espectral(g, b, lam), M.tikhonov_espectral(svd, b, lam), rtol=1e-6, atol=1e-9)


def test_punto_estacionario(problema):
    _, A, _, _, _, b = problema
    L = G.primeras_diferencias(N)
    lam = 1e-2
    x = G.tikhonov_general_aumentado(A, L, b, lam)
    gradiente = 2 * A.T @ (A @ x - b) + 2 * lam * L.T @ (L @ x)
    assert np.linalg.norm(gradiente) < 1e-9


def test_limite_lambda_grande_es_la_media(problema):
    """Con L1 y lambda -> inf, x tiende a la constante c que minimiza ||A c1 - b||; como A1 = 1, c = media(b)."""
    _, A, _, _, _, b = problema
    g = G.GSVD.de(A, G.primeras_diferencias(N))
    x = G.tikhonov_general_espectral(g, b, 1e12)
    assert np.allclose(x, b.mean(), atol=1e-6)


def test_factores_filtro_generalizados(problema):
    _, A, *_ = problema
    g = G.GSVD.de(A, G.primeras_diferencias(N))
    f = G.factores_filtro_general(g, 1e-3)
    assert np.all((f > 0) & (f <= 1))
    assert np.isclose(f[np.argmax(g.c)], 1.0)                   # dirección del núcleo de L: s = 0, f = 1
    assert np.isinf(g.gamma).sum() == 1                          # un solo gamma infinito (dim ker L1 = 1)


def test_residuo_no_decreciente_en_lambda(problema):
    _, A, _, _, _, b = problema
    L = G.primeras_diferencias(N)
    g = G.GSVD.de(A, L)
    res = [residuo(A, G.tikhonov_general_espectral(g, b, lam), b) for lam in np.logspace(-12, 4, 40)]
    assert np.all(np.diff(res) >= -1e-10)


def test_discrepancia_general_no_usa_la_verdad_y_cumple_su_definicion(problema):
    import inspect
    assert "x_true" not in inspect.signature(G.discrepancia_general).parameters
    _, A, _, _, e, b = problema
    L = G.segundas_diferencias(N)
    lambdas = np.logspace(-14, 4, 109)
    delta, tau = np.linalg.norm(e), 1.1
    lam, i = G.discrepancia_general(A, L, b, lambdas, delta, tau)
    assert residuo(A, G.tikhonov_general_aumentado(A, L, b, lam), b) <= tau * delta
    assert residuo(A, G.tikhonov_general_aumentado(A, L, b, lambdas[i + 1]), b) > tau * delta


def test_penalizacion_de_diferencias_mejora_la_senal_suave():
    t = datos.malla(N)
    A = datos.matriz_suavizado(t, 0.03)
    x = datos.senal_suave(t)
    b_limpio = A @ x
    e = datos.ruido(b_limpio, 1e-3, datos.direccion_ruido(N, 2026))
    b = b_limpio + e
    lambdas = np.logspace(-14, 4, 109)
    errores = {}
    for nombre, L in (("I", np.eye(N)), ("L2", G.segundas_diferencias(N))):
        g = G.GSVD.de(A, L)
        lam, _ = G.discrepancia_general(A, L, b, lambdas, np.linalg.norm(e), gsvd=g)
        xr = G.tikhonov_general_espectral(g, b, lam)
        errores[nombre] = np.linalg.norm(xr - x) / np.linalg.norm(x)
    assert errores["L2"] < errores["I"] / 5
