"""Comprobaciones automáticas de las fórmulas y de la construcción del problema.

Ejecutar con:  python -m pytest -q
"""

import inspect

import numpy as np
import pytest

from reconstruccion import datos
from reconstruccion import metodos as M
from reconstruccion import seleccion as SEL
from reconstruccion.metricas import error_relativo, residuo


# --------------------------------------------------------------------------
# Ejemplo diagonal (sección 3.1, 4 y 5 del artículo)
# --------------------------------------------------------------------------

@pytest.fixture
def diagonal():
    A = np.diag([1.0, 1e-3])
    x_true = np.array([1.0, 1.0])
    b = np.array([1.0, 0.011])
    return A, x_true, b, M.SVD.de(A)


def test_diagonal_pseudoinversa(diagonal):
    A, x_true, b, svd = diagonal
    x = M.pseudoinversa_completa(svd, b)
    assert np.allclose(x, [1.0, 11.0])
    assert residuo(A, x, b) == pytest.approx(0.0, abs=1e-14)       # residuo cero ...
    assert np.linalg.norm(x - x_true) == pytest.approx(10.0)       # ... y error 10


def test_diagonal_svd_truncada(diagonal):
    A, x_true, b, svd = diagonal
    x = M.svd_truncada(svd, b, 1)
    assert np.allclose(x, [1.0, 0.0])
    assert np.linalg.norm(x - x_true) == pytest.approx(1.0)
    assert residuo(A, x, b) == pytest.approx(0.011)


def test_diagonal_tikhonov(diagonal):
    A, x_true, b, svd = diagonal
    x = M.tikhonov_espectral(svd, b, 1e-4)
    assert x == pytest.approx([0.99990, 0.10891], abs=1e-5)
    assert np.linalg.norm(x - x_true) == pytest.approx(0.8911, abs=1e-4)


# --------------------------------------------------------------------------
# Construcción del problema (sección 8)
# --------------------------------------------------------------------------

def test_matriz_conserva_constantes_y_no_es_simetrica():
    t = datos.malla(128)
    A = datos.matriz_suavizado(t, 0.03)
    assert np.allclose(A @ np.ones(128), np.ones(128))
    assert not np.allclose(A, A.T)
    assert np.all(A >= 0)


def test_ruido_tiene_nivel_relativo_exacto():
    t, A, x_true, b_limpio, e, b = datos.problema(eta=1e-3, semilla=7)
    assert np.linalg.norm(e) / np.linalg.norm(A @ x_true) == pytest.approx(1e-3, rel=1e-12)
    assert np.allclose(b, b_limpio + e)


def test_semilla_reproducible():
    z1 = datos.direccion_ruido(128, 2026)
    z2 = datos.direccion_ruido(128, 2026)
    assert np.array_equal(z1, z2)


def test_senal_verdadera_valores_puntuales():
    t = np.array([0.3, 0.65, 0.9])
    x = datos.senal_verdadera(t)
    assert x[0] == pytest.approx(1.0 + 0.0)               # pico de la gaussiana, fuera del pulso
    assert x[1] == pytest.approx(np.exp(-80 * 0.35 ** 2) + 0.6)
    assert x[2] == pytest.approx(np.exp(-80 * 0.6 ** 2))


# --------------------------------------------------------------------------
# Tikhonov: fórmula espectral, sistema aumentado y propiedades (secciones 5 y 6)
# --------------------------------------------------------------------------

@pytest.fixture
def aleatorio():
    rng = np.random.default_rng(0)
    A = rng.standard_normal((30, 20)) @ np.diag(np.logspace(0, -6, 20))
    b = rng.standard_normal(30)
    return A, b, M.SVD.de(A)


@pytest.mark.parametrize("lam", [1e-8, 1e-4, 1e-1, 10.0])
def test_tikhonov_tres_formas_coinciden(aleatorio, lam):
    A, b, svd = aleatorio
    xe = M.tikhonov_espectral(svd, b, lam)
    xa = M.tikhonov_aumentado(A, b, lam)
    assert np.allclose(xe, xa, rtol=1e-6, atol=1e-8)
    if lam >= 1e-4:   # las ecuaciones normales solo se comparan con lambda moderado
        assert np.allclose(xe, M.tikhonov_normales(A, b, lam), rtol=1e-6)


def test_tikhonov_es_punto_estacionario(aleatorio):
    A, b, svd = aleatorio
    lam = 1e-3
    x = M.tikhonov_espectral(svd, b, lam)
    gradiente = 2 * A.T @ (A @ x - b) + 2 * lam * x
    assert np.linalg.norm(gradiente) < 1e-10


def test_factores_filtro_entre_0_y_1_y_decrecientes(aleatorio):
    _, _, svd = aleatorio
    f = M.factores_filtro(svd.sigma, 1e-3)
    assert np.all((f > 0) & (f < 1))
    assert np.all(np.diff(f) <= 0)


def test_cota_de_estabilidad(aleatorio):
    A, b, svd = aleatorio
    rng = np.random.default_rng(1)
    for lam in [1e-6, 1e-3, 1.0]:
        for _ in range(20):
            e = rng.standard_normal(b.size)
            dif = np.linalg.norm(M.tikhonov_espectral(svd, b + e, lam) - M.tikhonov_espectral(svd, b, lam))
            assert dif <= np.linalg.norm(e) / (2 * np.sqrt(lam)) * (1 + 1e-12)


def test_limites_de_lambda():
    rng = np.random.default_rng(2)
    A = rng.standard_normal((15, 10))   # bien condicionada
    b = rng.standard_normal(15)
    svd = M.SVD.de(A)
    assert np.allclose(M.tikhonov_espectral(svd, b, 1e-12), np.linalg.pinv(A) @ b, atol=1e-8)
    assert np.linalg.norm(M.tikhonov_espectral(svd, b, 1e12)) < 1e-10


def test_svd_truncada_completa_es_pseudoinversa(aleatorio):
    A, b, svd = aleatorio
    assert np.allclose(M.svd_truncada(svd, b, len(svd.sigma)), M.pseudoinversa_completa(svd, b))


def test_residuo_tsvd_no_crece_con_k(aleatorio):
    A, b, svd = aleatorio
    res = [residuo(A, M.svd_truncada(svd, b, k), b) for k in range(0, 21)]
    assert np.all(np.diff(res) <= 1e-12)


def test_residuo_tikhonov_crece_con_lambda(aleatorio):
    A, b, svd = aleatorio
    res = [residuo(A, M.tikhonov_espectral(svd, b, lam), b) for lam in np.logspace(-10, 2, 30)]
    assert np.all(np.diff(res) >= -1e-12)


# --------------------------------------------------------------------------
# Pseudoinversa numérica y selección del parámetro
# --------------------------------------------------------------------------

def test_pseudoinversa_numerica_informa_rango_y_umbral():
    t = datos.malla(128)
    A = datos.matriz_suavizado(t, 0.03)
    _, rango, umbral = M.pseudoinversa_numerica(A, np.ones(128))
    sigma = np.linalg.svd(A, compute_uv=False)
    assert rango == int((sigma > umbral).sum())
    assert rango < 128     # la pinv de SciPy ya descarta direcciones


def test_discrepancia_no_recibe_la_verdad():
    for f in (SEL.discrepancia_tikhonov, SEL.discrepancia_tsvd):
        assert "x_true" not in inspect.signature(f).parameters


def test_discrepancia_cumple_su_definicion():
    t, A, x_true, _, e, b = datos.problema(eta=1e-3)
    svd = M.SVD.de(A)
    lambdas = svd.sigma[0] ** 2 * np.logspace(-10, 0, 61)
    delta, tau = np.linalg.norm(e), 1.1
    lam, i = SEL.discrepancia_tikhonov(A, b, svd, lambdas, delta, tau)
    assert residuo(A, M.tikhonov_espectral(svd, b, lam), b) <= tau * delta
    if i + 1 < len(lambdas):  # el siguiente lambda de la malla ya no es viable
        assert residuo(A, M.tikhonov_espectral(svd, b, lambdas[i + 1]), b) > tau * delta


def test_discrepancia_no_viable_devuelve_none():
    t, A, x_true, _, e, b = datos.problema(eta=1e-3)
    svd = M.SVD.de(A)
    lambdas = svd.sigma[0] ** 2 * np.logspace(-10, 0, 61)
    assert SEL.discrepancia_tikhonov(A, b, svd, lambdas, delta=1e-12) == (None, None)


def test_regularizacion_mejora_a_la_pseudoinversa():
    t, A, x_true, _, e, b = datos.problema(eta=1e-3)
    svd = M.SVD.de(A)
    x_pinv, _, _ = M.pseudoinversa_numerica(A, b)
    lambdas = svd.sigma[0] ** 2 * np.logspace(-10, 0, 61)
    lam, _ = SEL.discrepancia_tikhonov(A, b, svd, lambdas, np.linalg.norm(e))
    assert error_relativo(M.tikhonov_espectral(svd, b, lam), x_true) < 0.2
    assert error_relativo(x_pinv, x_true) > 1e6
