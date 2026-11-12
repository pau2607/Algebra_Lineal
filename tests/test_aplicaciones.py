"""Comprobaciones de los ejemplos prácticos (imagen y audio).

La idea de todas las pruebas es la misma: para tamaños pequeños se FORMA la
matriz grande (Kronecker o circulante) y se compara con las fórmulas que la
evitan. Si coinciden, los resultados de 256 x 256 píxeles y 16 384 muestras
salen de las mismas fórmulas del caso 1D.

Ejecutar con:  python -m pytest -q
"""

import inspect
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as sla
from PIL import Image

from reconstruccion import circulante as CIRC
from reconstruccion import datos
from reconstruccion import metodos as M
from reconstruccion import separable as SEP


def vec(X):
    """Apila las columnas de X (vectorización por columnas)."""
    return X.reshape(-1, order="F")


def desvec(x, forma):
    return x.reshape(forma, order="F")


# --------------------------------------------------------------------------
# Operador separable B = A_v X A_h^T
# --------------------------------------------------------------------------

@pytest.fixture
def separable():
    rng = np.random.default_rng(3)
    # factores con espectros distintos para que los productos sigma_i sigma_j no empaten
    A_v = rng.standard_normal((9, 9)) @ np.diag(np.logspace(0, -5, 9))
    A_h = rng.standard_normal((7, 7)) @ np.diag(np.logspace(0, -3.3, 7))
    op = SEP.OperadorSeparable.de(A_v, A_h)
    K = np.kron(A_h, A_v)
    B = rng.standard_normal((9, 7))
    return op, K, B


def test_aplicar_coincide_con_kronecker(separable):
    op, K, _ = separable
    X = np.random.default_rng(4).standard_normal((9, 7))
    assert np.allclose(vec(op.aplicar(X)), K @ vec(X))


def test_valores_singulares_son_productos(separable):
    op, K, _ = separable
    assert np.allclose(np.sort(op.S.ravel())[::-1], np.linalg.svd(K, compute_uv=False))


@pytest.mark.parametrize("lam", [1e-8, 1e-4, 1e-1])
def test_tikhonov_separable_coincide_con_kronecker(separable, lam):
    op, K, B = separable
    X = op.tikhonov(B, lam)
    x_esp = M.tikhonov_espectral(M.SVD.de(K), vec(B), lam)
    x_aum = M.tikhonov_aumentado(K, vec(B), lam)
    assert np.allclose(vec(X), x_esp, rtol=1e-6, atol=1e-9)
    assert np.allclose(vec(X), x_aum, rtol=1e-6, atol=1e-9)


@pytest.mark.parametrize("k", [1, 5, 20, 63])
def test_svd_truncada_separable_coincide_con_kronecker(separable, k):
    op, K, B = separable
    assert np.allclose(vec(op.svd_truncada(B, k)), M.svd_truncada(M.SVD.de(K), vec(B), k), rtol=1e-6, atol=1e-9)


def test_pseudoinversa_separable_coincide_con_kronecker(separable):
    op, K, B = separable
    assert np.allclose(vec(op.pseudoinversa_completa(B)), np.linalg.solve(K, vec(B)), rtol=1e-6)


def test_residuo_espectral_igual_al_directo(separable):
    op, _, B = separable
    for lam in [1e-6, 1e-3, 1.0]:
        assert op.residuo_tikhonov_espectral(B, lam) == pytest.approx(op.residuo(op.tikhonov(B, lam), B), rel=1e-9)


def test_audio_por_bloques_es_el_caso_1d_columna_a_columna():
    """Con A_h = I, Tikhonov separable = Tikhonov 1D aplicado a cada columna."""
    t = datos.malla(32)
    A = datos.matriz_suavizado(t, 0.03)
    B = np.random.default_rng(5).standard_normal((32, 6))
    op = SEP.OperadorSeparable.de(A, np.eye(6))
    svd = M.SVD.de(A)
    X = op.tikhonov(B, 1e-4)
    for j in range(6):
        assert np.allclose(X[:, j], M.tikhonov_espectral(svd, B[:, j], 1e-4))


def test_ruido_matricial_nivel_exacto():
    B = np.random.default_rng(6).random((20, 30))
    E = SEP.ruido_matricial(B, 1e-2, 2026)
    assert np.linalg.norm(E) / np.linalg.norm(B) == pytest.approx(1e-2, rel=1e-12)


def test_discrepancia_separable_no_recibe_la_verdad():
    for f in (SEP.discrepancia_tikhonov, SEP.discrepancia_tsvd):
        assert not {"X", "X_true", "x_true"} & set(inspect.signature(f).parameters)


def test_imagen_regularizada_mejora_a_la_observada():
    """Gato reducido a 64 x 64: Tikhonov por discrepancia mejora a la imagen observada;
    la pseudoinversa completa es inservible."""
    ruta = Path(__file__).resolve().parents[1] / "recursos" / "gato_256.png"
    X = np.asarray(Image.open(ruta).resize((64, 64), Image.LANCZOS), dtype=float) / 255
    A = datos.matriz_suavizado(datos.malla(64), 0.03)
    op = SEP.OperadorSeparable.de(A, A)
    B_limpio = op.aplicar(X)
    E = SEP.ruido_matricial(B_limpio, 1e-2, 2026)
    B = B_limpio + E
    lambdas = op.sigma1 ** 2 * np.logspace(-10, 0, 61)
    lam, _ = SEP.discrepancia_tikhonov(op, B, lambdas, np.linalg.norm(E))
    error = lambda Y: np.linalg.norm(Y - X) / np.linalg.norm(X)  # noqa: E731
    assert op.residuo(op.tikhonov(B, lam), B) <= 1.1 * np.linalg.norm(E)
    assert error(op.tikhonov(B, lam)) < error(B)
    assert error(op.pseudoinversa_completa(B)) > 1e3


# --------------------------------------------------------------------------
# Operador circulante (audio)
# --------------------------------------------------------------------------

@pytest.fixture
def circular():
    h = CIRC.nucleo_gaussiano_circular(48, 2.5)
    return CIRC.OperadorCirculante(h), sla.circulant(h)


def test_nucleo_suma_uno_y_es_simetrico():
    h = CIRC.nucleo_gaussiano_circular(100, 3.8)
    assert h.sum() == pytest.approx(1.0)
    assert np.allclose(h[1:], h[1:][::-1])


def test_circulante_aplicar_es_producto_matricial(circular):
    op, C = circular
    x = np.random.default_rng(7).standard_normal(48)
    assert np.allclose(op.aplicar(x), C @ x)


def test_valores_singulares_circulante_son_fourier(circular):
    op, C = circular
    assert np.allclose(np.sort(op.sigma)[::-1], np.linalg.svd(C, compute_uv=False), atol=1e-14)


@pytest.mark.parametrize("lam", [1e-8, 1e-4, 1e-1])
def test_tikhonov_fourier_coincide_con_svd(circular, lam):
    op, C = circular
    b = np.random.default_rng(8).standard_normal(48)
    x = op.tikhonov(b, lam)
    assert np.allclose(x, M.tikhonov_espectral(M.SVD.de(C), b, lam), rtol=1e-6, atol=1e-9)
    assert np.allclose(x, M.tikhonov_aumentado(C, b, lam), rtol=1e-6, atol=1e-9)


def test_discrepancia_circulante_cumple_su_definicion(circular):
    op, C = circular
    rng = np.random.default_rng(9)
    x = rng.standard_normal(48)
    e = 1e-3 * rng.standard_normal(48)
    b = op.aplicar(x) + e
    lambdas = np.logspace(-10, 0, 61)
    lam, i = CIRC.discrepancia_tikhonov(op, b, lambdas, np.linalg.norm(e))
    assert op.residuo(op.tikhonov(b, lam), b) <= 1.1 * np.linalg.norm(e)
    if i + 1 < len(lambdas):
        assert op.residuo(op.tikhonov(b, lambdas[i + 1]), b) > 1.1 * np.linalg.norm(e)
