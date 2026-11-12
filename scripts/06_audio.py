"""Ejemplo práctico 2: sonido "apagado" (suavizado) y con ruido.

Modelo:  b = A x + e,  con A la convolución circular con un núcleo gaussiano del
mismo ancho que el experimento principal (0.03 * 127 ≈ 3.8 muestras). A 8000
muestras por segundo, ese suavizado deja pasar los graves y borra los agudos,
como escuchar a través de una pared. La SVD de A la da la transformada de
Fourier (ver reconstruccion/circulante.py), así que se puede trabajar con
N = 16 384 muestras sin formar la matriz.

Genera archivos .wav en resultados/audio/ para escuchar el antes y el después.

Uso:
    python scripts/06_audio.py                    # melodía sintética
    python scripts/06_audio.py --audio voz.wav    # agrega una grabación propia
"""

import argparse

import numpy as np
from comun import (C1, C2, C3, LAMBDA_REL, RAIZ, REJILLA, SEMILLA, SUPERFICIE, TAU, TENUE,
                   TINTA, TINTA_2, guardar, guardar_tabla, plt)
from scipy.io import wavfile
from scipy.signal import resample_poly

from reconstruccion import circulante as CIRC
from reconstruccion import datos

FS = 8000                  # muestras por segundo
N_AUDIO = 16384            # 2.048 s
SIGMA_MUESTRAS = 0.03 * 127
ETAS_AUDIO = [1e-3, 1e-2]
COLOR_ETA_AUDIO = {1e-3: C2, 1e-2: C3}
DIR_AUDIO = RAIZ / "resultados" / "audio"
DIR_AUDIO.mkdir(parents=True, exist_ok=True)


def melodia(N=N_AUDIO, fs=FS):
    """Do-Mi-Sol-Do'-Sol-Mi-Do-Sol con 8 armónicos y ataque/caída tipo cuerda pulsada."""
    t = np.arange(N) / fs
    notas = [261.63, 329.63, 392.00, 523.25, 392.00, 329.63, 261.63, 196.00]
    x = np.zeros(N)
    dur = N / len(notas)
    for i, f in enumerate(notas):
        i0 = int(i * dur)
        tt = t[i0:] - t[i0]
        envolvente = np.exp(-4.0 * tt) * (1 - np.exp(-400.0 * tt))
        tono = sum(np.sin(2 * np.pi * h * f * tt) / h ** 1.2 for h in range(1, 9))
        x[i0:] += envolvente * tono
    return 0.8 * x / np.max(np.abs(x))


def cargar_audio(ruta, fs=FS, max_seg=4.0):
    fs_orig, y = wavfile.read(ruta)
    y = y.astype(float)
    if y.ndim > 1:
        y = y.mean(axis=1)
    if fs_orig != fs:
        g = np.gcd(int(fs_orig), fs)
        y = resample_poly(y, fs // g, int(fs_orig) // g)
    y = y[: int(max_seg * fs)]
    return 0.8 * y / np.max(np.abs(y))


def guardar_wav(x, nombre, normalizar=False):
    if normalizar:
        x = 0.3 * x / np.max(np.abs(x))
    wavfile.write(DIR_AUDIO / nombre, FS, (np.clip(x, -1, 1) * 32767).astype(np.int16))


parser = argparse.ArgumentParser()
parser.add_argument("--audio", help="ruta a un .wav propio (opcional)")
args, _ = parser.parse_known_args()

sonidos = {"melodia": melodia()}
if args.audio:
    sonidos["propio"] = cargar_audio(args.audio)

filas = []
res = {}
for nombre, x in sonidos.items():
    N = x.size
    op = CIRC.OperadorCirculante(CIRC.nucleo_gaussiano_circular(N, SIGMA_MUESTRAS))
    lambdas = op.sigma1 ** 2 * LAMBDA_REL
    b_limpio = op.aplicar(x)
    z = datos.direccion_ruido(N, SEMILLA)
    print(f"{nombre}: N = {N} muestras ({N / FS:.2f} s a {FS} Hz); suavizado de {SIGMA_MUESTRAS:.2f} muestras")
    print(f"  sigma_max = {op.sigma.max():.4g}, sigma_min = {op.sigma.min():.2e}")
    guardar_wav(x, f"{nombre}_original.wav")

    for eta in ETAS_AUDIO:
        e = datos.ruido(b_limpio, eta, z)
        b = b_limpio + e
        delta = float(np.linalg.norm(e))
        x_pc = op.pseudoinversa_completa(b)
        lam_d, i_d = CIRC.discrepancia_tikhonov(op, b, lambdas, delta, TAU)
        x_tik = op.tikhonov(b, lam_d)
        err_l = np.array([np.linalg.norm(op.tikhonov(b, lam) - x) / np.linalg.norm(x) for lam in lambdas])
        i_o = int(np.argmin(err_l))
        res[(nombre, eta)] = dict(x=x, b=b, tik=x_tik, op=op, lam=lam_d, i_d=i_d, i_o=i_o, err_l=err_l)

        def fila(metodo, parametro, xr):
            return {"sonido": nombre, "eta": eta, "método": metodo, "parámetro": parametro,
                    "error relativo": float(np.linalg.norm(xr - x) / np.linalg.norm(x)),
                    "SNR (dB)": CIRC.snr_db(xr, x), "residuo / delta": op.residuo(xr, b) / delta}

        filas += [
            fila("Observado (sin reconstruir)", "-", b),
            fila("Pseudoinversa completa", f"{int((op.h_hat != 0).sum())} frecuencias (ĥ_k ≠ 0)", x_pc),
            fila("Tikhonov · discrepancia", f"λ/σ1² = {LAMBDA_REL[i_d]:.2e}", x_tik),
            fila("Tikhonov · oráculo", f"λ/σ1² = {LAMBDA_REL[i_o]:.2e}", op.tikhonov(b, lambdas[i_o])),
        ]
        print(f"  eta = {eta:g}")
        for f in filas[-4:]:
            print(f"    {f['método']:<30} {f['parámetro']:<20} error = {f['error relativo']:.3g}  "
                  f"SNR = {f['SNR (dB)']:8.1f} dB")

        sufijo = f"eta1e{int(np.log10(eta))}"
        guardar_wav(b, f"{nombre}_observado_{sufijo}.wav")
        guardar_wav(x_tik, f"{nombre}_tikhonov_{sufijo}.wav")
        if eta == ETAS_AUDIO[-1]:
            guardar_wav(x_pc, f"{nombre}_pseudoinversa_{sufijo}_VOLUMEN_REDUCIDO.wav", normalizar=True)

guardar_tabla(filas, "06_audio")
print("  audios -> resultados/audio/")

# --------------------------------------------------------------------------------
# Figura: forma de onda, espectro y curva de error (melodía)
# --------------------------------------------------------------------------------
eta = ETAS_AUDIO[-1]
r = res[("melodia", eta)]
x, b, x_tik, op = r["x"], r["b"], r["tik"], r["op"]
fig, axs = plt.subplots(1, 3, figsize=(15, 4.0))
CAJA = dict(frameon=True, facecolor=SUPERFICIE, edgecolor=REJILLA, framealpha=0.95)

i0 = int(3 * N_AUDIO / 8) + 40           # justo después del ataque de la cuarta nota (Do agudo)
ventana = slice(i0, i0 + int(0.02 * FS))
tms = np.arange(ventana.stop - ventana.start) / FS * 1000
axs[0].plot(tms, x[ventana], color=TINTA, lw=2.2, label="original")
axs[0].plot(tms, b[ventana], color=C1, lw=1.4, label="observado (apagado + ruido)")
axs[0].plot(tms, x_tik[ventana], color=C2, lw=1.4, label="Tikhonov (discrepancia)")
axs[0].set_xlabel("tiempo (ms)")
axs[0].set_title("Forma de onda, 20 ms")
axs[0].set_ylim(-0.8, 0.98)
axs[0].set_yticks([-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75])
axs[0].legend(loc="upper center", fontsize=7.5, ncol=3, **CAJA)

f = np.fft.rfftfreq(N_AUDIO, 1 / FS)


def espectro_db(y):
    return 20 * np.log10(np.abs(np.fft.rfft(y)) + 1e-12)


ref = espectro_db(x).max()
axs[1].plot(f, espectro_db(x) - ref, color=TINTA, lw=0.8, label="original")
axs[1].plot(f, espectro_db(b) - ref, color=C1, lw=0.8, alpha=0.8, label="observado")
axs[1].plot(f, espectro_db(x_tik) - ref, color=C2, lw=0.8, alpha=0.9, label="Tikhonov")
sig = op.sigma[: f.size]
axs[1].plot(f, 20 * np.log10(sig ** 2 / (sig ** 2 + r["lam"]) + 1e-300), color=C3, lw=2,
            label=r"factor de filtro $f_k(\lambda)$")
axs[1].plot(f, 20 * np.log10(sig + 1e-300), color=TENUE, lw=1.5, ls="--", label=r"$\sigma_k=|\hat h_k|$")
axs[1].set_xlim(0, 2500)
axs[1].set_ylim(-110, 10)
axs[1].set_xlabel("frecuencia (Hz)")
axs[1].set_ylabel("dB")
axs[1].set_title(rf"Espectro, $\eta=10^{{{int(np.log10(eta))}}}$")
axs[1].legend(loc="upper right", fontsize=7.5, **CAJA)

for e_ in ETAS_AUDIO:
    rr = res[("melodia", e_)]
    col = COLOR_ETA_AUDIO[e_]
    axs[2].loglog(LAMBDA_REL, rr["err_l"], color=col, label=rf"$\eta=10^{{{int(np.log10(e_))}}}$")
    axs[2].plot(LAMBDA_REL[rr["i_o"]], rr["err_l"][rr["i_o"]], "o", color=col, ms=7, mfc="white", mew=2)
    axs[2].plot(LAMBDA_REL[rr["i_d"]], rr["err_l"][rr["i_d"]], "s", color=col, ms=7)
axs[2].plot([], [], "o", color=TENUE, mfc="white", mew=2, label="oráculo")
axs[2].plot([], [], "s", color=TENUE, label="discrepancia")
axs[2].set_ylim(0.1, 10)
axs[2].set_xlabel(r"$\lambda/\sigma_1^2$")
axs[2].set_title("Error relativo del audio reconstruido")
axs[2].legend(loc="upper right", fontsize=8)
for ax in axs:
    ax.title.set_color(TINTA)
    ax.xaxis.label.set_color(TINTA_2)
guardar(fig, "14_audio.png")
