"""Verificacion de la secuencia 3D tsSOS-QISS bSSFP a 0.55T.

Todo lo que se comprueba aqui se MIDE sobre la secuencia ya ensamblada (bloques y
formas de onda), no sobre las variables con las que se construyo.  La version
anterior declaraba "balanceado" sumando las areas de los eventos, que es una
tautologia, y por eso no vio que faltaba la prefase antes del primer alpha.

  1. Periodicidad: intervalos RF-RF identicos, alpha/2 a TR/2, TE = TR/2.
  2. Momento de gradiente nulo en el CENTRO de cada RF del tren (x, y, z).
  3. Codificacion: kz en cada eco = indice * delta_kz; k in-plane = 0 en el eco.
  4. QI medido: centro de la inversion -> eco del TR que muestrea kz = 0.
  5. Periodo del disparo = TR_QISS.
  6. Limites de gradiente por eje y magnitud vectorial (PNS).
  7. Eficiencia de inversion en musculo frente a la constante usada para el nulo.
  8. Reparto angular de las vistas.
  9. check_timing() y test_report() de PyPulseq.

Sale con codigo 1 si alguna comprobacion falla.
"""

import sys

import numpy as np
import pypulseq as pp

import qiss3d_tssos_bssfp_055T as q
from bloch_tools import bloch_rf

fallos = []


def check(nombre, ok, detalle=""):
    print(f"  [{'OK' if ok else 'FALLA'}] {nombre}" + (f"  {detalle}" if detalle else ""))
    if not ok:
        fallos.append(nombre)


N_VIEWS_TEST = 3
print(f"Construyendo secuencia de prueba (1 slab, {N_VIEWS_TEST} vistas)...")
seq, core = q.build_sequence(n_slabs=1, n_views=N_VIEWS_TEST, verbose=False)

# ----------------------------------------------------------------- recorrido de bloques
t = 0.0
exc_centers, half_centers, inv_centers, adc_centers, shot_starts = [], [], [], [], []
for i in range(1, len(seq.block_events) + 1):
    b = seq.get_block(i)
    rf = getattr(b, "rf", None)
    if rf is not None:
        dur = rf.shape_dur
        if abs(dur - q.INV_DURATION) < 1e-9:
            inv_centers.append(t + rf.delay + q.INV_DURATION / 2)
            shot_starts.append(t)
        elif abs(dur - q.EXC_DURATION) < 1e-9:
            c = t + rf.delay + pp.calc_rf_center(rf)[0]
            # el alpha/2 es el de mitad de amplitud
            if np.abs(rf.signal).max() < 0.75 * np.abs(core.rf.signal).max():
                half_centers.append(c)
            else:
                exc_centers.append(c)
    adc = getattr(b, "adc", None)
    if adc is not None:
        adc_centers.append(t + adc.delay + adc.num_samples * adc.dwell / 2)
    t += seq.block_durations[i]
t_total = t

N = q.N_KZ_ACQUIRED
exc = np.array(exc_centers).reshape(N_VIEWS_TEST, N)
ech = np.array(adc_centers).reshape(N_VIEWS_TEST, N)
half = np.array(half_centers)
inv = np.array(inv_centers)

# ----------------------------------------------------------------- integral exacta de G
waveforms = seq.waveforms_and_times()[0]


def moment(axis, t_query):
    """Integral exacta de la forma de onda lineal a tramos hasta t_query [1/m]."""
    tt, g = np.asarray(waveforms[axis])
    cum = np.concatenate([[0.0], np.cumsum(np.diff(tt) * (g[:-1] + g[1:]) / 2)])
    tq = np.atleast_1d(t_query)
    idx = np.clip(np.searchsorted(tt, tq, side="right") - 1, 0, len(tt) - 2)
    dt = np.clip(tq - tt[idx], 0, None)
    seg = np.diff(tt)[idx]
    slope = np.where(seg > 0, (g[idx + 1] - g[idx]) / np.where(seg > 0, seg, 1), 0.0)
    inside = tq <= tt[-1]
    val = cum[idx] + g[idx] * dt + 0.5 * slope * dt ** 2
    return np.where(inside, val, cum[-1])


US = 1e6
print("\n" + "=" * 74 + "\n1. PERIODICIDAD DEL TREN\n" + "=" * 74)
d_rf = np.diff(exc, axis=1)
check("intervalos RF-RF identicos", np.ptp(d_rf) < 1e-9,
      f"{d_rf.min()*US:.1f}..{d_rf.max()*US:.1f} us (TR = {core.tr*US:.1f} us)")
check("intervalo RF-RF = TR declarado", abs(d_rf.mean() - q.TR_BSSFP) < 1e-9)
d_half = exc[:, 0] - half
check("alpha/2 a TR/2 del primer alpha", np.allclose(d_half, core.tr / 2, atol=1e-9),
      f"{d_half[0]*US:.1f} us")
te = ech - exc
check("TE = TR/2 en todos los TR", np.allclose(te, core.tr / 2, atol=1e-9),
      f"TE {te.min()*US:.1f}..{te.max()*US:.1f} us")

print("\n" + "=" * 74 + "\n2. MOMENTO NULO EN EL CENTRO DE CADA RF (condicion bSSFP)\n" + "=" * 74)
dk_xy, dk_z = q.DELTA_KZ * 0 + 1 / q.FOV_XY, q.DELTA_KZ
for ax, nombre, dk in [(0, "x", dk_xy), (1, "y", dk_xy), (2, "z", dk_z)]:
    k_half = moment(ax, half)[:, None]
    res = np.abs(moment(ax, exc.ravel()).reshape(exc.shape) - k_half).max()
    check(f"eje {nombre}: momento alpha/2 -> cada alpha", res < 1e-3 * dk,
          f"residuo max {res:.2e} 1/m ({res/dk*100:.4f} % de delta_k)")

print("\n" + "=" * 74 + "\n3. CODIFICACION\n" + "=" * 74)
kz = (moment(2, ech.ravel()).reshape(ech.shape) - moment(2, exc.ravel()).reshape(exc.shape))
kz_idx = kz / q.DELTA_KZ
check("kz en el eco = indice * delta_kz",
      np.allclose(kz_idx, q.KZ_INDICES[None, :], atol=1e-3),
      f"{kz_idx[0].round(2).min():+.0f}..{kz_idx[0].round(2).max():+.0f}")
kxy = np.hypot(moment(0, ech.ravel()) - moment(0, exc.ravel()),
               moment(1, ech.ravel()) - moment(1, exc.ravel())).max()
check("k in-plane = 0 en el centro del ADC", kxy < 1e-2 * dk_xy, f"{kxy:.2e} 1/m")
check("kz = 0 en el TR que el codigo cree", abs(kz_idx[0, q.IDX_KZ0]) < 1e-3,
      f"TR #{q.IDX_KZ0+1} de {N}")
kz_max = np.abs(q.KZ_INDICES).max() * q.DELTA_KZ
check("resolucion en z", abs(1 / (2 * kz_max) - q.PARTITION_THICKNESS) < 1e-9,
      f"{1/(2*kz_max)*1e3:.3f} mm")

print("\n" + "=" * 74 + "\n4-5. INTERVALO QUIESCENTE Y PERIODO\n" + "=" * 74)
qi_med = ech[:, q.IDX_KZ0] - inv
# El inicio del tren se cuantiza al raster de 10 us: se admite un paso de diferencia.
check("QI medido = QI de diseno", np.allclose(qi_med, q.QI, atol=10e-6),
      f"{qi_med[0]*1e3:.3f} ms vs {q.QI*1e3:.3f} ms")
periodo = np.diff(shot_starts)
check("periodo del disparo = TR_QISS", np.allclose(periodo, q.TR_QISS, atol=1e-9),
      f"{periodo[0]*1e3:.3f} ms")

print("\n" + "=" * 74 + "\n6. LIMITES DE GRADIENTE\n" + "=" * 74)
g_axis = np.array([np.abs(np.asarray(waveforms[a])[1]).max() for a in range(3)])
g_axis_mt = g_axis / q.GAMMA_BAR * 1e3
lim = q.MAX_GRAD_MT_M * q.GRAD_DERATE
check("amplitud por eje dentro del limite derateado", (g_axis_mt <= lim + 1e-6).all(),
      f"{g_axis_mt[0]:.2f} / {g_axis_mt[1]:.2f} / {g_axis_mt[2]:.2f} mT/m (limite {lim:.1f})")
t_common = np.arange(0, min(t_total, q.TR_QISS), 2e-6)
g_i = np.vstack([np.interp(t_common, *np.asarray(waveforms[a])) for a in range(3)])
g_vec = np.linalg.norm(g_i, axis=0).max() / q.GAMMA_BAR * 1e3
s_vec = np.linalg.norm(np.gradient(g_i, t_common, axis=1), axis=0).max() / q.GAMMA_BAR
print(f"  [info] magnitud vectorial: {g_vec:.1f} mT/m, slew vectorial: {s_vec:.1f} T/m/s")
print("         El vector supera el nameplate por eje, normal en 3D.  El PNS si depende")
print("         del vector y PyPulseq no lo modela: pasar el .seq por el chequeo del")
print("         fabricante (el bore de 80 cm juega en contra).")

print("\n" + "=" * 74 + "\n7. EFICIENCIA DE INVERSION EN MUSCULO\n" + "=" * 74)
rf_inv, gz_inv = q.make_adiabatic_inversion(q.system, q.INV_SLAB_THICKNESS, 0.0)
f_med = -bloch_rf(rf_inv, gz_inv.amplitude, [0.0], q.T1_MUSCLE, q.T2_MUSCLE)[2][0]
check("coincide con INV_EFFICIENCY_MUSCLE", abs(f_med - q.INV_EFFICIENCY_MUSCLE) < 0.01,
      f"medida {f_med:.3f}, constante {q.INV_EFFICIENCY_MUSCLE:.3f}")

print("\n" + "=" * 74 + "\n8. REPARTO ANGULAR DE LAS VISTAS\n" + "=" * 74)


def gap_stats(angles_deg):
    a = np.sort(np.asarray(angles_deg) % 180.0)
    gaps = np.append(np.diff(a), 180 - a[-1] + a[0])
    return gaps.max() / (180 / len(a))


n = q.N_VIEWS
for orden in ["paper", "golden", "uniform"]:
    ang = np.degrees([q.view_angle(i, n, orden) for i in range(n)])
    marca = "  <- en uso" if orden == q.VIEW_ORDER else ""
    print(f"  {orden:8s}: hueco maximo = {gap_stats(ang):.2f} x el ideal{marca}")
print("  (uniforme reparte perfecto pero deja streaks coherentes; el aureo los vuelve")
print("   incoherentes, que es lo que aprovecha una reconstruccion iterativa)")

print("\n" + "=" * 74 + "\n9. PYPULSEQ\n" + "=" * 74)
ok, err = seq.check_timing()
check("check_timing()", ok, "" if ok else str(err[:3]))
print(seq.test_report())

print("=" * 74)
if fallos:
    print(f"RESULTADO: {len(fallos)} comprobacion(es) FALLAN: " + "; ".join(fallos))
    sys.exit(1)
print("RESULTADO: todas las comprobaciones pasan.")
