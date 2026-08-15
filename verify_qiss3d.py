"""Verificacion de la secuencia 3D tsSOS-QISS bSSFP a 0.55T.

Comprueba lo que check_timing() NO comprueba y que en bSSFP es lo critico:

  1. Balance de momentos de gradiente: en bSSFP la trayectoria de k DEBE volver al
     origen en cada TR.  Un residuo distinto de cero acumula fase entre TRs, rompe
     el estado estacionario y produce bandas/perdida de senal.
  2. Cobertura radial: compara el incremento fijo del paper de 3T, el angulo aureo
     y el reparto uniforme 180/N, porque el del paper esta afinado para SU numero
     de vistas y con las 32 de aqui reparte mal.
  2b. Limites de gradiente por eje frente a magnitud vectorial (relevante para PNS).
  3. Codificacion kz: que los indices de particion recorran el rango esperado con
     partial Fourier 6/8.
  4. test_report() de PyPulseq.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import qiss3d_tssos_bssfp_055T as q

print("Construyendo secuencia de prueba (1 slab, 8 vistas)...")
seq, core = q.build_sequence(n_slabs=1, n_views=8, verbose=False)

# =============================================================================
# 1. Balance de momentos de gradiente en bSSFP
# =============================================================================
print("\n" + "=" * 70)
print("1. BALANCE DE MOMENTOS (condicion bSSFP)")
print("=" * 70)

# Se integra la forma de onda REAL de gradiente (no la trayectoria interpolada) y se
# comprueba que el momento acumulado vuelve a cero al final de cada TR del tren.
# gradient_waveforms() entrega las tres formas de onda muestreadas en el raster de
# gradiente desde t=0, asi que la suma acumulada por el raster ES la trayectoria k.
waveforms, _t_exc, _t_ref, t_adc, _ = seq.waveforms_and_times()
delta_k_xy = 1 / q.FOV_XY

# Cada eje llega como (tiempo, amplitud) en su propia base temporal irregular.
# Se integra por trapecios para obtener k(t) = integral de G dt.
def k_at(axis_wf, t_query):
    t, g = np.asarray(axis_wf)[0], np.asarray(axis_wf)[1]
    k_cum = np.concatenate([[0.0], np.cumsum(np.diff(t) * (g[:-1] + g[1:]) / 2)])
    return np.interp(t_query, t, k_cum)

# Instante en que bSSFP exige k = 0: justo antes de cada RF del tren.  Se localiza
# retrocediendo desde el inicio de cada ADC el tramo RF + prefase.
adc_starts = np.array(sorted(set(np.asarray(t_adc).min(axis=0)
                                 if np.asarray(t_adc).ndim > 1 else np.asarray(t_adc))))
tr_starts = []
t_cursor = 0.0
for i in range(len(seq.block_events)):
    blk = seq.get_block(i + 1)
    if getattr(blk, "adc", None) is not None:
        tr_starts.append(t_cursor - (core.t_rf + core.t_prewind))
    t_cursor += seq.block_durations[i + 1]
tr_starts = np.array(tr_starts)
# Un epsilon antes del RF, para caer dentro del TR anterior ya cerrado.
t_check = np.clip(tr_starts - 1e-7, 0, None)

k_before_rf = np.vstack([k_at(waveforms[a], t_check) for a in range(3)])

# El valor ABSOLUTO de k no es cero: el spoiler que sigue a la inversion adiabatica
# es deliberadamente no balanceado y desplaza el origen para el resto del disparo.
# Lo que bSSFP exige es que el INCREMENTO de k a lo largo de cada TR sea nulo.
# Se compara TR contra TR dentro de cada disparo, excluyendo los saltos entre disparos.
n_tr_shot = q.N_KZ_ACQUIRED
n_shots = len(t_check) // n_tr_shot
dk = []
for s in range(n_shots):
    seg = k_before_rf[:, s * n_tr_shot:(s + 1) * n_tr_shot]
    dk.append(np.diff(seg, axis=1))
dk = np.concatenate(dk, axis=1)

res_xy_k = np.max(np.abs(dk[:2, :]))
res_z_k = np.max(np.abs(dk[2, :]))

print(f"  TRs comprobados                    : {len(t_check)} en {n_shots} disparos")
print(f"  Deriva maxima de |kx,ky| por TR    : {res_xy_k:.3e} 1/m "
      f"({res_xy_k/delta_k_xy*100:.4f} % de delta_k = {delta_k_xy:.2f})")
print(f"  Deriva maxima de |kz|    por TR    : {res_z_k:.3e} 1/m "
      f"({res_z_k/q.DELTA_KZ*100:.4f} % de delta_kz = {q.DELTA_KZ:.2f})")

# Comprobacion analitica independiente: las areas de los eventos deben cancelarse
# exactamente dentro del TR, sin depender del muestreo de la forma de onda.
suma_z = core.gz.area + 2 * core.gz_rephase.area
suma_inplane = core.g_read.area + 2 * core.g_read_pre.area
print(f"  Suma analitica de areas en z       : {suma_z:.3e} 1/m")
print(f"  Suma analitica de areas in-plane   : {suma_inplane:.3e} 1/m")

if abs(suma_z) < 1e-9 and abs(suma_inplane) < 1e-9:
    print("  -> BALANCEADO: los momentos se cancelan exactamente dentro de cada TR")
else:
    print("  -> NO BALANCEADO: hay momento residual, el bSSFP no alcanzara estado estacionario")

# =============================================================================
# 2. Cobertura radial
# =============================================================================
print("\n" + "=" * 70)
print("2. COBERTURA RADIAL")
print("=" * 70)

GOLDEN_ANGLE = 180 * (1 - 1 / ((1 + np.sqrt(5)) / 2))   # 68.7538 deg, periodicidad pi


def gap_stats(increment_deg, n_views):
    angles = np.sort((np.arange(n_views) * increment_deg) % 180.0)
    gaps = np.append(np.diff(angles), 180 - angles[-1] + angles[0])
    ideal = 180 / n_views
    return gaps.min(), gaps.max(), abs(gaps - ideal).max() / ideal * 100


print(f"  paper 3T (fijo) : {q.ANGLE_INCREMENT_PAPER_DEG} deg")
print(f"  angulo aureo    : {GOLDEN_ANGLE:.4f} deg")
print("  uniforme        : 180/N deg (depende de N)\n")
print("          paper 15.8245        angulo aureo         uniforme 180/N")
print("   N     min   max  desv.max    min   max  desv.max    min   max  desv.max")
for n_views in [8, 16, 24, 32, 48, 64]:
    a = gap_stats(q.ANGLE_INCREMENT_PAPER_DEG, n_views)
    b = gap_stats(GOLDEN_ANGLE, n_views)
    c = gap_stats(180 / n_views, n_views)
    print(f"  {n_views:3d}   {a[0]:5.2f} {a[1]:5.2f}  {a[2]:6.1f}%   "
          f"{b[0]:5.2f} {b[1]:5.2f}  {b[2]:6.1f}%   "
          f"{c[0]:5.2f} {c[1]:5.2f}  {c[2]:6.1f}%")
print("\n  desv.max = desviacion maxima del hueco angular respecto del reparto ideal.")
print("  El reparto uniforme es exacto para cualquier N, pero concentra el aliasing")
print("  en streaks coherentes. El aureo reparte peor y a cambio vuelve el aliasing")
print("  incoherente, que es lo que aprovecha una reconstruccion iterativa/CS.")

# =============================================================================
# 2b. Cumplimiento de limites de gradiente: por eje frente a vector
# =============================================================================
print("\n" + "=" * 70)
print("2b. LIMITES DE GRADIENTE: POR EJE vs VECTOR")
print("=" * 70)

g_per_axis = np.array([np.abs(np.asarray(waveforms[a])[1]).max() for a in range(3)])
g_per_axis_mt = g_per_axis / q.GAMMA_BAR * 1e3

# Magnitud vectorial: hay que remuestrear los tres ejes a una base comun.
t_common = np.linspace(0, min(np.asarray(waveforms[a])[0].max() for a in range(3)), 200000)
g_interp = np.vstack([np.interp(t_common, np.asarray(waveforms[a])[0],
                                np.asarray(waveforms[a])[1]) for a in range(3)])
g_vec_mt = np.linalg.norm(g_interp, axis=0).max() / q.GAMMA_BAR * 1e3
slew_vec = np.linalg.norm(np.gradient(g_interp, t_common, axis=1), axis=0).max()
slew_vec_t = slew_vec / q.GAMMA_BAR

print(f"  Nameplate Free.Max            : {q.MAX_GRAD_MT_M:.0f} mT/m, "
      f"{q.MAX_SLEW_T_M_S:.0f} T/m/s (por eje)")
print(f"  Limite derateado ({q.GRAD_DERATE:.0%})        : "
      f"{q.MAX_GRAD_MT_M*q.GRAD_DERATE:.1f} mT/m, {q.MAX_SLEW_T_M_S*q.GRAD_DERATE:.1f} T/m/s")
print(f"  Maximo por eje  (x, y, z)     : {g_per_axis_mt[0]:.2f}, "
      f"{g_per_axis_mt[1]:.2f}, {g_per_axis_mt[2]:.2f} mT/m")
if (g_per_axis_mt <= q.MAX_GRAD_MT_M * q.GRAD_DERATE + 1e-6).all():
    print("  -> CUMPLE por eje, que es como especifica el fabricante el gradiente")
else:
    print("  -> EXCEDE el limite por eje")
print(f"  Magnitud vectorial maxima     : {g_vec_mt:.2f} mT/m")
print(f"  Slew vectorial maximo         : {slew_vec_t:.1f} T/m/s")
print("\n  El vector supera el nameplate por eje al sumar ejes en diagonal, cosa normal")
print("  en 3D: los amplificadores son independientes por eje. Lo que SI depende del")
print("  vector es la estimulacion nerviosa periferica (PNS), que PyPulseq no modela.")
print("  ANTES DE IR AL SCANNER: pasar el .seq por el chequeo de PNS del fabricante.")
print("  El bore de 80 cm del Free.Max implica una bobina de gradiente grande, lo que")
print("  tiende a EMPEORAR el PNS por unidad de slew: no darlo por descontado.")

# =============================================================================
# 3. Codificacion de particion
# =============================================================================
print("\n" + "=" * 70)
print("3. CODIFICACION kz (partial Fourier)")
print("=" * 70)
print(f"  Grilla nominal        : {q.N_KZ_GRID} pasos")
print(f"  Adquiridos            : {q.N_KZ_ACQUIRED} (PF {q.PARTIAL_FOURIER_KZ:.3f})")
print(f"  Indices               : {q.KZ_INDICES[0]} .. {q.KZ_INDICES[-1]}")
print(f"  FOV_z (con OS 20%)    : {q.FOV_Z*1e3:.2f} mm")
print(f"  Espesor de particion  : {q.PARTITION_THICKNESS*1e3:.3f} mm")
print(f"  delta_kz              : {q.DELTA_KZ:.2f} 1/m")
kz_max_acquired = np.abs(q.KZ_INDICES).max() * q.DELTA_KZ
print(f"  |kz| maximo           : {kz_max_acquired:.1f} 1/m "
      f"-> resolucion z = {1/(2*kz_max_acquired)*1e3:.3f} mm")

# Comprobacion de coherencia: el lado negativo debe estar completo (es el que
# sostiene la reconstruccion hermitica del lado recortado).
if q.KZ_INDICES[0] == -q.N_KZ_GRID // 2:
    print("  -> lado negativo completo, recorte en el positivo: correcto para PF 6/8")
else:
    print("  -> ATENCION: el recorte de partial Fourier no esta donde se espera")

# =============================================================================
# 4. Graficos
# =============================================================================
k_traj_adc = seq.calculate_kspace()[0]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Trayectoria in-plane de las 8 vistas
axes[0].plot(k_traj_adc[0, :], k_traj_adc[1, :], ".", ms=0.4, alpha=0.5)
axes[0].set_xlabel("$k_x$ [1/m]"); axes[0].set_ylabel("$k_y$ [1/m]")
axes[0].set_title("Stack-of-stars: proyeccion in-plane")
axes[0].set_aspect("equal"); axes[0].grid(alpha=0.3)

# kz frente a numero de muestra
axes[1].plot(k_traj_adc[2, :], lw=0.5)
axes[1].set_xlabel("muestra ADC"); axes[1].set_ylabel("$k_z$ [1/m]")
axes[1].set_title("Codificacion de particion (inner slice loop)")
axes[1].grid(alpha=0.3)

# Distribucion angular
angles_full = np.degrees([q.view_angle(i, q.N_VIEWS) for i in range(q.N_VIEWS)])
axes[2].vlines(np.sort(angles_full), 0, 1, lw=1)
axes[2].set_xlabel("angulo de vista [deg]"); axes[2].set_xlim(0, 180)
axes[2].set_yticks([])
axes[2].set_title(f"Reparto de {q.N_VIEWS} vistas ({q.VIEW_ORDER})")

plt.tight_layout()
plt.savefig("qiss3d_verificacion.png", dpi=130)
print("\nGrafico guardado: qiss3d_verificacion.png")

# =============================================================================
# 5. test_report
# =============================================================================
print("\n" + "=" * 70)
print("4. test_report() DE PYPULSEQ")
print("=" * 70)
print(seq.test_report())
