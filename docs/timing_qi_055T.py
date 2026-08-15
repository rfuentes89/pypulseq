"""Eleccion del intervalo quiescente (QI) para tsSOS-QISS bSSFP de cuello a 0.55T.

El QI cumple DOS funciones simultaneas y en tension:
  (a) es el TI que nulea el fondo muscular tras la inversion adiabatica;
  (b) es el tiempo durante el cual sangre fresca no saturada entra al slab.
QI corto => mejor supresion de fondo pero menos inflow (critico en vertebrales, lentas).

Valores de relajacion a 0.55T aportados por el usuario:
    sangre  T1/T2 = 1122/263 ms
    grasa   T1/T2 =  187/93  ms
    musculo T1/T2 =  450/55  ms
"""

import numpy as np

# ------------------------------------------------------------------ relajacion a 0.55T
T1_BLOOD, T2_BLOOD = 1.122, 0.263      # s
T1_FAT, T2_FAT = 0.187, 0.093          # s
T1_MUSCLE, T2_MUSCLE = 0.450, 0.055    # s

# ------------------------------------------------------------------ geometria (Koktzoglou MRM 2020)
SLAB_THICKNESS = 19.5e-3               # m
# Velocidades medias promediadas en el ciclo, citadas por Koktzoglou MRM 2020
# a partir de literatura: carotida 20 cm/s, vertebral 10 cm/s.
V_CAROTID, V_VERTEBRAL = 0.20, 0.10    # m/s

# ------------------------------------------------------------------ readout ya verificado
TR_BSSFP = 6.28e-3                     # s, pone la grasa (79.6 Hz) en el stopband
N_KZ_ACQUIRED = 17                     # particiones por angulo radial (18 nominal, 20% OS, PF 6/8)
T_WINDOW = N_KZ_ACQUIRED * TR_BSSFP    # s, ventana de adquisicion por disparo
T_PREP = 25e-3                         # s, inversion adiabatica + spoiler (estimado, se ajusta luego)

print(f"Ventana de adquisicion por disparo: {T_WINDOW*1e3:.1f} ms\n")


def null_ti_steady_state(t1, tr_qiss):
    """TI que nulea un tejido de T1 dado, con recuperacion INCOMPLETA entre inversiones.

    En estado estacionario, la Mz justo antes de cada inversion vale
        x = (1 - E) / (1 + E),   E = exp(-TR/T1)
    y tras invertir (-x) la recuperacion cruza cero en TI = T1 * ln(1 + x).
    Con TR >> T1 esto tiende al clasico T1*ln(2).
    """
    e = np.exp(-tr_qiss / t1)
    x = (1 - e) / (1 + e)
    return t1 * np.log(1 + x)


def mz_after_inversion(t1, ti, tr_qiss):
    """Mz fraccional de un tejido en el instante TI tras la inversion (estado estacionario)."""
    e = np.exp(-tr_qiss / t1)
    x = (1 - e) / (1 + e)
    return 1 - (1 + x) * np.exp(-ti / t1)


print("TR_QISS  TI_null   inflow carotida  inflow vertebral   Mz grasa   Mz sangre")
print("  (ms)     (ms)         (mm)              (mm)            (%)        (%)")
print("-" * 78)

candidates = []
for tr_qiss in np.arange(0.30, 1.55, 0.05):
    ti = null_ti_steady_state(T1_MUSCLE, tr_qiss)

    # Coherencia: el disparo debe caber en el TR.  El centro de k (kz=0, orden ascendente)
    # cae a mitad de ventana, asi que el readout arranca en TI - T_WINDOW/2.
    t_readout_start = ti - T_WINDOW / 2
    if t_readout_start < T_PREP:
        continue                                  # el QI no alcanza ni para el modulo de prep
    if T_PREP + ti + T_WINDOW / 2 > tr_qiss:
        continue                                  # el disparo no cabe en el TR

    inflow_car = V_CAROTID * ti * 1e3             # mm
    inflow_vert = V_VERTEBRAL * ti * 1e3          # mm
    mz_fat = mz_after_inversion(T1_FAT, ti, tr_qiss)
    mz_blood = mz_after_inversion(T1_BLOOD, ti, tr_qiss)

    # Criterio de Koktzoglou MRM 2020: el inflow debe SUPERAR el espesor del slab.
    # Se reporta el margen para que la eleccion sea explicita y no un umbral oculto.
    margen_vert = inflow_vert / (SLAB_THICKNESS * 1e3)
    ok_vert = margen_vert >= 1.0
    mark = f"  {margen_vert:.2f}x" + (" OK" if ok_vert else " <-- insuficiente")
    print(f"{tr_qiss*1e3:6.0f}   {ti*1e3:6.1f}   {inflow_car:11.1f}   {inflow_vert:14.1f}   "
          f"{mz_fat*100:8.1f}   {mz_blood*100:8.1f}{mark}")
    if ok_vert:
        candidates.append((tr_qiss, ti))

print("-" * 78)
print(f"Espesor de slab a refrescar: {SLAB_THICKNESS*1e3:.1f} mm "
      f"(criterio: inflow vertebral >= 1.0x)\n")

if candidates:
    tr_sel, ti_sel = candidates[0]      # el TR mas corto que satisface la vertebral
    print(f"OPERACION ELEGIDA: TR_QISS = {tr_sel*1e3:.0f} ms, QI(TI) = {ti_sel*1e3:.1f} ms")
    print(f"  TI clasico T1*ln2 seria {T1_MUSCLE*np.log(2)*1e3:.1f} ms "
          f"-> la recuperacion incompleta lo acorta en {(T1_MUSCLE*np.log(2)-ti_sel)*1e3:.1f} ms")

    # ------------------------------------------------------------- tiempo de examen
    N_SLABS = 19                        # cobertura de cuello completo (paper 3T: 289 mm)
    n_read = 167                        # muestras por proyeccion a 1.2 mm / FOV 200 mm
    n_views_full = int(np.ceil(np.pi / 2 * n_read))
    print(f"\n  Vistas radiales para Nyquist completo: {n_views_full}")
    print("\n  vistas  submuestreo  t/slab   t total (19 slabs)")
    for n_views in [24, 32, 48, 64, 96]:
        t_slab = n_views * tr_sel
        t_total = t_slab * N_SLABS
        print(f"  {n_views:5d}   {n_views_full/n_views:8.1f}x  {t_slab:6.1f} s   "
              f"{int(t_total//60)}:{int(t_total%60):02d}")
else:
    print("Ningun TR cumple el criterio de inflow vertebral.")
