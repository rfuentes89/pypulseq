"""Estudio de factibilidad del readout bSSFP tsSOS-QISS en MAGNETOM Free.Max 0.55T.

Objetivo: determinar si el TR bSSFP puede alcanzar 6.25 ms (= 1/(2*80 Hz)), valor que
coloca la grasa en la banda de rechazo (stopband) de bSSFP a 0.55T, y con que ancho de
banda de receptor.  NO fija ningun parametro dependiente de T1 (QI/TI quedan fuera).

Hardware confirmado (Varghese et al., Front Cardiovasc Med 2023, doi:10.3389/fcvm.2023.1120982):
max_grad = 26 mT/m, max_slew = 45 T/m/s, bore 80 cm.
"""

import numpy as np
import pypulseq as pp

# ---------------------------------------------------------------- constantes fisicas
B0 = 0.55                      # T, campo del Free.Max
GAMMA_BAR = 42.576e6           # Hz/T, razon giromagnetica del proton
F_LARMOR = GAMMA_BAR * B0      # Hz
PPM_FAT_WATER = 3.4            # ppm, desplazamiento quimico grasa-agua (metileno vs agua)
DF_FAT = PPM_FAT_WATER * 1e-6 * F_LARMOR   # Hz, separacion grasa-agua a 0.55T

# En bSSFP con alternancia de fase RF, los nulos de la banda de paso caen en +-1/(2*TR).
# Poner la grasa exactamente en el nulo exige TR_fat = 1/(2*DF_FAT).
TR_FAT_STOPBAND = 1 / (2 * DF_FAT)

print(f"Larmor a {B0} T          : {F_LARMOR/1e6:.2f} MHz")
print(f"Separacion grasa-agua    : {DF_FAT:.1f} Hz")
print(f"TR que nulea grasa       : {TR_FAT_STOPBAND*1e3:.2f} ms")
print(f"(a 3T seria              : {1/(2*PPM_FAT_WATER*1e-6*GAMMA_BAR*3)*1e3:.2f} ms)\n")

# ---------------------------------------------------------------- limites del sistema
# Derate al 85% del maximo: margen para el redondeo al raster y para PNS.
DERATE = 0.85
system = pp.Opts(
    max_grad=26 * DERATE, grad_unit="mT/m",
    max_slew=45 * DERATE, slew_unit="T/m/s",
    rf_ringdown_time=20e-6,     # confirmado en el sistema real
    rf_dead_time=100e-6,        # confirmado en el sistema real
    adc_dead_time=10e-6,        # confirmado en el sistema real
    grad_raster_time=10e-6,
    rf_raster_time=1e-6,
    B0=B0,
)

# ---------------------------------------------------------------- geometria (Koktzoglou MRM 2020, cuello)
fov_xy = 200e-3        # m, FOV transversal del cuello
slab_thickness = 19.5e-3   # m, espesor de slab del paper de 3T
n_partitions = 18          # particiones nominales por slab (paper: 18, adquiere 14 con PF 6/8)

# Resolucion in-plane: a 0.55T se sacrifica resolucion por SNR respecto al sub-1mm de 3T.
for res_xy in [1.0e-3, 1.2e-3, 1.4e-3]:
    n_read = int(round(fov_xy / res_xy))
    k_max = 1 / (2 * res_xy)          # 1/m
    area_read = 2 * k_max             # 1/m, area total del lobulo de lectura (radial, centro-out completo)

    print(f"--- resolucion in-plane {res_xy*1e3:.1f} mm  ({n_read} muestras/proyeccion) ---")

    # Barrido de anchos de banda de receptor: a bajo campo conviene BW BAJO (mas SNR),
    # lo cual alarga el readout pero encaja bien con el TR largo que queremos igual.
    for bw_px in [300, 400, 500, 700]:
        t_acq = 1 / bw_px                       # s, duracion del readout plano
        t_acq = np.ceil(t_acq / system.grad_raster_time) * system.grad_raster_time
        g_read = area_read / t_acq              # Hz/m requerido en el flat-top
        g_read_mT = g_read / GAMMA_BAR * 1e3    # mT/m

        feasible = g_read <= system.max_grad
        flag = "OK " if feasible else "NO "
        print(f"  {flag} BW={bw_px:4d} Hz/px  T_acq={t_acq*1e3:5.2f} ms  "
              f"G_read={g_read_mT:5.2f} mT/m  (limite {26*DERATE:.1f})")
    print()

# ---------------------------------------------------------------- TR bSSFP realizable
print("=== Construccion de un TR bSSFP completo (res 1.2 mm, BW 400 Hz/px) ===")
res_xy = 1.2e-3
bw_px = 400
n_read = int(round(fov_xy / res_xy))
k_max = 1 / (2 * res_xy)

# --- RF de excitacion, selectiva de slab
rf, gz, gz_rephase = pp.make_sinc_pulse(
    flip_angle=100 * np.pi / 180,   # SAR bajo a 0.55T permite FA alto (Varghese 2023 usa 90-110 en bSSFP)
    duration=600e-6,
    slice_thickness=slab_thickness,
    apodization=0.5,
    time_bw_product=3,
    system=system,
    return_gz=True,
)

# --- lectura radial: prefase (mitad) + lectura completa + rewinder (bSSFP => momento neto cero)
t_acq = np.ceil((1 / bw_px) / system.grad_raster_time) * system.grad_raster_time
g_read = pp.make_trapezoid(channel="x", flat_area=2 * k_max, flat_time=t_acq, system=system)
g_read_pre = pp.make_trapezoid(channel="x", area=-g_read.area / 2, system=system)
adc = pp.make_adc(num_samples=n_read, duration=g_read.flat_time,
                  delay=g_read.rise_time, system=system)

# --- codificacion de particion (kz): peor caso = borde del espacio k
delta_kz = 1 / (n_partitions * (slab_thickness / n_partitions))   # = 1/slab_thickness
area_kz_max = (n_partitions / 2) * (1 / slab_thickness)
g_kz = pp.make_trapezoid(channel="z", area=area_kz_max, system=system)

# --- TR minimo.  Clave: el rephase de slab, la codificacion kz y la prefase de lectura
# viven en ejes distintos (z, z, x/y) y se tocan SIMULTANEAMENTE en un solo bloque.
# El rephase de slab y el kz sí comparten eje z, asi que se suman en un unico trapecio.
t_rf = pp.calc_duration(gz)
g_z_combined = pp.make_trapezoid(
    channel="z", area=gz_rephase.area + g_kz.area, system=system
)
t_prewind = max(pp.calc_duration(g_z_combined), pp.calc_duration(g_read_pre))
t_read = pp.calc_duration(g_read)

tr_min = t_rf + 2 * t_prewind + t_read   # simetrico: bSSFP rebobina todo
te_min = tr_min / 2

print(f"  RF (dur+ramp)         : {t_rf*1e3:.3f} ms")
print(f"  pre/rewind concurrente: {t_prewind*1e3:.3f} ms cada lado")
print(f"     - z (rephase+kz)   : {pp.calc_duration(g_z_combined)*1e3:.3f} ms "
      f"@ {g_z_combined.amplitude/GAMMA_BAR*1e3:.2f} mT/m")
print(f"     - x (prefase lect.): {pp.calc_duration(g_read_pre)*1e3:.3f} ms")
print(f"  lectura               : {t_read*1e3:.3f} ms  (flat {g_read.flat_time*1e3:.2f} ms)")
print(f"  -> TR minimo          : {tr_min*1e3:.3f} ms")
print(f"  -> TE (=TR/2)         : {te_min*1e3:.3f} ms")
print(f"  -> G lectura          : {g_read.amplitude/GAMMA_BAR*1e3:.2f} mT/m")

margen = TR_FAT_STOPBAND - tr_min
print(f"\n  TR objetivo (grasa en stopband): {TR_FAT_STOPBAND*1e3:.2f} ms")
if margen >= 0:
    print(f"  ALCANZABLE: sobran {margen*1e3:.2f} ms -> se rellena con delay o se baja el BW")
    bw_equiv = 1 / (t_acq + margen)
    print(f"  Si el margen se invierte en alargar el readout: BW ~ {bw_equiv:.0f} Hz/px (mas SNR)")
else:
    print(f"  NO alcanzable: faltan {-margen*1e3:.2f} ms")

# ---------------------------------------------------------------- ventana de adquisicion por QISS TR
n_kz_acquired = int(np.ceil(n_partitions * 1.2 * 0.75))   # 20% oversampling, PF 6/8 (paper 3T: 14 de 18)
print(f"\n=== Ventana de adquisicion por QISS TR (inner slice loop) ===")
for tr_bssfp in [tr_min, TR_FAT_STOPBAND]:
    win = n_kz_acquired * tr_bssfp
    print(f"  TR={tr_bssfp*1e3:.2f} ms x {n_kz_acquired} particiones = {win*1e3:.1f} ms por angulo radial")
