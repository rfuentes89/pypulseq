"""3D tsSOS-QISS MRA de cuello para MAGNETOM Free.Max 0.55T (PyPulseq 1.4.2).

Angiografia por resonancia sin contraste, thin-slab stack-of-stars QISS, readout bSSFP,
sin gating cardiaco.  Adaptacion a bajo campo de:

  * Koktzoglou et al., Magn Reson Med 2020;84:3316-3324, doi:10.1002/mrm.28339
    -> geometria de cuello (19 slabs solapados), esquema ungated, doble inversion
       (fondo sobre el slab + "tracking" venosa por encima), orden inner-slice-loop.
  * Edelman et al., Magn Reson Med 2020;83:1711-1720, doi:10.1002/mrm.28032
    -> readout bSSFP y estrategia de solape de slabs contra el venetian blind artifact.

Hardware (Varghese et al., Front Cardiovasc Med 2023, doi:10.3389/fcvm.2023.1120982):
    MAGNETOM Free.Max, 0.55T, bore 80 cm, gradientes 26 mT/m y 45 T/m/s.

TRES DESVIACIONES DELIBERADAS respecto de los papers, todas forzadas por el bajo campo:

  1. Supresion grasa por stopband de bSSFP, no por TEs en oposicion de fase.
     A 0.55T la separacion grasa-agua es 79.6 Hz, asi que TR = 1/(2*79.6) = 6.28 ms
     coloca la grasa exactamente en el nulo de la banda de paso.  A 3T ese truco
     exigiria TR = 1.15 ms (irrealizable), y por eso Koktzoglou usa TEs multiples.
     Feliz coincidencia: 6.28 ms es tambien ~ el TR minimo que permiten estos gradientes.

  2. QISS TR de 600 ms en vez de los 1500 ms del paper de 3T.
     Con readout bSSFP de flip angle alto no se puede sostener una ventana de
     adquisicion de 891 ms sin saturar la sangre entrante y destruir el contraste
     de inflow.  Se usa una ventana corta (un angulo radial por disparo) y se
     recupera tiempo de examen acortando el TR.

  3. El QI se calcula por nulling de estado ESTACIONARIO, no por T1*ln2.
     Con TR 600 ms y T1 muscular 450 ms el fondo no recupera del todo entre
     inversiones: el TI de nulo baja de 312 ms a 207 ms.  Ver docs/qiss_timing.md.

Unidades SI en todo el script: Hz, s, T/m, m.  PyPulseq trabaja en Hz/m para
gradientes; toda conversion desde mT/m es explicita.
"""

from types import SimpleNamespace

import numpy as np
import pypulseq as pp

# =============================================================================
# 1. CONSTANTES FISICAS Y DEL SISTEMA
# =============================================================================

B0 = 0.55                       # T
GAMMA_BAR = 42.576e6            # Hz/T, razon giromagnetica del proton
F_LARMOR = GAMMA_BAR * B0       # Hz -> 23.42 MHz

# Desplazamiento quimico grasa-agua (metileno). 3.4 ppm es el valor convencional
# para el pico principal de grasa; a 0.55T da 79.6 Hz.
PPM_FAT_WATER = 3.4
DF_FAT = PPM_FAT_WATER * 1e-6 * F_LARMOR

# --- Tiempos de relajacion a 0.55T (aportados por el usuario, no extrapolados de 1.5T/3T)
T1_BLOOD, T2_BLOOD = 1.122, 0.263      # s
T1_FAT, T2_FAT = 0.187, 0.093          # s
T1_MUSCLE, T2_MUSCLE = 0.450, 0.055    # s

# --- PENDIENTE DE CONFIRMAR CON EL SISTEMA REAL -----------------------------
# Estos tres tiempos muertos no estan publicados para el Free.Max y NO deben
# darse por buenos sin leerlos del scanner.  Los valores de abajo son
# conservadores; si los reales son menores, el TR minimo baja y sobra margen.
RF_DEAD_TIME = 100e-6           # s
RF_RINGDOWN_TIME = 30e-6        # s
ADC_DEAD_TIME = 10e-6           # s
# ----------------------------------------------------------------------------

# Derate al 85% del maximo de placa: cubre el redondeo al raster y deja margen
# de PNS.  Con 45 T/m/s el Free.Max dificilmente alcanza el umbral de PNS, pero
# el derate evita que check_timing falle por un trapecio marginalmente ilegal.
GRAD_DERATE = 0.85
MAX_GRAD_MT_M = 26.0
MAX_SLEW_T_M_S = 45.0

system = pp.Opts(
    max_grad=MAX_GRAD_MT_M * GRAD_DERATE, grad_unit="mT/m",
    max_slew=MAX_SLEW_T_M_S * GRAD_DERATE, slew_unit="T/m/s",
    rf_dead_time=RF_DEAD_TIME,
    rf_ringdown_time=RF_RINGDOWN_TIME,
    adc_dead_time=ADC_DEAD_TIME,
    grad_raster_time=10e-6,
    rf_raster_time=1e-6,
    B0=B0,
)

# =============================================================================
# 2. PARAMETROS DE LA SECUENCIA
# =============================================================================

# --- Geometria (Koktzoglou MRM 2020, cuello completo)
FOV_XY = 200e-3                 # m, FOV transversal
RES_XY = 1.2e-3                 # m, resolucion in-plane. A 3T el paper baja de 1 mm;
                                # a 0.55T se cede resolucion por SNR de forma deliberada.
SLAB_THICKNESS = 19.5e-3        # m, espesor de slab excitado
SLAB_OVERLAP = 4.5e-3           # m, solape entre slabs contiguos
N_SLABS = 19                    # -> cobertura axial de 289 mm, igual que el paper
SLICE_OVERSAMPLING = 1.20       # 20% de oversampling en direccion de particion
N_KZ_GRID = 18                  # pasos de codificacion en la grilla ya sobremuestreada
PARTIAL_FOURIER_KZ = 6 / 8      # -> se adquieren 14 de 18 pasos

# --- Timing (derivado en docs/qiss_timing.md, no numeros magicos)
TR_BSSFP = 6.28e-3              # s, = 1/(2*DF_FAT): grasa en el stopband de bSSFP
TR_QISS = 600e-3                # s, periodo de repeticion del modulo QISS completo
FLIP_ANGLE_DEG = 100            # deg. El SAR cae con B0^2, asi que a 0.55T se puede
                                # usar 90-110 deg en bSSFP (Varghese 2023 lo hace).

# --- Muestreo radial
N_VIEWS = 32                    # angulos radiales por slab (submuestreo ~8x)

# Orden de vistas radiales.  El incremento de 15.8245 deg del paper de 3T se eligio
# para SU numero de vistas; con las 32 de aqui reparte mal (el hueco angular maximo
# llega a 1.8x el ideal).  Ver verify_qiss3d.py para la comparacion numerica.
#   "uniform" -> 180/N. Reparto exacto para cualquier N, mejor PSF con gridding/NUFFT
#                simple, pero el aliasing sale como streaks coherentes.
#   "golden"  -> angulo aureo. Reparte algo peor y a cambio vuelve el aliasing
#                incoherente, que es lo que explota una reconstruccion iterativa/CS.
#   "paper"   -> reproduce literalmente el incremento de Koktzoglou MRM 2020.
# Por defecto "golden": a 0.55T el submuestreo de 8x casi con seguridad exige recon
# iterativa (Varghese 2023 necesito CS para todo en este scanner), y ahi el aliasing
# incoherente rinde mejor que el reparto perfecto.
VIEW_ORDER = "golden"
ANGLE_INCREMENT_PAPER_DEG = 15.8245
GOLDEN_ANGLE_DEG = 180 * (1 - 1 / ((1 + np.sqrt(5)) / 2))   # 68.7539 deg


def view_angle(i_view, n_views, order=VIEW_ORDER):
    """Angulo de la vista radial i, en radianes."""
    if order == "uniform":
        inc = 180.0 / n_views
    elif order == "golden":
        inc = GOLDEN_ANGLE_DEG
    elif order == "paper":
        inc = ANGLE_INCREMENT_PAPER_DEG
    else:
        raise ValueError(f"VIEW_ORDER desconocido: {order}")
    return np.deg2rad((i_view * inc) % 180.0)

# --- Modulos de supresion
INV_SLAB_THICKNESS = SLAB_THICKNESS   # inversion de fondo, coextensiva con el slab
INV_VENOUS_THICKNESS = 100e-3         # m, inversion "tracking" venosa (paper: 10 cm)
INV_VENOUS_GAP = 5e-3                 # m, hueco entre borde superior del slab y la banda
INV_DURATION = 8e-3                   # s, duracion del pulso adiabatico
SPOILER_CYCLES_PER_VOXEL = 4          # ciclos de fase por voxel tras la inversion

# Direccion de flujo: en el cuello la sangre arterial sube (caudal -> craneal, +z),
# la venosa baja.  Por eso la banda venosa va POR ENCIMA del slab.
Z_SUPERIOR = +1


# =============================================================================
# 3. BLOQUES REUTILIZABLES
# =============================================================================

def make_adiabatic_inversion(sys, thickness, z_center, duration=INV_DURATION):
    """Inversion adiabatica selectiva de slab, insensible a inhomogeneidad de B1.

    Sustituye al pulso FOCI de los papers por un hyperbolic secant, que es lo que
    ofrece PyPulseq 1.4.2 de fabrica.  DIFERENCIA A TENER PRESENTE: el FOCI da un
    perfil de slab de bordes mas abruptos que el HS, y Edelman (MRM 2020) apoya su
    version mejorada justamente en esa nitidez para no saturar la sangre entrante.
    Con HS conviene verificar por simulacion el perfil antes de ir al scanner.

    El SAR de un adiabatico es alto en campo alto, pero a 0.55T escala con B0^2 y
    deja de ser el factor limitante.
    """
    rf, gz, _ = pp.make_adiabatic_pulse(
        pulse_type="hypsec",
        duration=duration,
        slice_thickness=thickness,
        system=sys,
        return_gz=True,
        use="inversion",
    )
    # Desplazar la banda a su posicion axial mediante offset de frecuencia:
    # df = gamma * G * z, con G la amplitud del gradiente selector en Hz/m.
    rf.freq_offset = gz.amplitude * z_center
    return rf, gz


def make_spoiler(sys, voxel_size_m, n_cycles=SPOILER_CYCLES_PER_VOXEL, axis="z"):
    """Spoiler que dispersa n_cycles ciclos completos de fase a lo largo de un voxel.

    n_cycles = 4 es el punto de partida habitual: suficiente para deshacer la
    coherencia transversal residual tras la inversion sin gastar tiempo de mas.
    El area va en 1/m, que es la convencion de PyPulseq para momentos de gradiente.
    """
    area = n_cycles / voxel_size_m
    return pp.make_trapezoid(channel=axis, area=area, system=sys)


def make_bssfp_core(sys, fov_xy, res_xy, slab_thickness, flip_angle_deg, tr_target):
    """Construye los eventos invariantes del TR bSSFP (todo menos el giro radial y kz).

    Devuelve un namespace con el RF, el gradiente selector, el trapecio de lectura,
    su prefase, el ADC y los tiempos de relleno necesarios para clavar tr_target.

    bSSFP exige momento de gradiente NETO CERO en los tres ejes dentro del TR.  Se
    consigue con: gz + 2*gz_rephase = 0 en z, y prefase = -mitad del area de lectura
    a cada lado en el plano.
    """
    # Multiplo de 4: los interpretes de Siemens exigen que el numero de muestras
    # del ADC sea divisible por 4.
    n_read = int(round(fov_xy / res_xy / 4)) * 4
    k_max = 1 / (2 * res_xy)                       # 1/m

    rf, gz, gz_rephase = pp.make_sinc_pulse(
        flip_angle=flip_angle_deg * np.pi / 180,
        duration=600e-6,          # RF corto: en bSSFP cada us de RF sale del presupuesto de TR
        slice_thickness=slab_thickness,
        apodization=0.5,
        time_bw_product=3,        # compromiso perfil/duracion; TBW mayor afila el slab
                                  # pero alarga el RF y sube el TR minimo
        system=sys,
        return_gz=True,
    )

    # El readout ocupa el tiempo que sobra tras RF y rebobinados.  A 0.55T interesa
    # que sea LARGO (ancho de banda bajo = mas SNR), asi que se calcula por diferencia
    # en vez de fijar un BW arbitrario.
    g_read_tmp = pp.make_trapezoid(channel="x", flat_area=2 * k_max,
                                   flat_time=2.5e-3, system=sys)
    g_read_pre_tmp = pp.make_trapezoid(channel="x", area=-g_read_tmp.area / 2, system=sys)

    t_rf = pp.calc_duration(gz)
    t_prewind = pp.calc_duration(g_read_pre_tmp)
    t_read_available = tr_target - t_rf - 2 * t_prewind

    flat_time = t_read_available - 2 * g_read_tmp.rise_time
    flat_time = np.floor(flat_time / sys.grad_raster_time) * sys.grad_raster_time
    if flat_time <= 0:
        raise ValueError(
            f"TR objetivo {tr_target*1e3:.2f} ms no alcanza: RF+rebobinados ya consumen "
            f"{(t_rf + 2*t_prewind)*1e3:.2f} ms"
        )

    g_read = pp.make_trapezoid(channel="x", flat_area=2 * k_max,
                               flat_time=flat_time, system=sys)
    g_read_pre = pp.make_trapezoid(channel="x", area=-g_read.area / 2, system=sys)

    # El dwell debe caer en el raster del ADC (100 ns).  Se redondea hacia abajo y el
    # ADC resultante, algo mas corto que el flat-top, se CENTRA en el: muestrear solo
    # dentro de la meseta evita tomar puntos durante las rampas, donde el gradiente
    # todavia no es constante y la trayectoria radial se deformaria.
    adc_raster = sys.adc_raster_time
    dwell = np.floor(g_read.flat_time / n_read / adc_raster) * adc_raster
    adc_duration = n_read * dwell
    adc_offset = np.round((g_read.flat_time - adc_duration) / 2 / adc_raster) * adc_raster
    adc = pp.make_adc(num_samples=n_read, dwell=dwell,
                      delay=g_read.rise_time + adc_offset, system=sys)

    bw_per_pixel = 1 / adc_duration

    return SimpleNamespace(
        rf=rf, gz=gz, gz_rephase=gz_rephase,
        g_read=g_read, g_read_pre=g_read_pre, adc=adc,
        n_read=n_read, k_max=k_max, bw_per_pixel=bw_per_pixel,
        t_rf=t_rf, t_prewind=pp.calc_duration(g_read_pre),
    )


def rotate_in_plane(grad_proto, theta, sys):
    """Proyecta un gradiente prototipo del eje x sobre (x, y) para el angulo radial theta."""
    gx = pp.scale_grad(grad_proto, np.cos(theta))
    gy = pp.scale_grad(pp.make_trapezoid(
        channel="y",
        amplitude=grad_proto.amplitude,
        flat_time=getattr(grad_proto, "flat_time", 0),
        rise_time=grad_proto.rise_time,
        fall_time=grad_proto.fall_time,
        system=sys,
    ), np.sin(theta))
    return gx, gy


# =============================================================================
# 4. TIMING DERIVADO
# =============================================================================

def null_ti_steady_state(t1, tr_qiss):
    """TI de nulo con recuperacion incompleta entre inversiones sucesivas.

    En estado estacionario la Mz justo antes de cada inversion vale
        x = (1 - E)/(1 + E),  con E = exp(-TR/T1),
    y tras invertir cruza cero en TI = T1*ln(1 + x).  Para TR >> T1 tiende a T1*ln2.
    Ignorar esto y usar T1*ln2 dejaria el fondo sin nulear por ~105 ms de exceso.
    """
    e = np.exp(-tr_qiss / t1)
    return t1 * np.log(1 + (1 - e) / (1 + e))


N_KZ_ACQUIRED = int(round(N_KZ_GRID * PARTIAL_FOURIER_KZ))     # 14 de 18
FOV_Z = SLAB_THICKNESS * SLICE_OVERSAMPLING                    # 23.4 mm
PARTITION_THICKNESS = FOV_Z / N_KZ_GRID                        # 1.3 mm
DELTA_KZ = 1 / FOV_Z                                           # 1/m

QI = null_ti_steady_state(T1_MUSCLE, TR_QISS)                  # s, ~207 ms
T_ACQ_WINDOW = N_KZ_ACQUIRED * TR_BSSFP                        # s, ~88 ms

# Con orden de particion ascendente, kz = 0 cae a mitad de la ventana.  El QI se
# mide desde la inversion hasta el centro del espacio k, asi que el readout arranca
# QI - ventana/2 despues de la inversion.
T_READOUT_START = QI - T_ACQ_WINDOW / 2

# Indices de particion con partial Fourier 6/8: se conserva el lado negativo completo
# y se recorta el positivo, que es donde el conjugado hermitico permite reconstruir.
KZ_INDICES = np.arange(-N_KZ_GRID // 2, -N_KZ_GRID // 2 + N_KZ_ACQUIRED)


# =============================================================================
# 5. ENSAMBLADO
# =============================================================================

def build_sequence(n_slabs=N_SLABS, n_views=N_VIEWS, verbose=True):
    seq = pp.Sequence(system=system)
    core = make_bssfp_core(system, FOV_XY, RES_XY, SLAB_THICKNESS,
                           FLIP_ANGLE_DEG, TR_BSSFP)

    spoiler = make_spoiler(system, PARTITION_THICKNESS)

    # Pulso alpha/2 para entrar al estado estacionario de bSSFP sin oscilacion
    # transitoria (Le Roux).  Sin el, las primeras vistas radiales traen artefacto.
    rf_half, gz_half, gz_half_rephase = pp.make_sinc_pulse(
        flip_angle=FLIP_ANGLE_DEG / 2 * np.pi / 180,
        duration=600e-6, slice_thickness=SLAB_THICKNESS,
        apodization=0.5, time_bw_product=3, system=system, return_gz=True,
    )

    # Posicion axial del centro de cada slab
    slab_step = SLAB_THICKNESS - SLAB_OVERLAP
    z_centers = (np.arange(n_slabs) - (n_slabs - 1) / 2) * slab_step

    for i_slab, z_c in enumerate(z_centers):
        # Offset de frecuencia que traslada RF y alpha/2 al slab i
        freq_slab = core.gz.amplitude * z_c
        core.rf.freq_offset = freq_slab
        rf_half.freq_offset = freq_slab

        # Banda venosa: por encima del slab, con hueco de seguridad
        z_venous = z_c + Z_SUPERIOR * (SLAB_THICKNESS / 2 + INV_VENOUS_GAP
                                       + INV_VENOUS_THICKNESS / 2)

        for i_view in range(n_views):
            theta = view_angle(i_view, n_views)

            # ---------------- modulo de preparacion ----------------
            rf_inv_bg, gz_inv_bg = make_adiabatic_inversion(
                system, INV_SLAB_THICKNESS, z_c)
            seq.add_block(rf_inv_bg, gz_inv_bg)

            rf_inv_ven, gz_inv_ven = make_adiabatic_inversion(
                system, INV_VENOUS_THICKNESS, z_venous)
            seq.add_block(rf_inv_ven, gz_inv_ven)

            seq.add_block(spoiler)

            # ---------------- intervalo quiescente ----------------
            t_prep = (pp.calc_duration(gz_inv_bg) + pp.calc_duration(gz_inv_ven)
                      + pp.calc_duration(spoiler))
            t_qi_delay = T_READOUT_START - t_prep - TR_BSSFP / 2   # alpha/2 ocupa TR/2
            if t_qi_delay < 0:
                raise ValueError(
                    f"QI insuficiente: la preparacion ({t_prep*1e3:.1f} ms) excede "
                    f"el inicio de readout ({T_READOUT_START*1e3:.1f} ms)"
                )
            seq.add_block(pp.make_delay(
                np.round(t_qi_delay / system.grad_raster_time) * system.grad_raster_time))

            # ---------------- preparacion alpha/2 ----------------
            rf_half.phase_offset = np.pi          # opuesta al primer alpha del tren
            seq.add_block(rf_half, gz_half)
            seq.add_block(gz_half_rephase, pp.make_delay(TR_BSSFP / 2
                                                          - pp.calc_duration(gz_half)))

            # ---------------- tren bSSFP: inner slice loop ----------------
            gx_read, gy_read = rotate_in_plane(core.g_read, theta, system)
            gx_pre, gy_pre = rotate_in_plane(core.g_read_pre, theta, system)

            for i_tr, kz_idx in enumerate(KZ_INDICES):
                phase = np.pi * (i_tr % 2)        # alternancia 0 / 180 de bSSFP
                core.rf.phase_offset = phase
                core.adc.phase_offset = phase

                seq.add_block(core.rf, core.gz)

                # kz y rephase de slab comparten eje z: se suman en un trapecio.
                # Van concurrentes con la prefase de lectura (x, y) -> ahorra ~1 ms
                # de TR, que es lo que hace alcanzable el TR de 6.28 ms.
                area_z = core.gz_rephase.area + kz_idx * DELTA_KZ
                gz_pre = pp.make_trapezoid(channel="z", area=area_z, system=system)
                seq.add_block(gz_pre, gx_pre, gy_pre)

                seq.add_block(gx_read, gy_read, core.adc)

                area_z_post = core.gz_rephase.area - kz_idx * DELTA_KZ
                gz_post = pp.make_trapezoid(channel="z", area=area_z_post, system=system)
                seq.add_block(gz_post, gx_pre, gy_pre)

            # ---------------- relleno hasta completar el QISS TR ----------------
            t_shot = (t_prep + t_qi_delay + TR_BSSFP / 2
                      + N_KZ_ACQUIRED * TR_BSSFP)
            t_fill = TR_QISS - t_shot
            if t_fill < 0:
                raise ValueError(
                    f"El disparo ({t_shot*1e3:.1f} ms) no cabe en TR_QISS "
                    f"({TR_QISS*1e3:.1f} ms)"
                )
            seq.add_block(pp.make_delay(
                np.round(t_fill / system.grad_raster_time) * system.grad_raster_time))

        if verbose:
            print(f"  slab {i_slab+1}/{n_slabs} en z = {z_c*1e3:+7.1f} mm  listo")

    return seq, core


# =============================================================================
# 6. INFORME Y VERIFICACION
# =============================================================================

def report(core, n_slabs, n_views):
    print("\n" + "=" * 72)
    print("  3D tsSOS-QISS bSSFP  |  MAGNETOM Free.Max 0.55T")
    print("=" * 72)
    print(f"  Larmor                      : {F_LARMOR/1e6:.2f} MHz")
    print(f"  Separacion grasa-agua       : {DF_FAT:.1f} Hz")
    print(f"  TR bSSFP (grasa en stopband): {TR_BSSFP*1e3:.2f} ms  "
          f"[objetivo 1/(2*df) = {1/(2*DF_FAT)*1e3:.2f} ms]")
    print(f"  TE (= TR/2)                 : {TR_BSSFP/2*1e3:.2f} ms")
    print(f"  Flip angle                  : {FLIP_ANGLE_DEG} deg")
    print(f"  Ancho de banda receptor     : {core.bw_per_pixel:.0f} Hz/px")
    print(f"  Muestras por proyeccion     : {core.n_read}")
    print("-" * 72)
    print(f"  QISS TR                     : {TR_QISS*1e3:.0f} ms  (ungated)")
    print(f"  QI / TI de nulo             : {QI*1e3:.1f} ms")
    print(f"     T1*ln2 ingenuo seria     : {T1_MUSCLE*np.log(2)*1e3:.1f} ms "
          f"(exceso de {(T1_MUSCLE*np.log(2)-QI)*1e3:.1f} ms si se ignora el estado estacionario)")
    print(f"  Ventana de adquisicion      : {T_ACQ_WINDOW*1e3:.1f} ms")
    print(f"  Inicio de readout post-inv. : {T_READOUT_START*1e3:.1f} ms")
    print("-" * 72)
    print(f"  Slabs                       : {n_slabs} x {SLAB_THICKNESS*1e3:.1f} mm, "
          f"solape {SLAB_OVERLAP*1e3:.1f} mm")
    print(f"  Cobertura axial             : "
          f"{((n_slabs-1)*(SLAB_THICKNESS-SLAB_OVERLAP) + SLAB_THICKNESS)*1e3:.0f} mm")
    print(f"  Particiones                 : {N_KZ_ACQUIRED} de {N_KZ_GRID} "
          f"(PF {PARTIAL_FOURIER_KZ:.3f}, OS {(SLICE_OVERSAMPLING-1)*100:.0f}%)")
    print(f"  Resolucion adquirida        : {RES_XY*1e3:.2f} x {RES_XY*1e3:.2f} x "
          f"{PARTITION_THICKNESS*1e3:.2f} mm")
    print(f"  Vistas radiales             : {n_views} "
          f"(Nyquist pleno = {int(np.ceil(np.pi/2*core.n_read))}, "
          f"submuestreo {np.ceil(np.pi/2*core.n_read)/n_views:.1f}x)")
    print(f"  Orden de vistas             : {VIEW_ORDER}")
    t_total = n_slabs * n_views * TR_QISS
    print(f"  Tiempo de examen            : {int(t_total//60)}:{int(t_total%60):02d} "
          f"({t_total:.0f} s)")
    print("-" * 72)
    print("  Comportamiento esperado de la magnetizacion en el QI:")
    for name, t1 in [("musculo (fondo)", T1_MUSCLE), ("grasa", T1_FAT),
                     ("sangre estatica/venosa", T1_BLOOD)]:
        e = np.exp(-TR_QISS / t1)
        x = (1 - e) / (1 + e)
        mz = 1 - (1 + x) * np.exp(-QI / t1)
        print(f"    {name:24s}: Mz = {mz*100:+6.1f} %")
    print("    sangre arterial entrante: Mz = +100.0 %  (nunca vio la inversion)")
    print("=" * 72 + "\n")


if __name__ == "__main__":
    import sys as _sys

    quick = "--quick" in _sys.argv
    n_slabs = 2 if quick else N_SLABS
    n_views = 4 if quick else N_VIEWS
    if quick:
        print(">> modo rapido: 2 slabs x 4 vistas (solo para verificacion de timing)\n")

    print("Construyendo secuencia...")
    seq, core = build_sequence(n_slabs=n_slabs, n_views=n_views, verbose=quick)
    report(core, n_slabs, n_views)

    print("Verificando timing...")
    ok, error_report = seq.check_timing()
    if ok:
        print("  check_timing(): OK")
    else:
        print("  check_timing(): PROBLEMAS")
        for e in error_report[:20]:
            print("   ", e)

    seq.set_definition("FOV", [FOV_XY, FOV_XY,
                               (n_slabs - 1) * (SLAB_THICKNESS - SLAB_OVERLAP)
                               + SLAB_THICKNESS])
    seq.set_definition("Name", "qiss3d_tssos_bssfp_055T")
    seq.set_definition("TR_bSSFP_ms", TR_BSSFP * 1e3)
    seq.set_definition("TR_QISS_ms", TR_QISS * 1e3)
    seq.set_definition("QI_ms", QI * 1e3)
    seq.set_definition("FlipAngle_deg", FLIP_ANGLE_DEG)
    seq.set_definition("RadialViews", n_views)
    seq.set_definition("ViewOrder", VIEW_ORDER)
    seq.set_definition("Slabs", n_slabs)
    seq.set_definition("B0_T", B0)

    out = "qiss3d_tssos_bssfp_055T_quick.seq" if quick else "qiss3d_tssos_bssfp_055T.seq"
    seq.write(out)
    print(f"\nExportado: {out}")
