"""3D tsSOS-QISS MRA de cuello para MAGNETOM Free.Max 0.55T (PyPulseq 1.4.2).

Angiografia sin contraste, thin-slab stack-of-stars QISS, readout bSSFP, sin gating.
Adaptacion a bajo campo de:

  * Koktzoglou et al., Magn Reson Med 2020;84:3316-3324, doi:10.1002/mrm.28339
    -> geometria de cuello (19 slabs solapados), esquema ungated, inner-slice-loop.
  * Edelman et al., Magn Reson Med 2020;83:1711-1720, doi:10.1002/mrm.28032
    -> readout bSSFP, slab excitado mas grueso que el conservado.

Hardware (Varghese et al., Front Cardiovasc Med 2023, doi:10.3389/fcvm.2023.1120982):
    MAGNETOM Free.Max, 0.55T, bore 80 cm, gradientes 26 mT/m y 45 T/m/s.

Disparo:  inversion de fondo (WURST) -> saturacion venosa -> QI -> fat-sat
          -> alpha/2 -> 15 TR bSSFP (un angulo radial, todas las particiones).

DESVIACIONES DELIBERADAS respecto de los papers (detalle en docs/qiss_timing.md):

  1. QISS TR de 700 ms, no 1500.  Con bSSFP de flip alto una ventana de 891 ms
     saturaria la sangre entrante; se usa una ventana corta y un TR corto.
  2. El QI sale del nulo en estado ESTACIONARIO con inversion imperfecta, no de
     T1*ln2: 207 ms en vez de 312.
  3. Venas por SATURACION, no por inversion: con un QI de 200 ms una inversion
     unica deja la sangre que fluye en Mz = -0.66, brillante en magnitud.
  4. Slab excitado 1.2x mas grueso y 33 % de oversampling en kz: el RF corto que
     exige bSSFP tiene un perfil blando.
  5. Grasa: fat-sat espectral + TR que la acerca al nulo del stopband.  Los TE en
     oposicion de fase del paper de 3T aqui quedarian casi en fase (79.6 Hz).

Esta version corrige los errores de una auditoria independiente; las cifras de
contraste salen de sim_bssfp_contrast.py y los perfiles de sim_inversion_profile.py.
Unidades SI: Hz, s, T/m, m (los gradientes de PyPulseq van en Hz/m).
"""

import os
import re
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

# --- Tiempos muertos CONFIRMADOS en el sistema real -------------------------
# Leidos del Free.Max del usuario (no son valores de catalogo ni estimaciones).
RF_DEAD_TIME = 100e-6           # s
RF_RINGDOWN_TIME = 20e-6        # s
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

# Version del interprete Pulseq instalado en el scanner del usuario.  PyPulseq 1.4.2
# escribe formato 1.4.2, asi que hay coincidencia exacta y no hace falta degradar
# ninguna feature.  Se comprueba en tiempo de ejecucion para que un futuro upgrade
# de pypulseq no rompa la compatibilidad en silencio.
SCANNER_INTERPRETER_VERSION = (1, 4, 2)

# --- B1 maximo de la bobina de transmision ----------------------------------
# PENDIENTE DE CONFIRMAR.  PyPulseq no modela un limite de B1 (pp.Opts no tiene
# campo para el), asi que si se excede no salta aqui: lo rechaza el interprete
# del scanner, ya con el paciente en la camilla.  Por eso la auditoria de abajo.
#
# Como averiguarlo en un Siemens (ver docs/qiss_timing.md, seccion 10):
#   B1_max [uT] = 11.74 * V_amplificador_max / V_referencia
# porque la "reference voltage" de Siemens se define como la que produce un
# pulso rectangular de 180 grados en 1 ms, que son exactamente 11.74 uT.
MAX_B1_UT = None                # poner el valor real en uT; None = solo informa

# Los WURST degradan de forma SUAVE al bajar el B1 (no hay un piso abrupto): la
# eficiencia de inversion en el centro cae a 0.95 en torno a 0.78x del B1 de diseno
# (barrido en sim_inversion_profile.py).  Por debajo de ese valor el nulo del fondo
# se corre y hay que rediseñar el pulso, no solo aceptar la perdida.
B1_FRACCION_MINIMA_ADIABATICO = 0.78

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
PARTITION_THICKNESS = 1.3e-3    # m, igual que el paper (19.5 mm x 1.2 / 18)
N_KZ_GRID = 20                  # pasos en la grilla sobremuestreada -> FOV_z = 26 mm.
                                # El paper usa 18 (20 % de oversampling), pero su FLASH
                                # admite un RF largo y selectivo.  En bSSFP el RF tiene
                                # que ser corto y su perfil es mas blando: con 18 pasos
                                # el tejido excitado fuera del FOV_z se repliega dentro
                                # del slab (ver EXC_* y sim_inversion_profile.py).
PARTIAL_FOURIER_KZ = 6 / 8      # -> se adquieren 15 de 20 pasos

# --- Timing (derivado en docs/qiss_timing.md, no numeros magicos)
READOUT_BW_TARGET = 500         # Hz/px.  El TR bSSFP se DERIVA de este valor.
                                # Compromiso medido en sim_bssfp_contrast.py:
                                #   356 Hz/px -> TR 7.68 ms, arteria:grasa  9:1
                                #   500 Hz/px -> TR 7.00 ms, arteria:grasa 21:1, -14 % SNR
                                #   700 Hz/px -> TR 6.64 ms, arteria:grasa 35:1, -27 % SNR
                                # Un TR mas corto acerca la grasa (79.6 Hz) al nulo del
                                # stopband, en +-1/(2*TR).  500 es el punto medio.
TR_QISS = 700e-3                # s, periodo de repeticion del modulo QISS completo.
                                # 700 y no 600 ms: con la eficiencia real de inversion
                                # (0.88, no 1) el nulo del musculo se adelanta, y a
                                # 600 ms el QI ya no deja entrar los 19.5 mm de sangre
                                # vertebral.  Ver docs/qiss_timing.md.
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

# --- Excitacion bSSFP
#
# Parametros elegidos por simulacion del perfil (sim_inversion_profile.py).  El sinc
# de 600 us / TBW 3 de la primera version dejaba el flip en 44 grados en el borde del
# slab y todavia 17 grados a 13.8 mm del centro.
#   * "slice_thickness" en PyPulseq es el ancho a media altura: para tener flip util
#     en los 19.5 mm que se conservan hay que excitar mas grueso (Edelman MRM 2020 usa
#     un multiplicador de 1.5 por la misma razon).
#   * TBW 6 lleva el flip por debajo de 3 grados a 15.6 mm del centro, que con
#     FOV_z = 26 mm queda fuera de la zona que se repliega sobre el slab conservado.
#   * 2 ms de duracion piden 19.6 uT de pico.  El B1 escala como 1/duracion: si el
#     scanner da mas, acortar el pulso devuelve tiempo al readout.
EXC_DURATION = 2.0e-3                 # s
EXC_TBW = 6                           # producto tiempo-ancho de banda
EXC_APODIZATION = 0.5
EXC_THICKNESS = 1.2 * SLAB_THICKNESS  # m, ancho a media altura del perfil excitado

# --- Inversion de fondo (WURST)
#
# 6 kHz / 24 ms mantiene el B1 de pico en 9.4 uT.  El espesor nominal es MAYOR que el
# slab a proposito: el WURST solo invierte bien (Mz <= -0.8) en el 81 % central de su
# ancho nominal.  Con espesor = slab, las particiones de borde quedaban sin suprimir
# (bandas cada 15 mm).  Con 24.1 mm la zona bien invertida cubre los 19.5 mm y la
# inversion se pasa 2.35 mm por debajo del slab: 6 % del bolo carotideo y 11 % del
# vertebral llegan atenuados.  Es el precio de no tener FOCI.
INV_PULSE_TYPE = "wurst"
INV_DURATION = 24e-3                  # s
INV_BANDWIDTH = 6000                  # Hz, ancho de barrido
INV_SLAB_THICKNESS = 24.1e-3          # m

# Eficiencia de inversion en musculo (T2 = 55 ms relaja durante los 24 ms del pulso).
# Medida por Bloch en el centro del slab; verify_qiss3d.py la vuelve a medir y avisa
# si deja de coincidir.  Entra en el calculo del nulo.
INV_EFFICIENCY_MUSCLE = 0.88

# --- Saturacion venosa
#
# SATURACION y no inversion.  La sangre yugular fluye: recorre 7-14 cm por TR_QISS,
# asi que cada disparo ve UN solo pulso partiendo de M0.  Tras una inversion unica
# la sangre tiene Mz = -0.66 en el centro de k, que en una imagen de magnitud es casi
# tan brillante como la arteria.  Tras una saturacion tiene Mz = +0.16.  (Los papers
# invierten, pero con QI de 500-700 ms y T1 mas largos; aqui el QI es de ~200 ms.)
#
# La banda solo tiene que cubrir lo que la sangre venosa alcanza a BAJAR durante el
# QI: 41 mm a 20 cm/s, 62 mm a 30 cm/s.  Los 100 mm del paper responden a su QI de
# 583 ms.  Una banda mas fina permite un gradiente selector mas fuerte y un borde
# mas abrupto en milimetros.
VSAT_THICKNESS = 60e-3                # m
VSAT_DURATION = 10e-3                 # s
VSAT_TBW = 16                         # borde: Mz vuelve a 0.95 a 0.9 mm fuera de la banda
VSAT_GAP = 2e-3                       # m entre el tope del slab y el borde de la banda.
                                      # Si la banda tocara el slab, el fondo saturado
                                      # ahi no pasaria por el nulo y saldria brillante.
SPOILER_CYCLES_PER_VOXEL = 4          # ciclos de fase por voxel

# --- Saturacion grasa espectral
#
# Aporta MENOS de lo que se creyo al principio.  En el eco de kz=0 (TR #11 del
# tren) la grasa ya casi olvido la saturacion y su senal la fija sobre todo la
# cercania al nulo del stopband; el fat-sat baja la grasa en los primeros TRs
# (0.116 -> 0.074 en el primero) y poco en kz=0 (0.036 -> 0.031).  Se conserva
# porque cuesta 19 ms de un intervalo que estaria muerto y 0.5 uT de B1.
#
# Duracion elegida por simulacion del perfil espectral: a 0.55T solo hay 79.6 Hz entre
# grasa y agua.  15 ms con TBW 1.0 deja el agua en +0.996 y la grasa en +0.070 en
# resonancia, y aguanta +-20 Hz de B0 con el agua en >= 0.939.
FATSAT_DURATION = 15e-3               # s
FATSAT_TBW = 1.0

# Direccion de flujo: en el cuello la sangre arterial sube (caudal -> craneal) y la
# venosa baja, asi que la banda venosa va POR ENCIMA del slab.
# A VERIFICAR EN EL SCANNER: que +z logico de Pulseq sea craneal depende de la posicion
# del paciente registrada en el protocolo.  Si esta al reves la banda cae DEBAJO del
# slab y satura el bolo arterial: no hay angiograma.  Probar con un fantoma marcado.
Z_SUPERIOR = +1


# =============================================================================
# 3. BLOQUES REUTILIZABLES
# =============================================================================

def make_excitation(sys, flip_angle_deg):
    """Sinc selectivo de slab para el tren bSSFP (y su alpha/2, misma forma)."""
    return pp.make_sinc_pulse(
        flip_angle=np.deg2rad(flip_angle_deg),
        duration=EXC_DURATION,
        slice_thickness=EXC_THICKNESS,
        apodization=EXC_APODIZATION,
        time_bw_product=EXC_TBW,
        system=sys,
        return_gz=True,
    )


def make_adiabatic_inversion(sys, thickness, z_center, duration=INV_DURATION):
    """Inversion adiabatica selectiva de slab (sustituye al FOCI de los papers).

    El perfil se verifico por integracion de Bloch (sim_inversion_profile.py).
    El SAR de un adiabatico limita en campo alto; a 0.55T escala con B0^2 y no.
    """
    kwargs = dict(pulse_type=INV_PULSE_TYPE, duration=duration,
                  slice_thickness=thickness, system=sys, return_gz=True,
                  use="inversion")
    if INV_PULSE_TYPE == "wurst":
        kwargs["bandwidth"] = INV_BANDWIDTH
    rf, gz, _ = pp.make_adiabatic_pulse(**kwargs)
    # df = G * z traslada la banda a su posicion axial (G en Hz/m).
    rf.freq_offset = gz.amplitude * z_center
    return rf, gz


def make_venous_saturation(sys, z_center):
    """Saturacion de 90 grados selectiva en espacio para la banda venosa."""
    rf, gz, _ = pp.make_sinc_pulse(
        flip_angle=np.pi / 2, duration=VSAT_DURATION,
        slice_thickness=VSAT_THICKNESS, apodization=0.5,
        time_bw_product=VSAT_TBW, system=sys, return_gz=True, use="saturation",
    )
    rf.freq_offset = gz.amplitude * z_center
    return rf, gz


def make_fatsat(sys, duration=FATSAT_DURATION, tbw=FATSAT_TBW):
    """Pulso de 90 grados selectivo en FRECUENCIA (no en espacio) centrado en la grasa.

    Sin gradiente durante el RF: no puede saturar la sangre entrante por su posicion,
    que es la objecion por la que Koktzoglou evita el fat-sat a 3T.
    """
    rf = pp.make_gauss_pulse(
        flip_angle=np.pi / 2, duration=duration, time_bw_product=tbw,
        freq_offset=-DF_FAT,      # la grasa resuena por debajo del agua
        system=sys, use="saturation",
    )
    return rf[0] if isinstance(rf, tuple) else rf


def make_spoiler(sys, voxel_size_m, n_cycles=SPOILER_CYCLES_PER_VOXEL, axis="z"):
    """Spoiler que dispersa n_cycles ciclos completos de fase a lo largo de un voxel.

    4 ciclos es el punto de partida habitual.  El area va en 1/m.
    """
    return pp.make_trapezoid(channel=axis, area=n_cycles / voxel_size_m, system=sys)


def _raster(t, sys, mode=np.round):
    return mode(t / sys.grad_raster_time + (1e-9 if mode is np.round else 0)) \
        * sys.grad_raster_time


def make_bssfp_core(sys, bw_target=None):
    """Eventos invariantes del TR bSSFP y su reparto exacto de tiempos.

    Estructura de cada TR (todos los TR duran exactamente lo mismo):

        A  RF + selector de slab                              t_rf
        B  rephase de slab + kz   ||  prefase de lectura      t_b
        C  lectura + ADC                                      t_read
        D  rebobinado kz + prefase de slab || rebobinado      t_d

    B y D tienen duracion FIJA, dimensionada para el kz mas exigente; asi el TR no
    depende del indice de particion.  bSSFP exige momento neto cero por TR en los
    tres ejes y, ademas, momento cero de RF a RF *en el centro de cada pulso*.

    El TR no se impone: sale de sumar el RF, los rebobinados minimos y una lectura
    del ancho de banda pedido.  (La primera version imponia 6.28 ms para poner la
    grasa en el nulo del stopband; con un RF de perfil aceptable ese TR ya no cabe
    en estos gradientes, y el fat-sat hace innecesario perseguirlo.)
    """
    bw_target = READOUT_BW_TARGET if bw_target is None else bw_target
    n_read = int(round(FOV_XY / RES_XY / 4)) * 4     # Siemens: multiplo de 4
    k_max = 1 / (2 * RES_XY)

    rf, gz, gz_rephase = make_excitation(sys, FLIP_ANGLE_DEG)
    rf_half, gz_half, _ = make_excitation(sys, FLIP_ANGLE_DEG / 2)

    # Centro del RF dentro de su bloque.  No siempre cae a mitad: si la rampa del
    # selector es mas corta que el dead time del RF, PyPulseq retrasa el gradiente.
    t_rf = pp.calc_duration(rf, gz)
    t_pre = rf.delay + pp.calc_rf_center(rf)[0]      # inicio de bloque -> centro del RF
    t_post = t_rf - t_pre                            # centro del RF -> fin de bloque

    # Lobulo z mas exigente: rephase de slab + |kz| maximo
    area_z_max = abs(gz_rephase.area) + np.abs(KZ_INDICES).max() * DELTA_KZ
    t_z = pp.calc_duration(pp.make_trapezoid(channel="z", area=area_z_max, system=sys))

    # A 0.55T interesa una lectura larga: menos ancho de banda = mas SNR.
    # Dwell en multiplos de 0.5 us: con n_read multiplo de 4 la ventana del ADC dura
    # un numero PAR de microsegundos y se puede centrar en la meseta con un retardo
    # entero en us.  (Con 0.1 us de resolucion el centrado dejaba un retardo de
    # 260.4 us, que el check_timing del pypulseq 1.4.2 publicado rechaza.)
    dwell = np.round(1 / (bw_target * n_read) / 0.5e-6) * 0.5e-6
    flat = _raster(n_read * dwell, sys, np.ceil)
    for _ in range(3):
        g_read = pp.make_trapezoid(channel="x", flat_area=2 * k_max, system=sys,
                                   flat_time=flat)
        g_read_pre = pp.make_trapezoid(channel="x", area=-g_read.area / 2, system=sys)
        t_pw = max(t_z, pp.calc_duration(g_read_pre))
        # Para que TE = TR/2 exacto hace falta (t_post + t_b) = (t_d + t_pre)
        t_b = t_pw + max(t_pre - t_post, 0.0)
        t_d = t_pw + max(t_post - t_pre, 0.0)
        t_read = pp.calc_duration(g_read)
        tr = t_rf + t_b + t_read + t_d
        # TR/2 tambien tiene que caer en el raster (lo usan el eco y el alpha/2):
        # si el TR sale impar en unidades de raster se alarga la lectura un paso.
        if int(round(tr / sys.grad_raster_time)) % 2 == 0:
            break
        flat += sys.grad_raster_time
    tr_real = tr

    adc_duration = n_read * dwell
    adc_offset = (g_read.flat_time - adc_duration) / 2     # entero en us, ver arriba
    adc = pp.make_adc(num_samples=n_read, dwell=dwell,
                      delay=g_read.rise_time + adc_offset, system=sys)

    # Entre el alpha/2 y el primer alpha: rephase del alpha/2 MAS prefase del alpha.
    # Son -gz.area en total.  (Poner solo el rephase, como en la primera version,
    # dejaba 1.9 ciclos de fase a traves del slab y arruinaba la catalizacion.)
    g_half_gap = pp.make_trapezoid(channel="z", area=-gz.area, system=sys)
    t_half_gap = tr / 2 - t_post - t_pre
    if pp.calc_duration(g_half_gap) > t_half_gap + 1e-9:
        raise ValueError("El lobulo entre alpha/2 y alpha no cabe en TR/2")

    return SimpleNamespace(
        rf=rf, gz=gz, gz_rephase=gz_rephase, rf_half=rf_half, gz_half=gz_half,
        g_half_gap=g_half_gap, t_half_gap=_raster(t_half_gap, sys),
        g_read=g_read, g_read_pre=g_read_pre, adc=adc,
        n_read=n_read, k_max=k_max, bw_per_pixel=1 / adc_duration,
        t_rf=t_rf, t_pre=t_pre, t_post=t_post, t_b=t_b, t_d=t_d, t_read=t_read,
        tr=tr_real, te=t_post + t_b + t_read / 2,
    )


def rotate_in_plane(grad_proto, theta, sys):
    """Proyecta un gradiente prototipo del eje x sobre (x, y) para el angulo theta."""
    gx = pp.scale_grad(grad_proto, np.cos(theta))
    gy = pp.scale_grad(pp.make_trapezoid(
        channel="y", amplitude=grad_proto.amplitude,
        flat_time=getattr(grad_proto, "flat_time", 0),
        rise_time=grad_proto.rise_time, fall_time=grad_proto.fall_time, system=sys,
    ), np.sin(theta))
    return gx, gy


# =============================================================================
# 4. TIMING DERIVADO
# =============================================================================

def null_ti_steady_state(t1, tr_qiss, efficiency=1.0):
    """TI de nulo con recuperacion incompleta e inversion imperfecta.

    Con eficiencia f (la inversion lleva Mz a -f*Mz), la Mz justo antes de cada
    inversion en estado estacionario es
        M = (1 - E) / (1 + f*E),   E = exp(-TR/T1),
    y cruza cero en  TI = T1 * ln(1 + f*M).
    Con f = 1 y TR >> T1 tiende al T1*ln2 de libro, que aqui daria un error de 100 ms.
    """
    e = np.exp(-tr_qiss / t1)
    m = (1 - e) / (1 + efficiency * e)
    return t1 * np.log(1 + efficiency * m)


N_KZ_ACQUIRED = int(round(N_KZ_GRID * PARTIAL_FOURIER_KZ))     # 15 de 20
FOV_Z = N_KZ_GRID * PARTITION_THICKNESS                        # 26 mm
SLICE_OVERSAMPLING = FOV_Z / SLAB_THICKNESS                    # 1.33
DELTA_KZ = 1 / FOV_Z                                           # 1/m

# Partial Fourier: lado negativo completo, positivo recortado (conjugado hermitico).
KZ_INDICES = np.arange(-N_KZ_GRID // 2, -N_KZ_GRID // 2 + N_KZ_ACQUIRED)
IDX_KZ0 = int(np.where(KZ_INDICES == 0)[0][0])   # TR que pasa por kz = 0 (NO es el central)

# QI: del CENTRO de la inversion de fondo al ECO del TR que muestrea kz = 0.
QI = null_ti_steady_state(T1_MUSCLE, TR_QISS, INV_EFFICIENCY_MUSCLE)

# TR del tren, derivado (no impuesto) del RF, los rebobinados y el ancho de banda.
TR_BSSFP = make_bssfp_core(system).tr
T_ACQ_WINDOW = N_KZ_ACQUIRED * TR_BSSFP


# =============================================================================
# 5. ENSAMBLADO
# =============================================================================

def build_sequence(n_slabs=N_SLABS, n_views=N_VIEWS, verbose=True):
    seq = pp.Sequence(system=system)
    core = make_bssfp_core(system)

    spoiler = make_spoiler(system, PARTITION_THICKNESS)
    rf_fatsat = make_fatsat(system)
    t_spoil = pp.calc_duration(spoiler)
    t_fatsat = pp.calc_duration(rf_fatsat) + t_spoil

    slab_step = SLAB_THICKNESS - SLAB_OVERLAP
    z_centers = (np.arange(n_slabs) - (n_slabs - 1) / 2) * slab_step
    timing = None

    for i_slab, z_c in enumerate(z_centers):
        freq_slab = core.gz.amplitude * z_c
        core.rf.freq_offset = freq_slab
        core.rf_half.freq_offset = freq_slab

        z_venous = z_c + Z_SUPERIOR * (SLAB_THICKNESS / 2 + VSAT_GAP + VSAT_THICKNESS / 2)

        # Dependen del slab y no del angulo radial: se construyen una vez por slab.
        rf_inv, gz_inv = make_adiabatic_inversion(system, INV_SLAB_THICKNESS, z_c)
        rf_vsat, gz_vsat = make_venous_saturation(system, z_venous)
        t_inv = pp.calc_duration(rf_inv, gz_inv)
        t_vsat = pp.calc_duration(rf_vsat, gz_vsat)

        # El centro del slab cruza la resonancia a mitad del barrido del WURST.
        t_inv_center = rf_inv.delay + INV_DURATION / 2

        # eco(kz=0) = inicio del alpha/2 + t_pre + (IDX_KZ0 + 1) * TR
        t_half_start = _raster(t_inv_center + QI - (IDX_KZ0 + 1) * core.tr - core.t_pre,
                               system)
        t_qi_delay = t_half_start - (t_inv + t_vsat + t_spoil + t_fatsat)
        if t_qi_delay < 0:
            raise ValueError(f"QI insuficiente: faltan {-t_qi_delay*1e3:.1f} ms para "
                             "alojar la preparacion antes del tren")
        t_shot = t_half_start + core.t_rf + core.t_half_gap + N_KZ_ACQUIRED * core.tr
        t_fill = _raster(TR_QISS - t_shot, system)
        if t_fill < 0:
            raise ValueError(f"El disparo ({t_shot*1e3:.1f} ms) no cabe en TR_QISS")

        if timing is None:
            timing = SimpleNamespace(
                t_inv_center=t_inv_center, t_vsat_center=t_inv + rf_vsat.delay
                + pp.calc_rf_center(rf_vsat)[0],
                t_fatsat_start=t_half_start - t_fatsat, t_half_start=t_half_start,
                t_first_alpha=t_half_start + core.t_pre + core.tr / 2,
                t_echo_kz0=t_half_start + core.t_pre + (IDX_KZ0 + 1) * core.tr,
                t_shot=t_shot)

        for i_view in range(n_views):
            theta = view_angle(i_view, n_views)

            # ---------------- preparacion ----------------
            seq.add_block(rf_inv, gz_inv)             # nulo del fondo en el slab
            seq.add_block(rf_vsat, gz_vsat)           # sangre venosa por encima
            seq.add_block(spoiler)
            seq.add_block(pp.make_delay(t_qi_delay))  # intervalo quiescente
            seq.add_block(rf_fatsat)                  # pegado al tren: la grasa no
            seq.add_block(spoiler)                    # tiene tiempo de recuperarse

            # ---------------- catalizacion alpha/2 (Le Roux) ----------------
            core.rf_half.phase_offset = np.pi         # opuesta al primer alpha
            seq.add_block(core.rf_half, core.gz_half)
            seq.add_block(core.g_half_gap, pp.make_delay(core.t_half_gap))

            # ---------------- tren bSSFP: inner slice loop ----------------
            gx_read, gy_read = rotate_in_plane(core.g_read, theta, system)
            gx_pre, gy_pre = rotate_in_plane(core.g_read_pre, theta, system)

            for i_tr, kz_idx in enumerate(KZ_INDICES):
                phase = np.pi * (i_tr % 2)            # alternancia 0 / 180
                core.rf.phase_offset = phase
                core.adc.phase_offset = phase

                seq.add_block(core.rf, core.gz)

                # kz y rephase de slab comparten eje: un solo trapecio, concurrente
                # con la prefase de lectura.  Es lo que mantiene corto el TR.
                gz_b = pp.make_trapezoid(channel="z", system=system,
                                         area=core.gz_rephase.area + kz_idx * DELTA_KZ)
                seq.add_block(gz_b, gx_pre, gy_pre, pp.make_delay(core.t_b))

                seq.add_block(gx_read, gy_read, core.adc)

                gz_d = pp.make_trapezoid(channel="z", system=system,
                                         area=core.gz_rephase.area - kz_idx * DELTA_KZ)
                seq.add_block(gz_d, gx_pre, gy_pre, pp.make_delay(core.t_d))

            seq.add_block(pp.make_delay(t_fill))      # completa el TR_QISS

        if verbose:
            print(f"  slab {i_slab+1}/{n_slabs} en z = {z_c*1e3:+7.1f} mm  listo")

    core.timing = timing
    return seq, core


# =============================================================================
# 6. INFORME Y VERIFICACION
# =============================================================================

def pulse_table(core):
    """(nombre, evento RF) de todos los pulsos de la secuencia, para auditar el B1."""
    rf_inv, _ = make_adiabatic_inversion(system, INV_SLAB_THICKNESS, 0.0)
    rf_vsat, _ = make_venous_saturation(system, 0.0)
    return [("excitacion bSSFP", core.rf), ("inversion de fondo WURST", rf_inv),
            ("saturacion venosa", rf_vsat), ("fat-sat espectral", make_fatsat(system))]


def report(core, n_slabs, n_views):
    tm = core.timing
    print("\n" + "=" * 74)
    print("  3D tsSOS-QISS bSSFP  |  MAGNETOM Free.Max 0.55T")
    print("=" * 74)
    print(f"  Larmor / grasa-agua         : {F_LARMOR/1e6:.2f} MHz / {DF_FAT:.1f} Hz")
    print(f"  TR / TE bSSFP               : {core.tr*1e3:.2f} / {core.te*1e3:.3f} ms "
          f"(nulos del stopband en +-{1/(2*core.tr):.1f} Hz)")
    print(f"  Flip angle                  : {FLIP_ANGLE_DEG} deg")
    print(f"  Excitacion                  : sinc {EXC_DURATION*1e3:.1f} ms, TBW {EXC_TBW}, "
          f"FWHM {EXC_THICKNESS*1e3:.1f} mm")
    print(f"  Ancho de banda receptor     : {core.bw_per_pixel:.0f} Hz/px "
          f"({core.n_read} muestras)")
    print("-" * 74)
    print(f"  QISS TR                     : {TR_QISS*1e3:.0f} ms  (ungated)")
    print(f"  QI (centro inversion -> eco de kz=0): {QI*1e3:.1f} ms  "
          f"[medido en la secuencia: {(tm.t_echo_kz0 - tm.t_inv_center)*1e3:.1f} ms]")
    print(f"     eficiencia de inversion usada    : {INV_EFFICIENCY_MUSCLE:.2f}")
    print(f"     T1*ln2 de libro daria            : {T1_MUSCLE*np.log(2)*1e3:.1f} ms")
    print(f"  kz = 0 en el TR             : #{IDX_KZ0+1} de {N_KZ_ACQUIRED}")
    print(f"  Ventana de adquisicion      : {N_KZ_ACQUIRED*core.tr*1e3:.1f} ms")
    print(f"  Cronologia del disparo (ms desde su inicio):")
    print(f"     centro inversion fondo {tm.t_inv_center*1e3:6.1f} | centro sat. venosa "
          f"{tm.t_vsat_center*1e3:6.1f} | fat-sat {tm.t_fatsat_start*1e3:6.1f}")
    print(f"     primer alpha {tm.t_first_alpha*1e3:6.1f} | eco kz=0 {tm.t_echo_kz0*1e3:6.1f} "
          f"| fin del tren {tm.t_shot*1e3:6.1f}")
    print("-" * 74)
    cobertura = (n_slabs - 1) * (SLAB_THICKNESS - SLAB_OVERLAP) + SLAB_THICKNESS
    print(f"  Slabs                       : {n_slabs} x {SLAB_THICKNESS*1e3:.1f} mm, "
          f"solape {SLAB_OVERLAP*1e3:.1f} mm -> {cobertura*1e3:.1f} mm")
    print(f"  Particiones                 : {N_KZ_ACQUIRED} de {N_KZ_GRID} "
          f"({N_KZ_ACQUIRED/N_KZ_GRID:.2f}), FOV_z {FOV_Z*1e3:.1f} mm "
          f"(oversampling {(SLICE_OVERSAMPLING-1)*100:.0f} %)")
    print(f"  Resolucion adquirida        : {RES_XY*1e3:.2f} x {RES_XY*1e3:.2f} x "
          f"{PARTITION_THICKNESS*1e3:.2f} mm")
    n_nyq = int(np.ceil(np.pi / 2 * core.n_read))
    print(f"  Vistas radiales             : {n_views} ({VIEW_ORDER}; Nyquist {n_nyq}, "
          f"submuestreo {n_nyq/n_views:.1f}x)")
    t_total = n_slabs * n_views * TR_QISS
    print(f"  Tiempo de examen            : {int(t_total//60)}:{int(t_total%60):02d}")
    print("-" * 74)
    t_inflow = tm.t_echo_kz0 - tm.t_inv_center
    print("  Inflow arterial hasta el eco de kz=0 (velocidades medias de Koktzoglou 2020):")
    for vaso, v in [("carotida 20 cm/s", 0.20), ("vertebral 10 cm/s", 0.10)]:
        d = v * t_inflow
        print(f"     {vaso:18s}: {d*1e3:5.1f} mm = {d/SLAB_THICKNESS:.2f} x el slab")
    print("  El contraste por tejido NO se estima aqui: ver sim_bssfp_contrast.py.")
    print("=" * 74 + "\n")


def audit_b1(core):
    print("Auditoria de B1 de pico por pulso:")
    peak = 0.0
    for nombre, r in pulse_table(core):
        b1 = np.abs(r.signal).max() / GAMMA_BAR * 1e6
        peak = max(peak, b1)
        aviso = ""
        if MAX_B1_UT is not None and b1 > MAX_B1_UT:
            aviso = f"   EXCEDE el limite de {MAX_B1_UT} uT"
        print(f"  {nombre:26s}: {b1:6.2f} uT{aviso}")
    print(f"  {'PICO DE LA SECUENCIA':26s}: {peak:6.2f} uT")
    if MAX_B1_UT is None:
        print("  -> MAX_B1_UT sin definir.  PyPulseq no comprueba el B1: si el scanner no")
        print("     llega al pico, lo rechaza el interprete.  La excitacion se puede alargar")
        print("     (B1 ~ 1/duracion); los WURST degradan suavemente por debajo de "
              f"{B1_FRACCION_MINIMA_ADIABATICO:.2f}x su B1.")
    return peak


if __name__ == "__main__":
    import sys as _sys

    quick = "--quick" in _sys.argv
    n_slabs = 2 if quick else N_SLABS
    n_views = 4 if quick else N_VIEWS
    if quick:
        print(">> modo rapido: 2 slabs x 4 vistas (prueba de carga en el scanner)\n")

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
    seq.set_definition("TR_bSSFP_ms", core.tr * 1e3)
    seq.set_definition("TR_QISS_ms", TR_QISS * 1e3)
    seq.set_definition("QI_ms", QI * 1e3)
    seq.set_definition("FlipAngle_deg", FLIP_ANGLE_DEG)
    seq.set_definition("RadialViews", n_views)
    seq.set_definition("ViewOrder", VIEW_ORDER)
    seq.set_definition("Slabs", n_slabs)
    seq.set_definition("KzGrid", N_KZ_GRID)
    seq.set_definition("KzAcquired", N_KZ_ACQUIRED)
    seq.set_definition("B0_T", B0)

    out = "qiss3d_tssos_bssfp_055T_quick.seq" if quick else "qiss3d_tssos_bssfp_055T.seq"
    seq.write(out)
    print(f"\nExportado: {out}")

    print("\nCompatibilidad con el interprete Pulseq "
          f"{'.'.join(map(str, SCANNER_INTERPRETER_VERSION))} del scanner:")
    with open(out) as fh:
        head = fh.read(400)
    ver = tuple(int(re.search(rf"{k}\s+(\d+)", head).group(1))
                for k in ("major", "minor", "revision"))
    print(f"  Version declarada en el .seq : {'.'.join(map(str, ver))}"
          f"  {'(coincide)' if ver == SCANNER_INTERPRETER_VERSION else '(NO COINCIDE)'}")
    n_blocks = len(seq.block_events)
    print(f"  Bloques                      : {n_blocks}")
    print(f"  Tamano de archivo            : {os.path.getsize(out)/1e6:.2f} MB")
    print(f"  Muestras por ADC             : {core.n_read} "
          f"({'multiplo de 4, OK' if core.n_read % 4 == 0 else 'NO es multiplo de 4'})")
    dwell = core.adc.dwell
    print(f"  Dwell del ADC                : {dwell*1e9:.0f} ns "
          f"({'en raster de 100 ns' if round(dwell*1e9) % 100 == 0 else 'FUERA de raster'})")
    print("  Solo RF arbitrario, trapecios, ADC y delays: sin labels, triggers ni extensions.")
    print(f"  A VERIFICAR EN EL SCANNER: que el interprete cargue {n_blocks} bloques.")
    print("  Probar antes con --quick.\n")

    audit_b1(core)
