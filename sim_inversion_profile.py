"""Perfil de slab de los pulsos adiabaticos de inversion, por integracion de Bloch.

Responde la unica pregunta que decide si el hyperbolic secant sirve como sustituto
del FOCI de los papers:

    ¿Cuantos milimetros por DEBAJO del slab llega realmente la inversion de fondo?

Importa porque QISS depende de que la sangre arterial que sube desde abajo no haya
visto el pulso.  Si el borde inferior es blando, la cabecera del bolo entrante llega
invertida y por tanto oscura, que es justo lo que Edelman (MRM 2020) evita usando
FOCI en vez de un sinc: "la region de supresion del FOCI no se extiende mas alla del
borde del slab y tiene un perfil mucho mas abrupto".

Metodo: se toma la forma de onda COMPLEJA real del pulso tal como la genera PyPulseq
y se integra la ecuacion de Bloch en el sistema rotante para una columna de isocromas
a lo largo de z.  No hay aproximacion de angulo pequeno ni linealizacion: se rota el
vector de magnetizacion en torno al campo efectivo en cada paso de tiempo, que es
exacto salvo por la discretizacion temporal (paso = raster de RF, 1 us).

Equivalente al "nivel 1" que se haria en KomaMRI; se hace aqui porque la politica de
egress de la sesion bloquea los hosts de descarga de Julia.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pypulseq as pp

import qiss3d_tssos_bssfp_055T as q


def bloch_profile(rf, gz_amplitude, z_positions, t1, t2, dt=None):
    """Integra Bloch a lo largo del pulso y devuelve Mz final para cada posicion z.

    rf.signal viene en Hz (es gamma*B1/2pi), complejo: la modulacion de frecuencia
    del adiabatico esta codificada en la fase.  El desfase por posicion es
    Delta_omega = 2*pi * G * z, con G en Hz/m.
    """
    t = np.asarray(rf.t)
    sig = np.asarray(rf.signal)              # Hz, complejo
    if dt is None:
        dt = float(np.min(np.diff(t)))

    # Remuestreo a paso uniforme para integrar
    t_uni = np.arange(t[0], t[-1], dt)
    sig_u = np.interp(t_uni, t, sig.real) + 1j * np.interp(t_uni, t, sig.imag)

    w1 = 2 * np.pi * sig_u                   # rad/s, complejo
    w1x, w1y = w1.real, w1.imag

    e1, e2 = np.exp(-dt / t1), np.exp(-dt / t2)

    mz_out = np.empty_like(z_positions)
    for i, z in enumerate(z_positions):
        dw = 2 * np.pi * gz_amplitude * z    # rad/s, desfase por posicion
        m = np.array([0.0, 0.0, 1.0])        # equilibrio

        for k in range(len(t_uni)):
            # Campo efectivo en el sistema rotante
            wx, wy, wz = w1x[k], w1y[k], dw
            wnorm = np.sqrt(wx * wx + wy * wy + wz * wz)
            if wnorm > 0:
                ax, ay, az = wx / wnorm, wy / wnorm, wz / wnorm
                th = -wnorm * dt             # rotacion levogira en torno al eje
                c, s = np.cos(th), np.sin(th)
                dot = ax * m[0] + ay * m[1] + az * m[2]
                # Formula de rotacion de Rodrigues
                m = (m * c
                     + np.array([ay * m[2] - az * m[1],
                                 az * m[0] - ax * m[2],
                                 ax * m[1] - ay * m[0]]) * s
                     + np.array([ax, ay, az]) * dot * (1 - c))
            # Relajacion
            m[0] *= e2
            m[1] *= e2
            m[2] = 1 + (m[2] - 1) * e1

        mz_out[i] = m[2]
    return mz_out


def transition_width(z, mz, edge="inferior"):
    """Ancho de transicion 10-90% en el borde indicado, en mm.

    Se define entre Mz = -0.8 (invertido, 90% del camino) y Mz = +0.8 (intacto).
    """
    if edge == "inferior":
        mask = z < 0
        zz, mm = z[mask][::-1], mz[mask][::-1]   # del centro hacia afuera
    else:
        mask = z > 0
        zz, mm = z[mask], mz[mask]

    try:
        z_inv = zz[np.where(mm <= -0.8)[0][-1]]   # ultimo punto bien invertido
        z_free = zz[np.where(mm >= 0.8)[0][0]]    # primer punto intacto
    except IndexError:
        return np.nan, np.nan, np.nan
    return abs(z_free - z_inv) * 1e3, z_inv * 1e3, z_free * 1e3


# =============================================================================
print("=" * 74)
print("  PERFIL DE INVERSION — pulso adiabatico vs FOCI de los papers")
print("=" * 74)

sys_ = q.system
T1B, T2B = q.T1_BLOOD, q.T2_BLOOD

# ---- Inversion de fondo: coextensiva con el slab de 19.5 mm
rf_bg, gz_bg = q.make_adiabatic_inversion(sys_, q.INV_SLAB_THICKNESS, 0.0)
rf_bg.freq_offset = 0.0     # perfil centrado en z=0 para medirlo

b1_peak_ut = np.abs(rf_bg.signal).max() / q.GAMMA_BAR * 1e6
print(f"\nPulso: {q.INV_PULSE_TYPE}, duracion {q.INV_DURATION*1e3:.1f} ms")
print(f"  Gradiente selector : {gz_bg.amplitude/q.GAMMA_BAR*1e3:.2f} mT/m")
print(f"  B1 de pico exigido : {b1_peak_ut:.2f} uT")
print(f"  Ancho de banda     : {gz_bg.amplitude*q.INV_SLAB_THICKNESS:.0f} Hz "
      f"sobre el slab de {q.INV_SLAB_THICKNESS*1e3:.1f} mm")

z = np.linspace(-40e-3, 40e-3, 321)          # +-40 mm, paso 0.25 mm
print("\nIntegrando Bloch sobre 321 posiciones...")
mz_bg = bloch_profile(rf_bg, gz_bg.amplitude, z, T1B, T2B)

w_inf, z_inv_i, z_free_i = transition_width(z, mz_bg, "inferior")
w_sup, z_inv_s, z_free_s = transition_width(z, mz_bg, "superior")

print(f"\n  Borde INFERIOR (por donde entra la sangre arterial):")
print(f"    ultimo punto invertido (Mz<=-0.8) : {z_inv_i:+.2f} mm")
print(f"    primer punto intacto   (Mz>=+0.8) : {z_free_i:+.2f} mm")
print(f"    ancho de transicion 10-90%        : {w_inf:.2f} mm")
print(f"  Borde superior: transicion {w_sup:.2f} mm")
print(f"  Borde nominal del slab            : {-q.INV_SLAB_THICKNESS/2*1e3:+.2f} mm")
print(f"  Fuga por debajo del borde nominal : "
      f"{abs(z_free_i) - q.INV_SLAB_THICKNESS/2*1e3:.2f} mm")

# ---- Impacto sobre el inflow
print("\n" + "-" * 74)
print("  IMPACTO SOBRE EL INFLOW ARTERIAL")
print("-" * 74)
fuga = abs(z_free_i) - q.INV_SLAB_THICKNESS / 2 * 1e3   # mm de zona contaminada
for vaso, v in [("carotida", 0.20), ("vertebral", 0.10)]:
    inflow_mm = v * q.QI * 1e3
    frac = fuga / inflow_mm * 100
    print(f"  {vaso:10s}: recorre {inflow_mm:5.1f} mm en el QI; "
          f"la zona contaminada son {fuga:.1f} mm -> {frac:.0f} % del bolo")

# ---- Banda venosa
print("\n" + "-" * 74)
print(f"  BANDA VENOSA ({q.INV_VENOUS_THICKNESS*1e3:.0f} mm, {q.INV_VENOUS_GAP*1e3:.0f} mm por encima del slab)")
print("-" * 74)
rf_ven, gz_ven = q.make_adiabatic_inversion(sys_, q.INV_VENOUS_THICKNESS, 0.0)
rf_ven.freq_offset = 0.0
z_ven = np.linspace(-80e-3, 80e-3, 321)
mz_ven = bloch_profile(rf_ven, gz_ven.amplitude, z_ven, T1B, T2B)
w_ven, _, z_free_ven = transition_width(z_ven, mz_ven, "inferior")
print(f"  Gradiente selector : {gz_ven.amplitude/q.GAMMA_BAR*1e3:.3f} mT/m")
print(f"  B1 de pico exigido : {np.abs(rf_ven.signal).max()/q.GAMMA_BAR*1e6:.2f} uT")
print(f"  Transicion en borde inferior : {w_ven:.2f} mm")
invasion = abs(z_free_ven) - q.INV_VENOUS_THICKNESS / 2 * 1e3
margen = q.INV_VENOUS_GAP * 1e3 - invasion
print(f"  Fuga bajo el borde nominal   : {invasion:.2f} mm")
print(f"  Hueco disponible             : {q.INV_VENOUS_GAP*1e3:.1f} mm")
print(f"  -> {'INVADE el slab por ' + f'{-margen:.2f} mm' if margen < 0 else f'no invade el slab (sobran {margen:.2f} mm)'}")

# ---- Grafico
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
axes[0].plot(z * 1e3, mz_bg, lw=1.6)
axes[0].axvspan(-q.INV_SLAB_THICKNESS / 2 * 1e3, q.INV_SLAB_THICKNESS / 2 * 1e3,
                color="tab:blue", alpha=0.12, label="slab nominal 19.5 mm")
axes[0].axhline(0, color="k", lw=0.5)
axes[0].axhline(-0.8, color="r", ls=":", lw=0.8)
axes[0].axhline(0.8, color="g", ls=":", lw=0.8)
axes[0].set_xlabel("z [mm]  (negativo = inferior, de donde viene la sangre)")
axes[0].set_ylabel("$M_z$ tras la inversion")
axes[0].set_title(f"Inversion de fondo ({q.INV_PULSE_TYPE} {q.INV_DURATION*1e3:.0f} ms)\n"
                  f"transicion inferior {w_inf:.1f} mm")
axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)

axes[1].plot(z_ven * 1e3, mz_ven, lw=1.6, color="tab:orange")
axes[1].axvspan(-q.INV_VENOUS_THICKNESS / 2 * 1e3, q.INV_VENOUS_THICKNESS / 2 * 1e3,
                color="tab:orange", alpha=0.12, label=f"banda nominal {q.INV_VENOUS_THICKNESS*1e3:.0f} mm")
axes[1].axhline(0, color="k", lw=0.5)
axes[1].set_xlabel("z [mm] relativo al centro de la banda")
axes[1].set_ylabel("$M_z$ tras la inversion")
axes[1].set_title(f"Banda venosa\ntransicion inferior {w_ven:.1f} mm")
axes[1].legend(fontsize=8); axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("perfil_inversion.png", dpi=130)
print("\nGrafico guardado: perfil_inversion.png")
