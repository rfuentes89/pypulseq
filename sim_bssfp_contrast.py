"""Contraste de la secuencia 3D tsSOS-QISS por simulacion de Bloch del disparo completo.

Simula el ciclo entero y lo repite hasta el estado estacionario:

    inversion de fondo -> recuperacion -> fat-sat -> alpha/2 -> tren bSSFP -> recuperacion

Los instantes salen de la secuencia ensamblada (core.timing), no de valores nominales.
La eficiencia de inversion, el residuo del fat-sat y el flip local salen de integrar
los pulsos reales (bloch_tools.bloch_rf).

Diferencias con la primera version de este script, que daba cifras equivocadas:
  * tomaba la Mz del instante de kz=0 como condicion inicial del tren (55 ms antes);
  * trataba la sangre venosa como estatica, cuando fluye y ve un solo pulso por disparo;
  * suponia inversion perfecta.

Limite del modelo: la sangre que fluye se trata como si entrara al slab ya en el
arranque del tren y viera todos los pulsos.  No hay simulacion de flujo espacial.
"""

import numpy as np

import qiss3d_tssos_bssfp_055T as q
from bloch_tools import bloch_rf, bssfp_train, flip_deg, relax, rot_rf  # noqa: F401


def pulse_effects(core, b0_offset=0.0, z=0.0):
    """Efecto de cada pulso real sobre agua y grasa en la posicion z del slab."""
    rf_inv, gz_inv = q.make_adiabatic_inversion(q.system, q.INV_SLAB_THICKNESS, 0.0)
    rf_fs = q.make_fatsat(q.system)
    core.rf.freq_offset = 0.0
    out = {}
    for nombre, t1, t2, cs in [("sangre", q.T1_BLOOD, q.T2_BLOOD, 0.0),
                               ("musculo", q.T1_MUSCLE, q.T2_MUSCLE, 0.0),
                               ("grasa", q.T1_FAT, q.T2_FAT, -q.DF_FAT)]:
        df = cs + b0_offset
        out[nombre] = dict(
            t1=t1, t2=t2, df=df,
            f_inv=-bloch_rf(rf_inv, gz_inv.amplitude, [z], t1, t2, df=df)[2][0],
            fatsat=bloch_rf(rf_fs, 0.0, [0.0], t1, t2, df=df)[2][0],
            flip=np.deg2rad(flip_deg(bloch_rf(core.rf, core.gz.amplitude, [z], df=df)[2][0])),
        )
    return out


def static_tissue(core, eff, n_cycles=12, fatsat=True, qi_shift=0.0):
    """Tejido que no se mueve: ciclo completo hasta el estado estacionario.

    qi_shift desplaza la inversion respecto del tren (positivo = QI mas largo).
    Devuelve (|Mxy| por TR en el ultimo ciclo, Mz al inicio del tren).
    """
    tm, tr, n = core.timing, core.tr, q.N_KZ_ACQUIRED
    t_inv = tm.t_inv_center - qi_shift
    t_fs = tm.t_fatsat_start + q.RF_DEAD_TIME + q.FATSAT_DURATION / 2
    t_half = tm.t_first_alpha - tr / 2
    t_end = tm.t_first_alpha + n * tr
    t1, t2, df = eff["t1"], eff["t2"], eff["df"]
    m = np.array([0.0, 0.0, 1.0])
    for _ in range(n_cycles):
        m = np.array([0.0, 0.0, -eff["f_inv"] * m[2]])        # inversion + spoiler
        m = relax(m, t_fs - t_inv, t1, t2)
        if fatsat:
            m = np.array([0.0, 0.0, eff["fatsat"] * m[2]])    # fat-sat + spoiler
        m = relax(m, t_half - t_fs, t1, t2)
        mz_start = m[2]
        sig, m = bssfp_train(m, t1, t2, df, eff["flip"], tr, n)
        m = relax(m, q.TR_QISS + t_inv - t_end, t1, t2, df)
    return sig, mz_start


def flowing_blood(core, eff, mz_at_train_start):
    """Sangre que llega fresca en cada disparo: no hay estado estacionario entre ciclos."""
    m = np.array([0.0, 0.0, eff["fatsat"] * mz_at_train_start])
    return bssfp_train(m, eff["t1"], eff["t2"], eff["df"], eff["flip"], core.tr,
                       q.N_KZ_ACQUIRED)[0]


def simulate(core, b0_offset=0.0, z=0.0):
    """Senal en el eco de kz=0 de cada tejido.  Devuelve dict nombre -> |Mxy|."""
    e = pulse_effects(core, b0_offset, z)
    tm, k0 = core.timing, q.IDX_KZ0
    t_half = tm.t_first_alpha - core.tr / 2
    t_sat = t_half - tm.t_vsat_center            # saturacion venosa -> inicio del tren
    t_inv = t_half - tm.t_inv_center             # inversion -> inicio del tren
    r = {}
    r["arteria (sangre entrante)"] = flowing_blood(core, e["sangre"], 1.0)[k0]
    r["vena que fluye, saturada"] = flowing_blood(
        core, e["sangre"], 1 - np.exp(-t_sat / q.T1_BLOOD))[k0]
    r["vena con inversion unica (diseno anterior)"] = flowing_blood(
        core, e["sangre"], 1 - 2 * np.exp(-t_inv / q.T1_BLOOD))[k0]
    r["sangre estancada en el slab"] = static_tissue(core, e["sangre"])[0][k0]
    r["musculo"] = static_tissue(core, e["musculo"])[0][k0]
    r["grasa con fat-sat"] = static_tissue(core, e["grasa"])[0][k0]
    r["grasa sin fat-sat"] = static_tissue(core, e["grasa"], fatsat=False)[0][k0]
    return r, e


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _, core = q.build_sequence(n_slabs=1, n_views=1, verbose=False)
    tm = core.timing

    print("=" * 78)
    print(f"  CONTRASTE EN EL ECO DE kz=0  |  TR {core.tr*1e3:.2f} ms, FA {q.FLIP_ANGLE_DEG} deg, "
          f"{q.N_KZ_ACQUIRED} TRs, TR_QISS {q.TR_QISS*1e3:.0f} ms")
    print("=" * 78)
    res, eff = simulate(core)
    print("  Efecto medido de cada pulso (centro del slab, en resonancia):")
    for k, v in eff.items():
        print(f"    {k:8s}: eficiencia de inversion {v['f_inv']:.3f} | Mz tras fat-sat x{v['fatsat']:+.3f}"
              f" | flip {np.degrees(v['flip']):.1f} deg")
    art = res["arteria (sangre entrante)"]
    print("\n  tejido                                        |Mxy|     arteria : tejido")
    print("  " + "-" * 72)
    for k, v in res.items():
        ratio = "" if k.startswith("arteria") else f"{art/v:8.1f} : 1"
        print(f"  {k:44s} {v:7.4f}   {ratio}")

    print("\n" + "=" * 78)
    print("  A LO ANCHO DEL SLAB (el flip y la inversion varian con z)")
    print("=" * 78)
    print("   z (mm)   flip   ef.inv   arteria   musculo   grasa(fs)   art:musc")
    for z_mm in [0, 3, 6, 8, 9.75]:
        r, e = simulate(core, z=z_mm * 1e-3)
        print(f"   {z_mm:5.2f}   {np.degrees(e['musculo']['flip']):5.1f}   {e['musculo']['f_inv']:.3f}"
              f"    {r['arteria (sangre entrante)']:.4f}    {r['musculo']:.4f}     "
              f"{r['grasa con fat-sat']:.4f}    {r['arteria (sangre entrante)']/r['musculo']:6.1f}")

    print("\n" + "=" * 78)
    print("  ¿ESTA EL QI EN EL NULO?  musculo en kz=0 al desplazar la inversion")
    print("=" * 78)
    shifts = np.arange(-40e-3, 41e-3, 5e-3)
    e0 = pulse_effects(core)
    mus = np.array([static_tissue(core, e0["musculo"], qi_shift=s)[0][q.IDX_KZ0] for s in shifts])
    for s, m in zip(shifts, mus):
        marca = "  <- diseno" if abs(s) < 1e-9 else ""
        print(f"   QI {1e3*(q.QI+s):6.1f} ms : |Mxy| musculo = {m:.4f}{marca}")
    best = shifts[np.argmin(mus)]
    print(f"  minimo en QI = {1e3*(q.QI+best):.0f} ms (diseno {q.QI*1e3:.0f} ms); "
          f"en el diseno el musculo esta {mus[len(mus)//2]/mus.min():.1f}x sobre el minimo")

    print("\n" + "=" * 78)
    print("  SENSIBILIDAD A B0 (el fat-sat es estrecho: 79.6 Hz entre grasa y agua)")
    print("=" * 78)
    print("   B0 off   arteria   grasa(fs)   art:grasa   Mz grasa tras fat-sat")
    for off in [-30, -20, -10, 0, 10, 20, 30]:
        r, e = simulate(core, b0_offset=off)
        print(f"   {off:+4d} Hz   {r['arteria (sangre entrante)']:.4f}    {r['grasa con fat-sat']:.4f}"
              f"     {r['arteria (sangre entrante)']/r['grasa con fat-sat']:6.1f}        x{e['grasa']['fatsat']:+.3f}")

    # ---------------------------------------------------------------- grafico
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    e = pulse_effects(core)
    curvas = {
        "arteria (entrante)": flowing_blood(core, e["sangre"], 1.0),
        "vena saturada": flowing_blood(core, e["sangre"],
                                       1 - np.exp(-(tm.t_first_alpha - core.tr/2 - tm.t_vsat_center) / q.T1_BLOOD)),
        "musculo": static_tissue(core, e["musculo"])[0],
        "grasa con fat-sat": static_tissue(core, e["grasa"])[0],
        "grasa sin fat-sat": static_tissue(core, e["grasa"], fatsat=False)[0],
    }
    for k, v in curvas.items():
        ax[0].plot(np.arange(1, len(v) + 1), v, "o-", ms=3, lw=1.2, label=k)
    ax[0].axvline(q.IDX_KZ0 + 1, color="k", ls=":", lw=0.8, label="kz = 0")
    ax[0].set_xlabel("TR dentro del disparo"); ax[0].set_ylabel("$|M_{xy}|$")
    ax[0].set_title("Evolucion durante la ventana de lectura"); ax[0].legend(fontsize=7)
    ax[0].grid(alpha=0.3)
    ax[1].plot(1e3 * (q.QI + shifts), mus, "o-", ms=3)
    ax[1].axvline(q.QI * 1e3, color="k", ls=":", lw=0.8, label="QI de diseno")
    ax[1].set_xlabel("QI [ms]"); ax[1].set_ylabel("$|M_{xy}|$ musculo en kz=0")
    ax[1].set_title("Posicion del nulo del fondo"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    plt.tight_layout(); plt.savefig("contraste_bssfp.png", dpi=130)
    print("\nGrafico guardado: contraste_bssfp.png")
