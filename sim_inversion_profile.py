"""Perfiles espaciales y espectrales de los pulsos reales, por integracion de Bloch.

Mide lo que decide si cada pulso sirve:
  * excitacion bSSFP: flip dentro del slab conservado y cuanto se excita fuera del
    FOV_z (lo que se repliega sobre el slab);
  * inversion de fondo: que la zona bien invertida CUBRA el slab (la primera version
    solo media la fuga hacia afuera y no vio que la transicion caia dentro), y cuanto
    se pasa por debajo, de donde viene la sangre arterial;
  * saturacion venosa: que no toque el slab;
  * fat-sat: selectividad espectral y tolerancia a B0;
  * WURST frente a B1: degradacion al bajar la amplitud.
"""

import copy

import numpy as np

import qiss3d_tssos_bssfp_055T as q
from bloch_tools import bloch_rf, flip_deg

HALF = q.SLAB_THICKNESS / 2


def first_beyond(z, cond):
    """Menor |z| a partir del cual cond se cumple para todo z mayor."""
    bad = np.where(~cond)[0]
    return z[min(bad[-1] + 1, len(z) - 1)] if len(bad) else z[0]


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    core = q.make_bssfp_core(q.system)
    z = np.linspace(0, 40e-3, 801)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    print("=" * 76 + "\n  1. EXCITACION bSSFP\n" + "=" * 76)
    fl = flip_deg(bloch_rf(core.rf, core.gz.amplitude, z)[2])
    for zz in [0, 3, 6, 8, HALF * 1e3]:
        print(f"   flip en z = {zz:5.2f} mm : {fl[np.argmin(abs(z - zz*1e-3))]:5.1f} deg")
    z3 = first_beyond(z, fl < 3.0)
    z_alias = q.FOV_Z - HALF          # a partir de aqui se repliega sobre el slab conservado
    print(f"   flip < 3 deg a partir de |z| = {z3*1e3:.1f} mm")
    print(f"   se repliega sobre el slab lo que este en |z| > FOV_z - {HALF*1e3:.2f} = {z_alias*1e3:.2f} mm")
    print(f"   -> {'OK: lo que se repliega esta excitado a < 3 deg' if z3 <= z_alias else 'ALIASING: subir N_KZ_GRID o EXC_TBW'}")
    ax[0].plot(z * 1e3, fl); ax[0].axvline(HALF * 1e3, color="k", ls=":", label="borde del slab")
    ax[0].axvline(z_alias * 1e3, color="r", ls="--", label="limite de aliasing")
    ax[0].set_xlabel("|z| [mm]"); ax[0].set_ylabel("flip [deg]"); ax[0].set_title("Excitacion")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    print("\n" + "=" * 76 + "\n  2. INVERSION DE FONDO (WURST)\n" + "=" * 76)
    rf_inv, gz_inv = q.make_adiabatic_inversion(q.system, q.INV_SLAB_THICKNESS, 0.0)
    for nombre, t1, t2 in [("sangre", q.T1_BLOOD, q.T2_BLOOD), ("musculo", q.T1_MUSCLE, q.T2_MUSCLE)]:
        mz = bloch_rf(rf_inv, gz_inv.amplitude, z, t1, t2)[2]
        print(f"   {nombre:8s}: Mz centro {mz[0]:+.3f} | Mz en el borde del slab "
              f"{mz[np.argmin(abs(z - HALF))]:+.3f}")
        if nombre == "sangre":
            z_in = z[np.where(mz <= -0.8)[0][-1]]
            z_out = first_beyond(z, mz >= 0.9)
            ax[1].plot(z * 1e3, mz, label="inversion de fondo")
    print(f"   zona bien invertida (Mz <= -0.8) hasta |z| = {z_in*1e3:.2f} mm "
          f"-> {'CUBRE' if z_in >= HALF - 1e-4 else 'NO CUBRE'} el slab ({HALF*1e3:.2f} mm)")
    fuga = z_out - HALF
    print(f"   sangre intacta (Mz >= 0.9) desde |z| = {z_out*1e3:.2f} mm: fuga de {fuga*1e3:.2f} mm bajo el slab")
    for vaso, v in [("carotida", 0.20), ("vertebral", 0.10)]:
        print(f"      {vaso:9s}: {fuga/(v*q.QI)*100:4.1f} % del bolo que entra en el QI llega atenuado")

    print("\n" + "=" * 76 + "\n  3. SATURACION VENOSA\n" + "=" * 76)
    rf_vs, gz_vs = q.make_venous_saturation(q.system, 0.0)
    zc = HALF + q.VSAT_GAP + q.VSAT_THICKNESS / 2          # centro de la banda desde el del slab
    zz = np.linspace(-10e-3, 90e-3, 1001)                  # z medido desde el centro del slab
    mz = bloch_rf(rf_vs, gz_vs.amplitude, zz - zc, q.T1_BLOOD, q.T2_BLOOD)[2]
    en_slab = mz[zz <= HALF].min()
    z_sat = zz[np.where(mz <= 0.2)[0][0]]
    print(f"   Mz minima dentro del slab: {en_slab:+.3f} -> {'no lo toca' if en_slab > 0.95 else 'INVADE el slab'}")
    print(f"   saturacion efectiva (Mz <= 0.2) desde {(z_sat-HALF)*1e3:.1f} mm por encima del slab")
    ax[1].plot(zz * 1e3, mz, label="saturacion venosa (por encima)")
    ax[1].axvspan(-HALF * 1e3, HALF * 1e3, alpha=0.12, label="slab conservado")
    ax[1].set_xlim(-12, 45); ax[1].set_xlabel("z desde el centro del slab [mm]")
    ax[1].set_ylabel("$M_z$"); ax[1].set_title("Preparacion"); ax[1].legend(fontsize=8)
    ax[1].grid(alpha=0.3)

    print("\n" + "=" * 76 + "\n  4. FAT-SAT ESPECTRAL\n" + "=" * 76)
    rf_fs = q.make_fatsat(q.system)
    print("   B0 off    Mz agua    Mz grasa")
    for off in [-20, -10, 0, 10, 20]:
        w = bloch_rf(rf_fs, 0, [0], q.T1_BLOOD, q.T2_BLOOD, df=off)[2][0]
        f = bloch_rf(rf_fs, 0, [0], q.T1_FAT, q.T2_FAT, df=off - q.DF_FAT)[2][0]
        print(f"   {off:+4d} Hz    {w:+.3f}     {f:+.3f}")

    print("\n" + "=" * 76 + "\n  5. WURST FRENTE A B1\n" + "=" * 76)
    b1 = np.abs(rf_inv.signal).max() / q.GAMMA_BAR * 1e6
    zi = np.linspace(-6e-3, 6e-3, 25)
    esc, efs = np.arange(0.5, 1.21, 0.1), []
    for k in esc:
        r = copy.deepcopy(rf_inv); r.signal = rf_inv.signal * k
        efs.append(-bloch_rf(r, gz_inv.amplitude, zi)[2].mean())
        print(f"   {k:4.2f} x ({k*b1:5.2f} uT): eficiencia {efs[-1]:.3f}")
    ax[2].plot(esc * b1, efs, "o-"); ax[2].axhline(0.95, color="k", ls=":")
    ax[2].set_xlabel("B1 de pico [uT]"); ax[2].set_ylabel("eficiencia de inversion")
    ax[2].set_title("WURST frente a B1 (sin relajacion)"); ax[2].grid(alpha=0.3)

    plt.tight_layout(); plt.savefig("perfil_inversion.png", dpi=130)
    print("\nGrafico guardado: perfil_inversion.png")
