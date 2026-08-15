"""Corrección del borde blando de la banda venosa.

sim_inversion_profile.py encontró que la banda venosa de 100 mm invade el slab de
imagen por 6 mm: su transicion inferior mide 22 mm y solo hay 5 mm de hueco.

La causa es geometrica, no del pulso: el ancho de transicion en Hz lo fija el pulso,
pero al pasarlo a milimetros se divide por la fuerza del gradiente selector.  Como la
banda venosa es 5 veces mas gruesa que el slab, su gradiente es 5 veces mas debil
(0.293 frente a 1.50 mT/m) y la MISMA transicion en Hz se estira 5 veces en mm.

Dos palancas: alargar el pulso (afila la transicion en Hz) o abrir el hueco entre
banda y slab (acepta la transicion y la aparta).  Se cuantifican ambas.
"""

import numpy as np
import pypulseq as pp

import qiss3d_tssos_bssfp_055T as q
from sim_inversion_profile import bloch_profile, transition_width

GAP_ACTUAL = q.INV_VENOUS_GAP * 1e3            # mm
T1B, T2B = q.T1_BLOOD, q.T2_BLOOD

print("=" * 78)
print("  BANDA VENOSA: ¿alargar el pulso o abrir el hueco?")
print("=" * 78)
print(f"\nSlab de imagen {q.SLAB_THICKNESS*1e3:.1f} mm | banda venosa "
      f"{q.INV_VENOUS_THICKNESS*1e3:.0f} mm | hueco actual {GAP_ACTUAL:.1f} mm\n")

print("  dur    tipo     B1pico   transicion   fuga bajo    hueco     veredicto")
print("  (ms)            (uT)        (mm)      borde (mm)  necesario")
print("  " + "-" * 72)

z = np.linspace(-80e-3, 80e-3, 321)
resultados = []

for dur_ms, ptype in [(8, "hypsec"), (12, "hypsec"), (16, "hypsec"),
                      (24, "hypsec"), (32, "hypsec"),
                      (8, "wurst"), (16, "wurst"), (24, "wurst")]:
    dur = dur_ms * 1e-3
    kwargs = dict(pulse_type=ptype, duration=dur,
                  slice_thickness=q.INV_VENOUS_THICKNESS,
                  system=q.system, return_gz=True, use="inversion")
    if ptype == "wurst":
        # El WURST toma el ancho de banda explicito; se pide el que cubre la banda.
        kwargs["bandwidth"] = 4000
    rf, gz, _ = pp.make_adiabatic_pulse(**kwargs)
    rf.freq_offset = 0.0

    mz = bloch_profile(rf, gz.amplitude, z, T1B, T2B)
    w, _, z_free = transition_width(z, mz, "inferior")
    if np.isnan(w):
        print(f"  {dur_ms:4.0f}   {ptype:7s}  "
              f"{np.abs(rf.signal).max()/q.GAMMA_BAR*1e6:6.2f}   "
              f"    ---        ---         ---     no invierte bien")
        continue

    fuga = abs(z_free) - q.INV_VENOUS_THICKNESS / 2 * 1e3
    hueco_nec = fuga
    ok = hueco_nec <= GAP_ACTUAL
    veredicto = "OK con hueco actual" if ok else f"exige hueco {hueco_nec:.0f} mm"
    b1 = np.abs(rf.signal).max() / q.GAMMA_BAR * 1e6
    print(f"  {dur_ms:4.0f}   {ptype:7s}  {b1:6.2f}   {w:8.2f}     {fuga:7.2f}    "
          f"{hueco_nec:7.1f}     {veredicto}")
    resultados.append((dur_ms, ptype, w, fuga, b1))

# ---------------------------------------------------------------- coste de abrir el hueco
print("\n" + "=" * 78)
print("  COSTE DE ABRIR EL HUECO: sangre venosa que se cuela")
print("=" * 78)
print("""
Abrir el hueco aparta la transicion del slab, pero deja una franja de sangre venosa
sin invertir justo encima del slab.  Esa sangre baja hacia el slab durante el QI.
La pregunta es si alcanza a entrar antes de la lectura.
""")
print(f"  QI = {q.QI*1e3:.0f} ms")
print("\n  hueco   sangre venosa sin suprimir     ¿alcanza a entrar al slab?")
print("  (mm)    a velocidad 5 / 10 / 20 cm/s")
print("  " + "-" * 62)
for hueco in [5, 8, 11, 14, 20]:
    entradas = []
    for v in [0.05, 0.10, 0.20]:
        recorrido = v * q.QI * 1e3          # mm que baja la sangre en el QI
        entra = recorrido > hueco
        entradas.append("si" if entra else "no")
    print(f"  {hueco:4d}    {' / '.join(f'{v*q.QI*1e3:5.1f}mm' for v in [0.05,0.10,0.20])}"
          f"      {' / '.join(entradas)}")

print("""
  A cualquier velocidad venosa realista la sangre recorre 10-40 mm en el QI, o sea
  que SIEMPRE cruza el hueco y entra al slab. Abrir el hueco no evita que entre
  sangre venosa: solo cambia CUANTA llega sin invertir. Con hueco H y velocidad v,
  la sangre que estaba entre 0 y H mm sobre el slab entra sin haber sido invertida.
""")

# ---------------------------------------------------------------- recomendacion
print("=" * 78)
print("  LECTURA DE LOS NUMEROS")
print("=" * 78)
if resultados:
    viables = [r for r in resultados if r[3] <= GAP_ACTUAL]
    if viables:
        mejor = min(viables, key=lambda r: r[0])
        print(f"\n  El pulso mas corto que cabe en el hueco actual de {GAP_ACTUAL:.0f} mm es")
        print(f"  {ptype} de {mejor[0]:.0f} ms: transicion {mejor[2]:.1f} mm, "
              f"fuga {mejor[3]:.1f} mm, B1 pico {mejor[4]:.1f} uT.")
    else:
        mejor = min(resultados, key=lambda r: r[3])
        print(f"\n  Ninguna variante cabe en el hueco actual de {GAP_ACTUAL:.0f} mm.")
        print(f"  La mejor es {mejor[1]} de {mejor[0]:.0f} ms con fuga de {mejor[3]:.1f} mm,")
        print(f"  que exigiria abrir el hueco a >= {mejor[3]:.0f} mm.")
print(f"""
  Presupuesto temporal: el modulo de preparacion ocupa hoy ~20 ms de los
  {q.QI*1e3:.0f} ms del QI, asi que alargar los adiabaticos es casi gratis en tiempo
  de examen. El limite real es el B1 de pico: hay que confirmar el maximo del cuerpo
  del Free.Max antes de dar por bueno cualquiera de estos pulsos.
""")
