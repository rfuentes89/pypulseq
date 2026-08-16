"""Contraste real del tren bSSFP de QISS a 0.55T, por simulacion de Bloch.

Verifica la afirmacion central del diseno, que hasta ahora estaba solo argumentada:

    con TR = 6.28 ms la grasa cae en el nulo del stopband de bSSFP y se suprime sola.

Hay un motivo concreto para dudar: el nulo del stopband es un fenomeno de ESTADO
ESTACIONARIO, pero la lectura de QISS son solo 14 TRs (88 ms) tras el intervalo
quiescente.  Si la grasa no alcanza su estado estacionario dentro de esa ventana, el
argumento no se sostiene y hay grasa brillante en las MIP.

Se simula el tren completo: estado tras la inversion -> preparacion alpha/2 -> 14 TRs
con alternancia de fase, muestreando en TE = TR/2.  Se compara con el estado
estacionario asintotico para ver cuanto se parece uno al otro.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import qiss3d_tssos_bssfp_055T as q

TR = q.TR_BSSFP
TE = TR / 2
ALPHA = np.deg2rad(q.FLIP_ANGLE_DEG)
N_TR = q.N_KZ_ACQUIRED
IDX_CENTRO_K = int(np.where(q.KZ_INDICES == 0)[0][0])   # TR que muestrea kz = 0


def rot_rf(alpha, phase):
    """Rotacion por un pulso duro de angulo alpha en torno a un eje transversal
    orientado segun phase (aproximacion de pulso instantaneo)."""
    c, s = np.cos(alpha), np.sin(alpha)
    cp, sp = np.cos(phase), np.sin(phase)
    # Rotacion de alpha en torno al eje (cos(phase), sin(phase), 0)
    return np.array([
        [cp * cp + sp * sp * c, cp * sp * (1 - c), sp * s],
        [cp * sp * (1 - c), sp * sp + cp * cp * c, -cp * s],
        [-sp * s, cp * s, c],
    ])


def free_precess(dt, t1, t2, df):
    """Precesion libre con relajacion durante dt, con desfase df en Hz."""
    phi = 2 * np.pi * df * dt
    e1, e2 = np.exp(-dt / t1), np.exp(-dt / t2)
    c, s = np.cos(phi), np.sin(phi)
    a = np.array([[e2 * c, e2 * s, 0],
                  [-e2 * s, e2 * c, 0],
                  [0, 0, e1]])
    b = np.array([0, 0, 1 - e1])
    return a, b


def bssfp_train(t1, t2, df, mz0, n_tr=N_TR, alpha=ALPHA, tr=TR):
    """Corre alpha/2 + n_tr TRs de bSSFP y devuelve |Mxy| en TE de cada TR."""
    m = np.array([0.0, 0.0, mz0])

    # Preparacion alpha/2 con fase opuesta al primer alpha, seguida de TR/2:
    # lleva la magnetizacion al centro de la elipse del estado estacionario y
    # evita la oscilacion transitoria de las primeras vistas.
    m = rot_rf(alpha / 2, np.pi) @ m
    a, b = free_precess(tr / 2, t1, t2, df)
    m = a @ m + b

    señales = []
    for n in range(n_tr):
        m = rot_rf(alpha, np.pi * (n % 2)) @ m
        a, b = free_precess(TE, t1, t2, df)
        m = a @ m + b
        señales.append(np.hypot(m[0], m[1]))
        a, b = free_precess(tr - TE, t1, t2, df)
        m = a @ m + b
    return np.array(señales)


def mz_tras_inversion(t1, tr_qiss=q.TR_QISS, ti=None):
    """Mz fraccional en el instante TI, en estado estacionario de inversiones."""
    ti = q.QI if ti is None else ti
    e = np.exp(-tr_qiss / t1)
    x = (1 - e) / (1 + e)
    return 1 - (1 + x) * np.exp(-ti / t1)


# =============================================================================
print("=" * 78)
print("  CONTRASTE DEL TREN bSSFP  |  TR = %.2f ms, FA = %d deg, %d TRs"
      % (TR * 1e3, q.FLIP_ANGLE_DEG, N_TR))
print("=" * 78)
print(f"  Nulo del stopband en +-1/(2*TR) = +-{1/(2*TR):.1f} Hz")
print(f"  Grasa a {-q.DF_FAT:.1f} Hz respecto del agua")
print(f"  kz = 0 se muestrea en el TR #{IDX_CENTRO_K + 1} de {N_TR}\n")

tejidos = [
    # nombre,                    T1,          T2,          df [Hz],   Mz inicial
    ("sangre arterial entrante", q.T1_BLOOD, q.T2_BLOOD, 0.0, 1.0),
    ("sangre estatica/venosa",   q.T1_BLOOD, q.T2_BLOOD, 0.0, mz_tras_inversion(q.T1_BLOOD)),
    ("musculo (fondo)",          q.T1_MUSCLE, q.T2_MUSCLE, 0.0, mz_tras_inversion(q.T1_MUSCLE)),
    # Mz de la grasa tras el fat-sat espectral de 15 ms.  El valor sale de la
    # simulacion del perfil espectral del pulso: +0.070 en resonancia, y hasta
    # +0.272 con una desviacion de B0 de 20 Hz (caso peor que tambien se reporta).
    ("grasa (con fat-sat)",      q.T1_FAT,   q.T2_FAT,   -q.DF_FAT, 0.070),
    ("grasa (sin fat-sat)",      q.T1_FAT,   q.T2_FAT,   -q.DF_FAT, mz_tras_inversion(q.T1_FAT)),
    ("grasa (fat-sat, B0 20Hz)", q.T1_FAT,   q.T2_FAT,   -q.DF_FAT, 0.272),
]

print("  tejido                     Mz inicial   |Mxy| en kz=0   |Mxy| medio")
print("  " + "-" * 68)
resultados = {}
for nombre, t1, t2, df, mz0 in tejidos:
    s = bssfp_train(t1, t2, df, mz0)
    resultados[nombre] = s
    print(f"  {nombre:26s}  {mz0:+7.3f}      {s[IDX_CENTRO_K]:9.4f}      {s.mean():9.4f}")

s_art = resultados["sangre arterial entrante"][IDX_CENTRO_K]
print("\n  Contraste respecto de la sangre arterial entrante (en kz = 0):")
for nombre in resultados:
    if nombre == "sangre arterial entrante":
        continue
    ratio = resultados[nombre][IDX_CENTRO_K] / s_art
    print(f"    arteria / {nombre:26s} = {1/ratio if ratio>0 else np.inf:6.1f} : 1")

# =============================================================================
print("\n" + "=" * 78)
print("  ¿LLEGA LA GRASA A SU ESTADO ESTACIONARIO EN 14 TRs?")
print("=" * 78)
s_fat_largo = bssfp_train(q.T1_FAT, q.T2_FAT, -q.DF_FAT,
                          mz_tras_inversion(q.T1_FAT), n_tr=400)
ss_fat = s_fat_largo[-50:].mean()
print(f"  |Mxy| de grasa en estado estacionario (400 TRs) : {ss_fat:.5f}")
print(f"  |Mxy| de grasa en el TR de kz=0 (#{IDX_CENTRO_K+1})            : "
      f"{resultados['grasa (con fat-sat)'][IDX_CENTRO_K]:.5f}")
print(f"  |Mxy| de grasa promediada en los {N_TR} TRs        : {resultados['grasa (con fat-sat)'].mean():.5f}")
factor = resultados['grasa (con fat-sat)'][IDX_CENTRO_K] / ss_fat if ss_fat > 0 else np.inf
print(f"  -> en kz=0 la grasa esta {factor:.1f}x por encima de su estado estacionario")
if resultados['grasa (con fat-sat)'][IDX_CENTRO_K] / s_art < 0.1:
    print("  -> aun asi queda por debajo del 10 % de la senal arterial: SUPRESION EFECTIVA")
else:
    print("  -> queda por ENCIMA del 10 % de la senal arterial: la supresion NO alcanza")

# =============================================================================
print("\n" + "=" * 78)
print("  SENSIBILIDAD A DESPLAZAMIENTOS DE B0")
print("=" * 78)
print("""  El nulo es estrecho: si el agua se desplaza de resonancia, la grasa sale del
  nulo y ademas el agua se acerca a el. Se barre el offset de B0.
""")
print("  offset B0    |Mxy| arteria   |Mxy| grasa   cociente")
print("  " + "-" * 54)
for off in [-40, -20, -10, 0, 10, 20, 40]:
    sa = bssfp_train(q.T1_BLOOD, q.T2_BLOOD, off, 1.0)[IDX_CENTRO_K]
    sf = bssfp_train(q.T1_FAT, q.T2_FAT, off - q.DF_FAT,
                     mz_tras_inversion(q.T1_FAT))[IDX_CENTRO_K]
    print(f"  {off:+5d} Hz     {sa:9.4f}     {sf:9.4f}    {sa/sf if sf>0 else np.inf:7.1f} : 1")

# =============================================================================
# Grafico: respuesta frente a off-resonancia
df_sweep = np.linspace(-160, 160, 321)
resp = {}
for nombre, t1, t2, _df, mz0 in tejidos:
    if nombre in ("sangre estatica/venosa", "grasa (fat-sat, B0 20Hz)"):
        continue
    resp[nombre] = np.array([bssfp_train(t1, t2, d, mz0)[IDX_CENTRO_K] for d in df_sweep])

fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
for nombre, curva in resp.items():
    axes[0].plot(df_sweep, curva, lw=1.5, label=nombre)
axes[0].axvline(-q.DF_FAT, color="tab:red", ls="--", lw=1,
                label=f"grasa a {-q.DF_FAT:.0f} Hz")
axes[0].axvline(1 / (2 * TR), color="k", ls=":", lw=0.8)
axes[0].axvline(-1 / (2 * TR), color="k", ls=":", lw=0.8,
                label=f"nulos a $\\pm${1/(2*TR):.0f} Hz")
axes[0].set_xlabel("off-resonancia [Hz]"); axes[0].set_ylabel("$|M_{xy}|$ en kz=0")
axes[0].set_title(f"Respuesta bSSFP, TR = {TR*1e3:.2f} ms")
axes[0].legend(fontsize=7); axes[0].grid(alpha=0.3)

for nombre, s in resultados.items():
    axes[1].plot(np.arange(1, N_TR + 1), s, "o-", ms=3, lw=1.2, label=nombre)
axes[1].axvline(IDX_CENTRO_K + 1, color="k", ls=":", lw=0.8, label="kz = 0")
axes[1].set_xlabel("TR dentro del disparo"); axes[1].set_ylabel("$|M_{xy}|$")
axes[1].set_title("Evolucion durante la ventana de lectura")
axes[1].legend(fontsize=7); axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("contraste_bssfp.png", dpi=130)
print("\nGrafico guardado: contraste_bssfp.png")
