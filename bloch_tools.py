"""Herramientas de simulacion de Bloch compartidas por los scripts sim_*.py.

Convencion de signos (la de Pulseq): un RF modulado como exp(+i*2*pi*f*t) resuena en
los espines cuya frecuencia en el sistema rotante es +f.  Con un gradiente G [Hz/m] un
espin en z precesa a G*z, de modo que rf.freq_offset = G*z0 selecciona el plano z0.
La version anterior de estos integradores tenia el eje z espejado; era inocuo con
perfiles simetricos pero intercambiaba las etiquetas inferior/superior.
"""

import numpy as np


def resample_rf(rf, dt=1e-6):
    """Devuelve (senal compleja en Hz, tiempos) remuestreados a paso uniforme dt.

    Incorpora rf.freq_offset y rf.phase_offset, que PyPulseq guarda como campos
    aparte y NO modula dentro de rf.signal.
    """
    t = np.asarray(rf.t, float)
    s = np.asarray(rf.signal)
    tu = np.arange(t[0], t[-1] + dt / 2, dt)
    su = np.interp(tu, t, s.real) + 1j * np.interp(tu, t, s.imag)
    f0 = float(getattr(rf, "freq_offset", 0.0) or 0.0)
    p0 = float(getattr(rf, "phase_offset", 0.0) or 0.0)
    return su * np.exp(1j * (2 * np.pi * f0 * tu + p0)), tu


def bloch_rf(rf, g_hz_m, z, t1=np.inf, t2=np.inf, df=0.0, m0=(0.0, 0.0, 1.0)):
    """Integra Bloch durante un pulso RF real, vectorizado en z.

    rf       : evento RF de PyPulseq (senal en Hz = gamma*B1/2pi)
    g_hz_m   : amplitud del gradiente selector durante el pulso [Hz/m] (0 = no selectivo)
    z        : posiciones [m]
    df       : off-resonancia adicional [Hz] (desplazamiento quimico, B0)
    Devuelve M con forma (3, nz).  Rotacion exacta en torno al campo efectivo en cada
    paso de 1 us (formula de Rodrigues); sin aproximacion de angulo pequeno.
    """
    s, t = resample_rf(rf)
    dt = t[1] - t[0]
    z = np.atleast_1d(np.asarray(z, float))
    m = np.tile(np.asarray(m0, float)[:, None], (1, z.size))
    # En el sistema rotante a la frecuencia portadora el espin precesa a +(G z + df);
    # la rotacion de M en torno a B_eff es levogira, de ahi el signo del eje z.
    wz = -2 * np.pi * (g_hz_m * z + df)
    e1 = np.exp(-dt / t1) if np.isfinite(t1) else 1.0
    e2 = np.exp(-dt / t2) if np.isfinite(t2) else 1.0
    for sk in s:
        wx, wy = 2 * np.pi * sk.real, -2 * np.pi * sk.imag
        wn = np.sqrt(wx * wx + wy * wy + wz * wz)
        wn = np.where(wn == 0, 1e-30, wn)
        ax, ay, az = wx / wn, wy / wn, wz / wn
        th = wn * dt
        c, sn = np.cos(th), np.sin(th)
        dot = ax * m[0] + ay * m[1] + az * m[2]
        cx = ay * m[2] - az * m[1]
        cy = az * m[0] - ax * m[2]
        cz = ax * m[1] - ay * m[0]
        m = np.vstack([m[0] * c + cx * sn + ax * dot * (1 - c),
                       m[1] * c + cy * sn + ay * dot * (1 - c),
                       m[2] * c + cz * sn + az * dot * (1 - c)])
        m[0] *= e2
        m[1] *= e2
        m[2] = 1 + (m[2] - 1) * e1
    return m


def flip_deg(mz):
    return np.degrees(np.arccos(np.clip(mz, -1, 1)))


# ---------------------------------------------------------------- pulso duro + precesion
def rot_rf(alpha, phase):
    """Rotacion de angulo alpha en torno al eje transversal (cos phase, sin phase, 0)."""
    c, s = np.cos(alpha), np.sin(alpha)
    cp, sp = np.cos(phase), np.sin(phase)
    return np.array([[cp * cp + sp * sp * c, cp * sp * (1 - c), sp * s],
                     [cp * sp * (1 - c), sp * sp + cp * cp * c, -cp * s],
                     [-sp * s, cp * s, c]])


def relax(m, dt, t1, t2, df=0.0):
    """Precesion libre con relajacion durante dt."""
    phi = 2 * np.pi * df * dt
    e1, e2 = np.exp(-dt / t1), np.exp(-dt / t2)
    c, s = np.cos(phi), np.sin(phi)
    return np.array([e2 * (c * m[0] + s * m[1]),
                     e2 * (-s * m[0] + c * m[1]),
                     1 + (m[2] - 1) * e1])


def bssfp_train(m_start, t1, t2, df, alpha, tr, n_tr, prep_half=True):
    """alpha/2 - TR/2 - (alpha - TR)xN con alternancia de fase.  Parte del vector
    m_start y devuelve (|Mxy| en TE=TR/2 de cada TR, M al final del ultimo TR)."""
    m = np.asarray(m_start, float).copy()
    if prep_half:
        m = rot_rf(alpha / 2, np.pi) @ m
        m = relax(m, tr / 2, t1, t2, df)
    sig = np.empty(n_tr)
    for n in range(n_tr):
        m = rot_rf(alpha, np.pi * (n % 2)) @ m
        m = relax(m, tr / 2, t1, t2, df)
        sig[n] = np.hypot(m[0], m[1])
        m = relax(m, tr / 2, t1, t2, df)
    return sig, m
