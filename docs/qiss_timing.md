# 3D tsSOS-QISS bSSFP a 0.55 T — de dónde sale cada número

Registro de diseño de `qiss3d_tssos_bssfp_055T.py`. Cada valor viene de una fuente
citada, de una medición sobre la secuencia ensamblada (`verify_qiss3d.py`) o de una
simulación de Bloch de los pulsos reales (`sim_inversion_profile.py`,
`sim_bssfp_contrast.py`). Lo que no está verificado figura en la última sección.

Esta versión corrige los errores encontrados por una auditoría independiente de la
primera; la sección 8 dice cuáles eran.

## 1. Hardware y datos de partida

| Parámetro | Valor | Fuente |
|---|---|---|
| Campo / bore | 0.55 T / 80 cm | Varghese et al. 2023 |
| Gradientes | 26 mT/m, 45 T/m/s por eje (se usa el 85 %) | Varghese et al. 2023 |
| ECG integrado | no → secuencia sin gating | Varghese et al. 2023 |
| Tiempos muertos | RF 100 µs, ringdown 20 µs, ADC 10 µs | leídos del scanner |
| Intérprete Pulseq | 1.4.2 | scanner |
| T1/T2 sangre, grasa, músculo | 1122/263, 187/93, 450/55 ms | aportados por el usuario |
| Grasa–agua | 3.4 ppm × 23.42 MHz = 79.6 Hz | — |

Varghese et al., *Front Cardiovasc Med* 2023, doi:10.3389/fcvm.2023.1120982.
Geometría y esquema: Koktzoglou et al., *MRM* 2020, doi:10.1002/mrm.28339;
Edelman et al., *MRM* 2020, doi:10.1002/mrm.28032.

## 2. El disparo

```
t = 12 ms   inversión de fondo, WURST 6 kHz / 24 ms, 24.1 mm        (centro del pulso)
t = 29 ms   saturación venosa, sinc 90° 10 ms, banda de 60 mm por encima del slab
            spoiler
            intervalo quiescente
t = 123 ms  fat-sat espectral, gauss 15 ms a −79.6 Hz, spoiler
t = 143 ms  α/2, luego rephase + prefase de slab
t = 146 ms  primer α; 15 TR de 7.02 ms, un ángulo radial, todas las particiones
t = 220 ms  eco del TR #11, que es el que muestrea kz = 0
t = 250 ms  fin del tren; relleno hasta TR_QISS = 700 ms
```

## 3. Intervalo quiescente

El QI se mide **del centro de la inversión al eco del TR que muestrea kz = 0**. Con
partial Fourier 6/8 ese TR es el #11 de 15, no el central.

El fondo no vuelve al equilibrio entre inversiones y la inversión no es perfecta: el
T2 de 55 ms del músculo relaja durante los 24 ms del WURST y la eficiencia medida es
f = 0.88. En estado estacionario:

```
M⁻ = (1 − E) / (1 + f·E),   E = exp(−TR_QISS/T1)        Mz justo antes de invertir
TI = T1 · ln(1 + f·M⁻)                                   cruce por cero
```

Con T1 = 450 ms, TR_QISS = 700 ms y f = 0.88: **QI = 207.4 ms** (T1·ln2 daría 312).
La simulación del ciclo completo, que incluye el efecto del tren, pone el mínimo de
señal muscular en 207 ms.

**Por qué TR_QISS = 700 ms.** Con 600 ms el nulo cae en 190 ms y la sangre vertebral
(10 cm/s) solo recorre 19.0 mm, menos que el slab. Con 700 ms recorre 20.7 mm (1.06×);
la carótida, a 20 cm/s, 41.5 mm (2.13×). Este margen cuenta la sangre que entra hasta
el eco de kz = 0, incluida la que entra ya con el tren en marcha.

## 4. Tren bSSFP

| | Valor | Motivo |
|---|---|---|
| Excitación | sinc 2.0 ms, TBW 6, FWHM 23.4 mm | perfil; ver abajo |
| B1 de pico | 19.65 µT | escala como 1/duración |
| Lectura | 168 muestras, dwell 12 µs, 496 Hz/px | compromiso SNR / grasa |
| TR / TE | 7.02 / 3.51 ms | **derivado**, no impuesto |
| Particiones | 15 de 20, 1.3 mm, FOV_z 26 mm | aliasing en kz |

**Perfil de excitación.** `slice_thickness` de PyPulseq es el ancho a media altura, así
que para tener flip útil en los 19.5 mm conservados se excita 1.2× más grueso. El flip
medido es 100° hasta 3 mm, 94° a 6 mm, 80° a 8 mm y 62° en el borde. Cae bajo 3° a
15.7 mm del centro; con FOV_z = 26 mm solo se repliega sobre el slab lo que está más
allá de 16.25 mm. Con las 18 particiones del paper (FOV_z 23.4 mm) habría aliasing: el
paper usa FLASH y puede permitirse un RF largo y selectivo, bSSFP no.

**El TR se deriva.** Suma del bloque de RF, los rebobinados dimensionados para el kz
más exigente y la lectura. Los bloques tienen duración fija, así que todos los TR son
idénticos. El dwell va en múltiplos de 0.5 µs para que el ADC quede centrado con un
retardo entero en µs.

**Ancho de banda.** Un TR más corto acerca la grasa al nulo del stopband, en ±1/(2·TR):

| Hz/px | TR | arteria : grasa | SNR relativo |
|---|---|---|---|
| 356 | 7.68 ms | 9 : 1 | 1.00 |
| **500** | **7.0 ms** | **21 : 1** | **0.86** |
| 700 | 6.64 ms | 35 : 1 | 0.73 |

**Catalización.** Entre el α/2 y el primer α va un lóbulo de área −(área del selector):
el rephase del α/2 más la prefase del α.

## 5. Preparación

**Inversión de fondo.** El WURST solo invierte bien (Mz ≤ −0.8) en el 81 % central de
su espesor nominal, así que se pide 24.1 mm para cubrir los 19.5 mm del slab. A cambio
se pasa 2.35 mm por debajo: el 6 % del bolo carotídeo y el 11 % del vertebral llegan
atenuados. Es el precio de no tener el FOCI de los papers.

**Venas: saturación, no inversión.** La sangre yugular fluye y ve un solo pulso por
disparo, partiendo de M0. Tras una inversión única tiene Mz = −0.66 en el eco de kz = 0
y en magnitud sale casi tan brillante como la arteria (1.4 : 1). Tras una saturación
tiene Mz ≈ +0.16 (7 : 1). Los papers invierten, pero con QI de 500–700 ms. La banda de
60 mm cubre lo que la sangre venosa baja durante el QI hasta 30 cm/s; no toca el slab
(Mz mínima dentro: 0.988) y satura de verdad desde 7.7 mm por encima.

**Fat-sat.** Deja el agua en 0.996 y la grasa en 0.070. Aporta menos de lo que se creyó:
en el eco de kz = 0 la grasa ya casi olvidó la saturación (0.037 → 0.032). Baja la grasa
sobre todo en los primeros TR. Se conserva porque ocupa tiempo muerto del QI.

## 6. Contraste simulado en el eco de kz = 0

Centro del slab, en resonancia, ciclo completo en estado estacionario:

| Tejido | \|Mxy\| | arteria : tejido |
|---|---|---|
| Arteria (sangre entrante) | 0.655 | — |
| Vena que fluye, saturada | 0.090 | 7 : 1 |
| Sangre estancada en el slab | 0.067 | 10 : 1 |
| Grasa | 0.032 | 21 : 1 |
| Músculo | 0.001 | > 500 : 1 |

En el borde del slab (9.75 mm) el flip baja a 62°: la arteria da 0.476 y el músculo
0.012 (39 : 1).

**Sensibilidad a B0.** El contraste con la grasa baja a entre 6 : 1 y 18 : 1 con ±10–20 Hz y a 5 : 1
con ±30 Hz. Es el punto más frágil: el shim sobre el cuello importa.

**Límites del modelo.** La sangre que fluye se trata como si viera el tren entero desde
el inicio; no hay simulación espacial de flujo. Las venas lentas (< 10 cm/s) no se
reemplazan del todo durante el QI y quedan entre el caso "saturada" y el "estancada".

## 7. Geometría y examen

19 slabs de 19.5 mm con 4.5 mm de solape → 289.5 mm. Resolución 1.2 × 1.2 × 1.3 mm.
32 vistas radiales por slab (submuestreo 8×) en ángulo áureo, pensado para una
reconstrucción iterativa. **7:05** (19 × 32 × 700 ms). 41 952 bloques, 2.0 MB.

## 8. Qué se corrigió respecto de la primera versión

| Error | Consecuencia que tenía | Corrección |
|---|---|---|
| Venas por inversión, evaluadas como sangre estática | venas brillantes (1.4 : 1, no 104 : 1) | saturación |
| Sin prefase de slab antes del primer α | 1.9 ciclos de fase en el slab, señal oscilando entre TR, ghost en kz | lóbulo completo tras el α/2 |
| Excitación de 600 µs / TBW 3 | flip de 44° en el borde y aliasing en kz | 2 ms / TBW 6, slab 1.2×, 20 particiones |
| Inversión del mismo espesor que el slab | fondo sin suprimir en los bordes, bandas cada 15 mm | 24.1 mm |
| kz = 0 supuesto a mitad de ventana; QI desde el inicio del bloque; TR real 6.23 ≠ 6.28 ms | nulo corrido 20–30 ms | tiempos medidos y TR derivado |
| "El TR de 6.28 ms suprime la grasa" y luego "el fat-sat da 22 : 1 frente a 3.6 : 1" | ambas cifras salían de condiciones iniciales erróneas | simulación del ciclo completo |
| Verificación que sumaba áreas de eventos | no detectó la prefase faltante | se mide el momento en el centro de cada RF |
| `adc_samples_divisor` en `Opts` | no corría con el pypulseq 1.4.2 de PyPI | eliminado; probado con ambas versiones |

## 9. Pendiente antes del scanner

1. **B1 máximo.** La excitación pide 19.65 µT y PyPulseq no lo comprueba. En Siemens,
   la *reference voltage* es la de un rectangular de 180° en 1 ms, que son 11.74 µT:
   `B1_max = 11.74 µT × V_amplificador_max / V_referencia`. Si no alcanza, alargar
   `EXC_DURATION`. Los WURST (9.37 µT) degradan suavemente: 0.96 de eficiencia a 0.8×
   su B1, 0.91 a 0.7×.
2. **Sentido de z.** Que +z lógico sea craneal depende de la posición del paciente en
   el protocolo. Si está invertido, la banda venosa cae bajo el slab y satura el bolo
   arterial. Probar con un fantoma marcado.
3. **PNS.** Por eje se respeta el límite (22.0 mT/m máx.), pero el vector llega a
   28 mT/m y 53 T/m/s. Pasar el `.seq` por el chequeo del fabricante.
4. **Carga en el intérprete.** 41 952 bloques; probar antes con `--quick` (552).
5. **Shim.** Ver sensibilidad a B0 en la sección 6.
