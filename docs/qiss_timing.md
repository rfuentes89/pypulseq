# Derivación del timing — 3D tsSOS-QISS bSSFP a 0.55T

Registro de dónde sale cada número de `qiss3d_tssos_bssfp_055T.py`. Ninguna constante
del script es un valor heredado de protocolos de campo alto: o se deriva acá, o viene
de una fuente citada, o está marcada como pendiente de confirmar.

## 1. Hardware de partida

| Parámetro | Valor | Fuente |
|---|---|---|
| B0 | 0.55 T | MAGNETOM Free.Max |
| Bore | 80 cm | Varghese et al. 2023 |
| Gradiente máximo | 26 mT/m (por eje) | Varghese et al. 2023 |
| Slew máximo | 45 T/m/s (por eje) | Varghese et al. 2023 |
| ECG integrado | **no lo trae** | Varghese et al. 2023 (usaron trigger externo) |

Varghese et al., *Front Cardiovasc Med* 2023, doi:10.3389/fcvm.2023.1120982.

La ausencia de ECG de fábrica es la razón de que esta implementación sea **ungated** y
siga el esquema de Koktzoglou (3T, cuello) en vez del de Edelman (1.5T, ECG-gated).

## 2. Tiempos de relajación a 0.55T

| Tejido | T1 (ms) | T2 (ms) |
|---|---|---|
| Sangre | 1122 | 263 |
| Grasa | 187 | 93 |
| Músculo | 450 | 55 |

Aportados por el usuario. **No** son extrapolaciones de 1.5T/3T, que a este campo
darían errores grandes: el T1 se acorta y el T2 se alarga al bajar B0.

Razón T2/T1 de la sangre = 0.234, que es lo que fija el contraste alcanzable en bSSFP.

## 3. TR de bSSFP: la grasa cae en el stopband

Larmor a 0.55 T: 42.576 MHz/T × 0.55 T = **23.42 MHz**.

Separación grasa-agua con el desplazamiento convencional de 3.4 ppm:

```
Δf = 3.4e-6 × 23.42e6 = 79.6 Hz
```

En bSSFP con alternancia de fase RF de 0/180°, la banda de paso se centra en resonancia
y los nulos caen en ±1/(2·TR). Poner la grasa exactamente en un nulo exige:

```
TR = 1 / (2 × 79.6 Hz) = 6.28 ms
```

**Por qué esto no funciona a 3T y sí acá:** a 3 T la separación es 434 Hz y el TR
requerido sería 1.15 ms, irrealizable. Por eso Koktzoglou et al. tuvieron que suprimir
grasa con un readout multi-eco muestreando en TEs de oposición de fase (1.6/3.7/5.7 ms).
Copiar esos TEs a 0.55 T daría cero supresión: acá la oposición de fase cae en 6.25 ms,
así que 1.6 ms está prácticamente *en fase*.

**Y encaja con el hardware.** El TR mínimo alcanzable con 26 mT/m y 45 T/m/s a 1.2 mm
in-plane es 5.96 ms (ver `feasibility_bssfp_055T.py`), siempre que el rephase de slab,
la codificación kz y la prefase de lectura se toquen **concurrentes** en un solo bloque:
encadenados en serie el TR se va a 8.12 ms y el truco deja de ser posible. Los 0.32 ms
de holgura se invierten en alargar el readout hasta un ancho de banda de **356 Hz/px**,
que a bajo campo es exactamente donde conviene gastarlos (menos BW = más SNR).

## 4. El intervalo quiescente cumple dos funciones a la vez

El QI es simultáneamente el TI que nulea el fondo y el tiempo de entrada de sangre
fresca. Son requisitos opuestos: acortarlo mejora la supresión de fondo y empeora el
inflow.

### 4.1 Nulling con recuperación incompleta

Con TR_QISS de 600 ms y T1 muscular de 450 ms, el fondo **no vuelve al equilibrio**
entre inversiones sucesivas. Usar el T1·ln2 de libro sería un error de 105 ms.

En estado estacionario, la Mz justo antes de cada inversión vale

```
x = (1 - E) / (1 + E),      E = exp(-TR/T1)
```

y tras invertir (−x) la recuperación cruza cero en

```
TI = T1 · ln(1 + x)
```

Para T1 = 450 ms y TR = 600 ms: **TI = 206.6 ms**, contra los 311.9 ms de T1·ln2.
Con TR ≫ T1 la expresión tiende al clásico T1·ln2, como debe ser.

### 4.2 Inflow

Velocidades medias promediadas en el ciclo cardiaco citadas por Koktzoglou et al.:
carótida 20 cm/s, vertebral 10 cm/s. Con TI = 206.6 ms:

| Vaso | Inflow | vs. slab de 19.5 mm |
|---|---|---|
| Carótida | 41.3 mm | 2.12× |
| Vertebral | 20.7 mm | 1.06× |

La vertebral es el caso que manda. Por debajo de TR_QISS ≈ 550 ms el inflow vertebral
ya no alcanza a refrescar el slab completo, y por eso el TR no se acorta más pese a que
el tiempo de examen lo agradecería.

### 4.3 Estado de cada tejido en el instante TI

| Tejido | Mz |
|---|---|
| Músculo (fondo) | 0.0 % — nuleado por diseño |
| Grasa | +36.3 % — **no** la nulea el TI |
| Sangre estática / venosa | −4.9 % |
| Sangre arterial entrante | +100 % — nunca vio la inversión |

Que la grasa quede en +36 % es la confirmación de que el nulo del stopband de bSSFP
(sección 3) no es un adorno: es lo único que la suprime.

Que la sangre venosa quede en −4.9 % no es casualidad buscada sino consecuencia de que,
con TR corto, el TI de nulo de un T1 largo baja mucho. Refuerza la banda "tracking".

## 5. Presupuesto del disparo

```
inversión de fondo (8 ms) + inversión venosa (8 ms) + spoiler
  ↓
delay del QI
  ↓
preparación α/2  (TR/2 = 3.14 ms)
  ↓
14 TRs bSSFP × 6.28 ms = 87.9 ms   ← un ángulo radial, todas las particiones
  ↓
relleno hasta completar TR_QISS = 600 ms
```

El centro de kz cae a mitad de ventana con orden ascendente, así que el readout arranca
en TI − 87.9/2 = **162.7 ms** tras la inversión.

**Por qué un solo ángulo radial por disparo.** Koktzoglou usa una ventana de 891 ms con
90 TRs, lo que es viable con FLASH de ángulo bajo. Con bSSFP a 100° una ventana así
satura la sangre entrante y destruye el contraste de inflow, que es el mecanismo entero
de QISS. Se mantiene la ventana corta y el tiempo de examen se recupera bajando el
TR_QISS de 1500 a 600 ms — algo que solo es posible porque, al ser ungated, el TR no
está atado al ritmo cardiaco.

## 6. Geometría y muestreo

| Parámetro | Valor | Origen |
|---|---|---|
| FOV in-plane | 200 mm | cuello |
| Resolución in-plane | 1.2 mm | cedida desde el sub-1 mm de 3T, por SNR |
| Espesor de slab | 19.5 mm | Koktzoglou 2020 |
| Solape de slabs | 4.5 mm | Koktzoglou 2020 |
| Número de slabs | 19 | → 290 mm de cobertura axial |
| Grilla kz | 18 pasos (20 % oversampling) | Koktzoglou 2020 |
| kz adquiridos | 14 (partial Fourier 6/8) | Koktzoglou 2020 |
| Espesor de partición | 1.30 mm | = 23.4 mm / 18 |
| Vistas radiales | 32 (submuestreo 8.2×) | elegido por tiempo de examen |
| Tiempo de examen | 6:04 | 19 × 32 × 600 ms |

Para referencia, el protocolo de 3T cubre 289 mm en 6:39.

### Orden de vistas

El incremento de 15.8245° del paper está afinado para *su* número de vistas. Con las 32
de acá reparte mal: el hueco angular máximo llega a 1.8× el ideal. Opciones medidas en
`verify_qiss3d.py`:

| N | paper 15.8245° | áureo 68.75° | uniforme 180/N |
|---|---|---|---|
| 8 | 207.7 % | 27.9 % | 0.0 % |
| 32 | 75.9 % | 78.3 % | 0.0 % |
| 64 | 40.9 % | 48.0 % | 0.0 % |

(desviación máxima del hueco respecto del reparto ideal)

El uniforme es exacto para cualquier N, pero concentra el aliasing en streaks coherentes.
El script usa **áureo por defecto**: a 8× de submuestreo la reconstrucción va a tener que
ser iterativa/CS —Varghese et al. necesitaron compressed sensing para prácticamente todo
en este scanner— y ahí el aliasing incoherente rinde mejor que el reparto perfecto.
Cambiar `VIEW_ORDER` a `"uniform"` si se reconstruye con gridding/NUFFT simple.

## 7. Lo que queda pendiente antes del scanner

1. ~~Tiempos muertos del sistema.~~ **Confirmados**: `rf_dead_time` = 100 µs,
   `rf_ringdown_time` = 20 µs, `adc_dead_time` = 10 µs. El ringdown resultó 10 µs más
   corto que la estimación conservadora, pero **el TR mínimo no se movió** (sigue en
   5.96 ms): el bloque de excitación lo limita la duración del gradiente selector de
   slab —600 µs de RF más 2 × 160 µs de rampas = 920 µs— y el evento de RF completo
   (100 + 600 + 20 = 720 µs) cabe holgado dentro. Los 10 µs se absorben en la rampa,
   no en el presupuesto de TR. El ancho de banda se mantiene en 356 Hz/px.
2. **PNS.** Los tres ejes cumplen por separado (19.7 / 19.6 / 22.0 mT/m contra el
   límta derateado de 22.1), que es como el fabricante especifica el gradiente. La
   magnitud **vectorial** llega a 28.0 mT/m y 53.3 T/m/s al sumar ejes en diagonal, cosa
   normal en 3D porque los amplificadores son independientes. Pero el PNS sí depende del
   vector y PyPulseq no lo modela. El bore de 80 cm implica una bobina de gradiente
   grande, lo que tiende a empeorar el PNS por unidad de slew: hay que pasar el `.seq`
   por el chequeo del fabricante, no darlo por descontado.
3. **Perfil del pulso adiabático.** El script usa hyperbolic secant, que es lo que trae
   PyPulseq 1.4.2. Los papers usan FOCI, que da bordes de slab más abruptos — y Edelman
   apoya su versión mejorada justamente en esa nitidez para no saturar sangre entrante.
   Conviene simular el perfil (KomaMRI o MRTwin) antes de ir al scanner.
4. **Intérprete Pulseq en el Free.Max.** El sistema corre syngo XA; la versión del
   intérprete condiciona qué features del `.seq` son utilizables.
5. **SAR.** A 0.55 T escala con B0² y no debería limitar ni con FA de 100° ni con los
   dos adiabáticos por disparo, pero conviene confirmarlo en el scanner.
