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

## 3. TR de bSSFP: por qué 6.28 ms

> **Corrección.** Este apartado sostenía que el TR de 6.28 ms suprime la grasa por sí
> solo al colocarla en el nulo del stopband. **La simulación lo refutó** — ver sección 9.
> El TR se mantiene en 6.28 ms porque es aproximadamente el mínimo que permiten estos
> gradientes, pero la supresión grasa la hace un pulso espectral dedicado.

### El argumento original (que resultó insuficiente)

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

Que la grasa quede en +36 % es lo que obliga a un fat-sat espectral dedicado: ni el TI
ni el nulo del stopband la bajan lo suficiente en una ventana de 14 TRs (sección 9).

Que la sangre venosa quede en −4.9 % no es casualidad buscada sino consecuencia de que,
con TR corto, el TI de nulo de un T1 largo baja mucho. Refuerza la banda "tracking".

## 5. Presupuesto del disparo

```
inversión de fondo WURST (16 ms) + inversión venosa WURST (16 ms) + spoiler
  ↓
delay del QI
  ↓
fat-sat espectral (15 ms) + spoiler      ← pegado al tren, para que la grasa no se recupere
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
3. ~~Perfil del pulso adiabático.~~ **Simulado y corregido** — ver sección 8.
4. ~~Intérprete Pulseq en el Free.Max.~~ **Confirmado: 1.4.2**, que coincide
   exactamente con la versión de PyPulseq usada, así que el `.seq` se escribe en
   formato 1.4.2 y no hay que degradar ninguna feature. La secuencia solo usa RF
   arbitrario (los adiabáticos), trapecios, ADC y delays — nada de labels, triggers
   ni extensions, que es donde más divergen las versiones de intérprete.
   Lo que **sí** queda por comprobar es que el intérprete cargue los **38 304
   bloques**: el `.seq` se carga entero en memoria. Conviene probar primero con
   `--quick` (2 slabs x 4 vistas, 504 bloques) y confirmar que carga antes del examen completo.
5. **SAR.** A 0.55 T escala con B0² y no debería limitar ni con FA de 100° ni con los
   dos adiabáticos por disparo, pero conviene confirmarlo en el scanner.


## 8. Perfil de los pulsos adiabáticos (simulación de Bloch)

Los papers usan FOCI; PyPulseq 1.4.2 no lo trae. En vez de asumir que un sustituto
sirve, se midió el perfil por integración directa de la ecuación de Bloch sobre la
forma de onda compleja real del pulso (`sim_inversion_profile.py`). Sin aproximación
de ángulo pequeño: rotación de Rodrigues en torno al campo efectivo, paso de 1 µs.

Importa porque QISS depende de que la sangre arterial que sube no haya visto la
inversión. Un borde inferior blando entrega el bolo con la cabecera oscura.

### El primer intento falló

Hyperbolic secant de 8 ms, banda venosa de 100 mm como el paper:

| | transición | fuga bajo el borde nominal | veredicto |
|---|---|---|---|
| Inversión de fondo (19.5 mm) | 4.5 mm | 2.25 mm | tolerable (5–11 % del bolo) |
| Banda venosa (100 mm) | 22.0 mm | 11.0 mm | **invade el slab 6 mm** |

La causa es geométrica, no del pulso. El ancho de transición **en Hz** lo fija el
pulso; al pasarlo a milímetros se divide por el gradiente selector. Una banda 5 veces
más gruesa que el slab necesita un gradiente 5 veces más débil (0.293 contra 1.50
mT/m), así que la misma transición se estira 5 veces. Con solo 5 mm de hueco, el borde
blando de la banda venosa suprimía los 6 mm superiores del slab de imagen.

### Dos correcciones

**Adelgazar la banda venosa de 100 a 60 mm.** La banda solo tiene que cubrir lo que la
sangre venosa alcanza a bajar durante el QI: 41 mm a 20 cm/s, 62 mm a 30 cm/s. Los
100 mm del paper responden a su QI de 583 ms (117 mm a 20 cm/s). Copiarlos con nuestro
QI de 207 ms era sobredimensionar, y el precio no era tiempo sino nitidez de borde.

**Cambiar hypsec por WURST de 16 ms.** Barrido en `sim_venous_band_fix.py`: el hypsec
satura en ~7.5 mm de fuga por más que se alargue, mientras el WURST sigue afilando.
Además exige **menos** B1 de pico.

### Resultado

| | transición | fuga | impacto |
|---|---|---|---|
| Inversión de fondo | 2.75 mm | 0.25 mm | 1 % del bolo, carótida y vertebral |
| Banda venosa (60 mm) | 8.5 mm | 1.0 mm | no invade; sobran 4 mm de margen |

B1 de pico exigido: **9.37 µT** (contra 13.24 µT del hypsec).

**Pendiente:** confirmar el B1 máximo del cuerpo del Free.Max. Un adiabático por
debajo de su umbral de B1 deja de invertir de forma uniforme y todo este análisis
se cae. Es el único número de esta sección que no está verificado.


## 9. Contraste del tren bSSFP (simulación) y la corrección de la supresión grasa

`sim_bssfp_contrast.py` simula el tren completo —estado post-inversión, preparación
α/2, 14 TRs con alternancia de fase, muestreo en TE = TR/2— para cada tejido.

### Lo que refutó

El nulo del stopband de bSSFP es un fenómeno de **estado estacionario**. La ventana de
lectura de QISS son 14 TRs, 88 ms. Con T1/T2 de grasa de 187/93 ms eso no alcanza:

| | \|Mxy\| |
|---|---|
| Grasa en estado estacionario (400 TRs) | 0.014 |
| Grasa en el centro de k (TR #10), sin fat-sat | 0.189 |

Trece veces y media por encima del asintótico. El contraste arteria/grasa quedaba en
**3.6:1**, con la grasa brillante en las MIP.

Se descartaron por medición dos alternativas más baratas: mover kz=0 al final del tren
sube el contraste solo a 5.1:1, y bajar el flip angle cuesta más señal arterial de la
que gana en supresión (a 50° el contraste incluso empeora, 2.1:1).

### La corrección

Pulso gaussiano de 90°, 15 ms, TBW 1.0, centrado a −79.6 Hz, **selectivo en frecuencia
y no en espacio** (sin gradiente durante el RF, así que no puede saturar sangre
entrante por su posición). Va pegado al tren para que la grasa no se recupere.

La duración se eligió simulando el perfil espectral, no por defecto:

| duración | Mz agua | Mz grasa | Mz agua peor caso a ±20 Hz |
|---|---|---|---|
| 10 ms | +0.878 | +0.047 | +0.660 |
| **15 ms** | **+0.996** | **+0.070** | **+0.939** |
| 20 ms | +0.982 | +0.092 | +0.994 (pero grasa sube a +0.413) |

B1 de pico: 0.496 µT, despreciable frente a los 9.37 µT de los adiabáticos.

### Contraste resultante en el centro de k

| tejido | \|Mxy\| | contra arteria |
|---|---|---|
| Sangre arterial entrante | 0.676 | — |
| Sangre estática / venosa | 0.007 | 104 : 1 |
| Grasa con fat-sat | 0.031 | 22 : 1 |
| Músculo | 0.049 | 14 : 1 |
| *Grasa sin fat-sat* | *0.189* | *3.6 : 1* |

El músculo pasa a ser el fondo limitante, no la grasa.

**Límite conocido:** con una desviación de B0 de 20 Hz la grasa sube a 0.140 y el
contraste cae a 4.8:1. A 0.55T la inhomogeneidad en Hz es pequeña, pero conviene un
buen shim sobre el cuello. Es el punto más frágil del diseño.


## 10. Cómo averiguar el B1 máximo del Free.Max

Es el único límite de hardware que PyPulseq **no** modela: `pp.Opts` no tiene campo de
B1, así que un pulso que exceda la bobina no falla al construir el `.seq` — lo rechaza
el intérprete en el scanner. Por eso el script hace una auditoría al exportar.

### Exigencia de la secuencia

| Pulso | B1 de pico |
|---|---|
| Excitación bSSFP (100°, 600 µs) | **33.36 µT** ← el que manda |
| Inversiones WURST (16 ms) | 9.37 µT |
| Fat-sat espectral (90°, 15 ms) | 0.50 µT |

**Piso duro: 8.43 µT.** Por debajo, los WURST dejan de invertir de forma uniforme
(eficiencia < 0.95) y el contraste de QISS se cae. No es una degradación suave.

### Ruta 1 — calcularlo desde la reference voltage (la más rápida)

Siemens define la *transmitter reference amplitude* como la tensión que produce un
pulso rectangular de 180° en 1 ms. Eso son exactamente 11.74 µT:

```
180° = 360 · γ̄ · B1 · T  →  B1 = 180 / (360 × 42.576e6 × 1e-3) = 11.74 µT
```

De ahí:

```
B1_max [µT] = 11.74 × V_amplificador_max / V_referencia
```

La reference voltage aparece tras la calibración de transmisor en la tarjeta de
ajustes del protocolo. La tensión máxima del amplificador es una constante del sistema
que hay que pedir a Siemens o al especialista de aplicaciones.

### Ruta 2 — empíricamente con Pulseq

Escribir un `.seq` mínimo con un pulso rectangular de amplitud creciente y ver a partir
de qué valor el intérprete lo rechaza o lo escala. Es la medida directa del límite
efectivo, incluyendo cualquier derate que aplique el sistema. Tiene la ventaja de medir
lo que realmente se puede tocar, no lo que dice la especificación.

### Ruta 3 — preguntar

Al especialista de aplicaciones de Siemens o en la comunidad de Pulseq. El dato de
bobina de cuerpo suele estar documentado en el archivo de configuración del sistema.

### Qué hacer con la respuesta

El B1 de pico escala aproximadamente como 1/duración del pulso, así que la excitación
bSSFP se puede alargar. El coste es indirecto: el bloque de excitación crece y, para
mantener el TR de 6.28 ms, el readout se acorta y sube el ancho de banda, lo que resta
SNR justo donde menos sobra.

| B1 máx | duración mínima del RF a 100° | ancho de banda resultante | veredicto |
|---|---|---|---|
| 8 µT | — | — | **inviable**: ni alcanza el piso adiabático |
| 10 µT | 2000 µs | 629 Hz/px | funciona, ~1.3× menos SNR |
| 15 µT | 1334 µs | ~450 Hz/px | funciona, penalización moderada |
| 20 µT | 1000 µs | 394 Hz/px | casi sin coste |
| ≥ 33 µT | 600 µs | 356 Hz/px | el diseño actual, sin tocar nada |

Cuando tengas el número, ponerlo en `MAX_B1_UT` y la auditoría avisa si algún pulso lo
excede.

### Margen adiabático: ajustado

Los WURST están diseñados a 9.37 µT sobre un umbral de 8.43 µT — solo **1.1× de
margen**. Eso es poco para un adiabático, cuya gracia es precisamente ser insensible al
B1 *por encima* del umbral. Si el B1 varía por el cuello (y con un bore de 80 cm varía),
algunas regiones pueden caer por debajo. Si el sistema da holgura, conviene subir la
amplitud de diseño de los adiabáticos para separarse del umbral.
