# Diagrama Hertzsprung-Russell con Gaia DR3

## 1. Descripción

Este proyecto construye un **diagrama Hertzsprung-Russell (HR)** a partir de
datos de la misión **Gaia DR3** de la ESA, obtenidos automáticamente mediante el
servicio **TAP de VizieR**. El objetivo es identificar de forma visual las dos
estructuras fundamentales del diagrama: la **secuencia principal** y la región
de **gigantes rojas**.

Todo el flujo de trabajo es reproducible y automático: un único script
(`pipeline.sh`) descarga los datos, construye una base de datos SQLite y genera
la gráfica final. No se requiere ninguna descarga manual de datos.

## 2. Fundamento astrofísico

Un **diagrama Hertzsprung-Russell** representa la luminosidad de las estrellas
frente a su temperatura superficial (o un indicador de ésta, como el color).
Es la herramienta básica para estudiar la evolución estelar, porque en él las
estrellas no se distribuyen al azar, sino que se agrupan en regiones bien
definidas.

- **Secuencia principal**: banda diagonal donde las estrellas pasan la mayor
  parte de su vida fusionando hidrógeno en helio. Las estrellas más calientes
  (azules) son más luminosas y se sitúan arriba a la izquierda; las más frías
  (rojas) son más débiles y se sitúan abajo a la derecha.
- **Gigantes rojas**: estrellas que han agotado el hidrógeno del núcleo y se
  han expandido y enfriado. Son muy luminosas (pequeño valor de `M_G`) a pesar
  de tener temperaturas superficiales bajas (color rojo), por lo que aparecen
  en la esquina superior derecha del diagrama.

### ¿Por qué `G_BP - G_RP` como indicador de color?

`G_BP` y `G_RP` son las magnitudes de Gaia en las bandas azul y roja. Su
diferencia (`G_BP - G_RP`) mide el **color** de la estrella: cuanto mayor es el
valor, más roja (y, en general, más fría) es la estrella. Sirve por tanto como
indicador de temperatura superficial.

### ¿Por qué `M_G` como magnitud absoluta?

La magnitud aparente `G` depende de la distancia, por lo que no permite comparar
luminosidades intrínsecas. La **magnitud absoluta** `M_G` (la magnitud que
tendría la estrella a 10 pársecs) elimina ese efecto y permite comparar la
luminosidad real de estrellas situadas a distancias distintas.

## 3. Datos

Los datos proceden de **Gaia DR3**, accesibles a través del servicio TAP de
VizieR (catálogo `I/355`, tabla `gaiadr3`):

```
https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync
```

La consulta se realiza en **ADQL**. Se utilizan únicamente las columnas
necesarias:

| Columna          | Descripción                                    |
|------------------|------------------------------------------------|
| `parallax`       | Paralaje (mas)                                 |
| `parallax_error` | Incertidumbre del paralaje (mas)               |
| `phot_g_mean_mag`| Magnitud aparente media en la banda G          |
| `phot_bp_mean_mag` | Magnitud aparente media en la banda BP       |
| `phot_rp_mean_mag` | Magnitud aparente media en la banda RP       |

En el catálogo de VizieR estas columnas se llaman `Plx`, `e_Plx`, `Gmag`,
`BPmag` y `RPmag`. El CSV descargado conserva esos nombres nativos, y
`constructor_db.py` los renombra a los nombres semánticos de la tabla anterior
(`Plx → p`, `e_Plx → p_error`, `Gmag → phot_g_mean_mag`, etc.).

## 4. Selección de datos

La consulta ADQL selecciona `TOP 50000` estrellas e impone las siguientes
condiciones de calidad:

- **`p > 20`** (es decir, `Plx > 20` mas): se conservan estrellas con paralaje
  mayor que 20 mas, lo que corresponde a distancias menores que ~50 pc. Este
  corte selecciona una muestra local casi completa y, además, hace que la
  consulta sobre la tabla completa de Gaia sea rápida.
- **`p / p_error > 10`**: se conservan únicamente las estrellas cuya relación
  señal/ruido del paralaje es alta, es decir, con distancias bien determinadas.

Además se exige que la magnitud G esté disponible (`IS NOT NULL`) y no se
emplea `SELECT *`, sino solo las columnas necesarias.

## 5. Cálculos

A lo largo del proyecto se utiliza la letra **`p`** para representar el
paralaje (en milisegundos de arco, mas). En el código, las columnas nativas de
VizieR se renombran como `Plx → p` y `e_Plx → p_error`.

- **Distancia** (en pársecs):

$$
d[\mathrm{pc}] = \frac{1000}{p[\mathrm{mas}]}
$$

- **Índice de color**:

$$
G_{BP} - G_{RP}
$$

- **Magnitud absoluta** en la banda G:

$$
M_G = G + 5\log_{10}(p) - 10
$$

donde `p` está en mas. Esta es la relación estándar que convierte la magnitud
aparente en absoluta usando el módulo de distancia.

## 6. Ejecución

El proyecto se ejecuta de forma totalmente automática:

```bash
chmod +x pipeline.sh
./pipeline.sh
```

El script comprueba las dependencias, descarga los datos (mostrando un
indicador giratorio de progreso), construye la base de datos SQLite y genera
el diagrama. No es necesario descargar ningún CSV ni editar ningún archivo
manualmente.

### Dependencias

- `wget`
- `Python 3`
- `pandas`, `numpy`, `matplotlib` (ver `requirements.txt`)

Instalación de los paquetes de Python:

```bash
python3 -m pip install -r requirements.txt
```

## 7. Estructura del proyecto

```text
gaia-dr3-hr-diagram/
├── README.md             # Documentación del proyecto
├── pipeline.sh           # Punto de entrada: descarga, construye DB y analiza
├── constructor_db.py     # Lee el CSV y construye la base de datos SQLite
├── analisis_visual.py    # Calcula variables y genera el diagrama HR
├── requirements.txt      # Dependencias de Python
├── .gitignore            # Archivos ignorados por Git
├── datos/
│   └── .gitkeep          # Carpeta para CSV y base de datos (generados)
└── resultados/
    └── .gitkeep          # Carpeta para el diagrama HR (generado)
```

Los datos (`datos/datos_gaia.csv`, `datos/datos_mision.db`) y la figura
(`resultados/diagrama_hr.png`) se generan en cada ejecución y **no** se
versionan en Git.

## 8. Resultado

El resultado final es el diagrama Hertzsprung-Russell:

```text
resultados/diagrama_hr.png
```

La figura se dibuja como **mapa de densidad** (histograma 2D en escala
logarítmica) sobre fondo negro, con cuatro ejes:

- **Inferior**: índice de color `G_BP - G_RP` (indicador de temperatura).
- **Superior**: tipo espectral (O, B, A, F, G, K, M) y temperatura efectiva
  `T_eff` (en K), mediante una calibración color–temperatura aproximada.
- **Izquierdo**: luminosidad en unidades solares `L/L_sol` (escala logarítmica).
- **Derecho**: magnitud absoluta `M_G` (invertida, de modo que las estrellas
  más luminosas quedan arriba).

Resolución: 300 dpi.

## 9. Interpretación

En el diagrama se distinguen claramente:

- **Secuencia principal**: la banda diagonal que cruza el diagrama desde la
  esquina superior izquierda (estrellas calientes y luminosas) hasta la
  esquina inferior derecha (estrellas frías y débiles). Es la etapa más larga
  de la vida estelar (fusión de hidrógeno en el núcleo).
- **Gigantes rojas**: la acumulación de estrellas luminosas y rojas en la
  esquina superior derecha, separada de la secuencia principal. Corresponden a
  estrellas que han agotado el hidrógeno del núcleo y se han expandido.
- **Enanas blancas**: una secuencia casi paralela a la principal en la zona
  inferior izquierda (azul y muy poco luminosa). Es el remanente final de las
  estrellas de masa baja e intermedia, y se hace nítida gracias al mapa de
  densidad. Su presencia cierra el relato evolutivo: de la secuencia principal
  a gigante roja y, finalmente, a enana blanca.

La identificación de estas regiones en el código es **geométrica y
aproximada** (cortes en color y magnitud absoluta), no una clasificación
espectroscópica.

## 10. Limitaciones

La interpretación del diagrama está sujeta a varias limitaciones:

- **Extinción interestelar**: el polvo atenúa y enrojece la luz, desplazando
  las estrellas hacia valores mayores de `G_BP - G_RP` y magnitudes más débiles.
  No se aplica corrección por extinción.
- **Incertidumbre del paralaje**: aunque se exige `p / p_error > 10`, el
  cálculo `d = 1000/p` puede sesgar la distancia (y por tanto `M_G`) para
  paralajes con error relativo no despreciable.
- **Selección de la muestra**: el corte `p > 20` mas limita la muestra a
  estrellas a menos de ~50 pc. Es una muestra local (casi completa en volumen)
  que favorece las estrellas débiles y cercanas, y no incluye las estrellas
  masivas más lejanas; además, no es representativa de toda la Galaxia.
- **Límite de magnitud de Gaia**: las estrellas más débiles no están incluidas,
  lo que recorta la parte inferior de la secuencia principal.
- **`M_G` sin corrección por extinción**: la magnitud absoluta calculada aquí
  no incorpora la extinción, por lo que está sistemáticamente afectada para
  estrellas en regiones con polvo.
