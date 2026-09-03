#!/usr/bin/env python3
"""
- Distancia:            d[pc] = 1000 / p[mas]
- Índice de color:      BP_RP = G_BP - G_RP
- Magnitud absoluta:    M_G = G + 5*log10(p) - 10   (con p en mas)
"""

import argparse
import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

# --- Regiones aproximadas del diagrama HR ------------------------------------
# Gigantes rojas: estrellas rojas (BP_RP grande) y luminosas (M_G pequeña),
GIGANTES_ROJAS = {
    "bp_rp_min": 1.0,
    "bp_rp_max": 3.5,
    "mg_min": -5.0,
    "mg_max": 2.5,
}
# Secuencia principal: de la forma M_G = a * BP_RP + b en el plano (BP_RP, M_G).
SECUENCIA_PRINCIPAL = {
    "a": 4.0,
    "b_brillante": -0.5,  # (estrellas más luminosas)
    "b_debil": 2.5,       # (estrellas más débiles)
    "bp_rp_min": -0.3,
    "bp_rp_max": 3.5,
}

# --- Curvas de clase de luminosidad (clasificación MKK) ------------------------
# Puntos ancla (BP_RP, M_G) de cada clase de luminosidad. Se interpolan con un
# polinomio de bajo grado para obtener curvas suaves. Edita estos valores a
# mano para afinar la posición de cada curva. La clase "V" se deriva de
# SECUENCIA_PRINCIPAL y no aparece aquí.
COLOR_CURVAS = "#e066ff"  # magenta de las curvas y de sus etiquetas

CURVAS_LUMINOSIDAD = {
    "Ia": [(-0.3, -6.0), (0.2, -6.2), (0.7, -6.4), (1.2, -6.5), (1.7, -6.7), (2.2, -6.9)],
    "Ib": [(-0.3, -4.2), (0.2, -4.4), (0.7, -4.6), (1.2, -4.8), (1.7, -5.0), (2.2, -5.2)],
    "II": [(-0.2, -2.2), (0.3, -2.4), (0.7, -2.6), (1.1, -2.8), (1.5, -3.0), (2.0, -3.3)],
    "III": [(0.5, 1.0), (0.8, 0.7), (1.0, 0.5), (1.4, 0.1), (1.8, -0.6), (2.2, -1.4)],
    "IV": [(0.1, 0.3), (0.4, 1.2), (0.7, 2.2), (1.0, 3.2), (1.3, 4.2), (1.6, 5.2)],
    "Enanas blancas": [(-0.3, 10.5), (-0.1, 11.3), (0.1, 12.2), (0.3, 13.3), (0.5, 14.5), (0.7, 15.6)],
}

# Orden de dibujado de las clases de luminosidad.
CLASES_LUMINOSIDAD = ["Ia", "Ib", "II", "III", "IV", "V", "Enanas blancas"]

# --- Ejes de luminosidad y magnitud absoluta ---------------------------------
# Relación estándar entre luminosidad y magnitud absoluta, tomando como
# L / L_sol = 10^{-0.4 * (M - M_sol)}
MAG_ABS_SOL = 4.74
LUMINOSIDADES = [1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

# --- Tipo espectral ------------------
CLASES_ESPECTRALES = [
    ("O", -0.35, "#9bb0ff"),
    ("B", -0.10, "#aabfff"),
    ("A",  0.10, "#f2f2ff"),
    ("F",  0.35, "#fff4e6"),
    ("G",  0.60, "#ffd9a3"),
    ("K",  1.00, "#ffa85c"),
    ("M",  1.80, "#ff6a3d"),
]

# Colormap continuo construido a partir de los colores de CLASES_ESPECTRALES.
# Se interpolan los 7 colores en la posición BP_RP de cada clase (O azul -> M
# rojo). Fuera del rango [-0.35, 1.80] el color queda fijo en el extremo.
_BP_MIN_ANCLA = CLASES_ESPECTRALES[0][1]
_BP_MAX_ANCLA = CLASES_ESPECTRALES[-1][1]
CMAP_ESPECTRAL = LinearSegmentedColormap.from_list(
    "espectral",
    [
        ((bp - _BP_MIN_ANCLA) / (_BP_MAX_ANCLA - _BP_MIN_ANCLA), color)
        for _, bp, color in CLASES_ESPECTRALES
    ],
)

# --- Límites y estilo del diagrama --------------------------------------------
X_MIN, X_MAX = -0.6, 4.0        # eje de color (BP_RP)
MG_MIN, MG_MAX = -10.0, 18.0    # magnitud absoluta (arriba y abajo)

COLOR_FONDO = "black"
COLOR_TEXTO = "white"

# --- Calibración aproximada color (BP_RP) <-> temperatura efectiva (Teff) ------
# Puntos de la relación color-temperatura para la secuencia principal. Sirven
# para situar las marcas del eje superior de temperatura. Se eligen pocos ticks
# bien repartidos para que sus etiquetas no se solapen en el rango BP_RP.
TEFF_TICKS = [30000, 12000, 8000, 6000, 5000, 4000, 3000]

TEFF_BP_RP = [
    (-0.50, 40000.0),
    (-0.25, 20000.0),
    (-0.05, 10000.0),
    (0.15, 7500.0),
    (0.40, 6000.0),
    (0.65, 5200.0),
    (1.00, 4400.0),
    (1.40, 3600.0),
    (2.00, 2900.0),
    (3.00, 2300.0),
    (4.00, 1900.0),
]
_BP_TEFF = [bp for bp, _ in TEFF_BP_RP]
_TEFF_VAL = [teff for _, teff in TEFF_BP_RP]


def teff_a_bp_rp(teff):
    """Convierte temperatura efectiva (K) en índice de color BP_RP."""
    return np.interp(teff, _TEFF_VAL[::-1], _BP_TEFF[::-1])


def luminosidad_a_mg(lum):
    """Convierte luminosidad relativa al Sol en magnitud absoluta."""
    return MAG_ABS_SOL - 2.5 * np.log10(lum)


def formatear_luminosidad(lum):
    """Devuelve una etiqueta para una luminosidad en potencias de 10."""
    exp = int(round(np.log10(lum)))
    if exp == 0:
        return "1"
    return rf"$10^{{{exp}}}$"


def parse_args():
    parser = argparse.ArgumentParser(description="Genera el diagrama HR.")
    parser.add_argument("--db", required=True, help="Ruta de la base SQLite.")
    parser.add_argument("--out", required=True, help="Ruta del PNG de salida.")
    return parser.parse_args()


def leer_datos(ruta_db):
    """Lee la tabla `estrellas` desde la base de datos SQLite."""
    ruta = Path(ruta_db)
    if not ruta.exists():
        raise FileNotFoundError(f"No se encontró la base de datos: {ruta}")
    with sqlite3.connect(ruta) as conexion:
        df = pd.read_sql_query("SELECT * FROM estrellas", conexion)
    if df.empty:
        raise ValueError("La tabla 'estrellas' está vacía.")
    return df


def calcular_variables(df):
    """Calcula distancia, color y magnitud absoluta a partir del paralaje `p`."""
    df = df.copy()

    # Distancia en pársecs a partir del paralaje en milisegundos de arco.
    df["distancia_pc"] = 1000.0 / df["p"]

    # Índice de color.
    df["BP_RP"] = df["phot_bp_mean_mag"] - df["phot_rp_mean_mag"]

    # Magnitud absoluta en la banda G (p en mas).
    df["M_G"] = df["phot_g_mean_mag"] + 5.0 * np.log10(df["p"]) - 10.0

    return df


def mascara_gigantes_rojas(bp_rp, mg):
    """Máscara booleana para la región aproximada de gigantes rojas."""
    return (
        (bp_rp >= GIGANTES_ROJAS["bp_rp_min"])
        & (bp_rp <= GIGANTES_ROJAS["bp_rp_max"])
        & (mg >= GIGANTES_ROJAS["mg_min"])
        & (mg <= GIGANTES_ROJAS["mg_max"])
    )


def mascara_secuencia_principal(bp_rp, mg):
    """Máscara booleana para la banda aproximada de la secuencia principal."""
    a = SECUENCIA_PRINCIPAL["a"]
    borde_brillante = a * bp_rp + SECUENCIA_PRINCIPAL["b_brillante"]
    borde_debil = a * bp_rp + SECUENCIA_PRINCIPAL["b_debil"]
    return (
        (bp_rp >= SECUENCIA_PRINCIPAL["bp_rp_min"])
        & (bp_rp <= SECUENCIA_PRINCIPAL["bp_rp_max"])
        & (mg >= borde_brillante)
        & (mg <= borde_debil)
    )


def curva_secuencia_principal():
    """Curva de la clase V (secuencia principal) derivada de SECUENCIA_PRINCIPAL."""
    a = SECUENCIA_PRINCIPAL["a"]
    b = (SECUENCIA_PRINCIPAL["b_brillante"] + SECUENCIA_PRINCIPAL["b_debil"]) / 2
    x = np.linspace(SECUENCIA_PRINCIPAL["bp_rp_min"], SECUENCIA_PRINCIPAL["bp_rp_max"], 300)
    return x, a * x + b


def ajustar_curva(puntos, grado=2):
    """Ajusta un polinomio de bajo grado a los puntos ancla (BP_RP, M_G)."""
    xs = np.array([p[0] for p in puntos], dtype=float)
    ys = np.array([p[1] for p in puntos], dtype=float)
    coefs = np.polyfit(xs, ys, grado)
    x_suave = np.linspace(xs.min(), xs.max(), 300)
    return x_suave, np.polyval(coefs, x_suave)


def mostrar_estadisticas(df):
    """Imprime estadísticas básicas de la muestra."""
    total = len(df)
    rango_color = (df["BP_RP"].min(), df["BP_RP"].max())
    rango_mg = (df["M_G"].min(), df["M_G"].max())
    distancias = df["distancia_pc"]

    n_gigantes = int(mascara_gigantes_rojas(df["BP_RP"], df["M_G"]).sum())
    n_ms = int(mascara_secuencia_principal(df["BP_RP"], df["M_G"]).sum())

    print("=== Estadísticas de la muestra ===")
    print(f"  Estrellas totales:                 {total:,}")
    print(f"  Rango de color (BP_RP):            [{rango_color[0]:.2f}, {rango_color[1]:.2f}] mag")
    print(f"  Rango de M_G:                      [{rango_mg[0]:.2f}, {rango_mg[1]:.2f}] mag")
    print(f"  Distancia mínima:                  {distancias.min():.1f} pc")
    print(f"  Distancia mediana (P50):           {distancias.median():.1f} pc")
    print(f"  Distancia máxima:                  {distancias.max():.1f} pc")
    print(f"  Percentil 5 / 95 de distancia:     {distancias.quantile(0.05):.1f} / {distancias.quantile(0.95):.1f} pc")
    print()
    print(f"  Estrellas en la secuencia principal: {n_ms:,}")
    print(f"  Estrellas en la región de gigantes:  {n_gigantes:,}")


def dibujar_regiones(ax):
    """Dibuja las curvas de clase de luminosidad (Ia, Ib, II, III, IV, V, enanas blancas)."""
    for nombre in CLASES_LUMINOSIDAD:
        if nombre == "V":
            x, y = curva_secuencia_principal()
        else:
            x, y = ajustar_curva(CURVAS_LUMINOSIDAD[nombre])

        ax.plot(x, y, color=COLOR_CURVAS, linewidth=1.2, zorder=3)
        # Etiqueta pegada al extremo derecho de cada curva (en vez de leyenda).
        ax.text(
            x[-1] + 0.05,
            y[-1],
            nombre,
            color=COLOR_CURVAS,
            fontsize=9,
            va="center",
            ha="left",
            zorder=4,
        )


def dibujar_densidad(ax, bp_rp, mg):
    """Dibuja el diagrama HR como campo estelar: puntos coloreados por temperatura.

    Cada estrella es un punto cuyo color depende del índice de color BP_RP
    (análogo aproximado del color Johnson B-V, no el mismo índice) a través del
    colormap continuo CMAP_ESPECTRAL. La transparencia hace que la densidad se
    aprecie por superposición de puntos, como en un campo estelar real.
    """
    ax.scatter(
        bp_rp,
        mg,
        s=2,
        c=bp_rp,
        cmap=CMAP_ESPECTRAL,
        vmin=_BP_MIN_ANCLA,
        vmax=_BP_MAX_ANCLA,
        alpha=0.7,
        edgecolors="none",
        linewidths=0,
        rasterized=True,
        zorder=1,
    )


def dibujar_ejes_superiores(ax):
    """Añade dos ejes superiores: tipo espectral (letras) y temperatura (Teff)."""
    # Eje de tipo espectral: letras coloreadas, con un pequeño margen sobre el
    # borde superior para que no se peguen a la figura.
    ax_spec = ax.twiny()
    ax_spec.set_xlim(ax.get_xlim())
    ax_spec.spines["top"].set_position(("axes", 1.03))
    ax_spec.spines["top"].set_color(COLOR_TEXTO)
    posiciones = [bp_rp for _, bp_rp, _ in CLASES_ESPECTRALES]
    etiquetas = [letra for letra, _, _ in CLASES_ESPECTRALES]
    ax_spec.set_xticks(posiciones)
    ax_spec.set_xticklabels(etiquetas, fontsize=12, fontweight="bold")
    for etiqueta, (_, _, color) in zip(ax_spec.get_xticklabels(), CLASES_ESPECTRALES):
        etiqueta.set_color(color)
    ax_spec.tick_params(axis="x", length=0)

    # Eje de temperatura efectiva, más separado del de tipo espectral. Su
    # etiqueta se coloca por encima de las marcas para no chocar con las letras.
    ax_teff = ax.twiny()
    ax_teff.set_xlim(ax.get_xlim())
    ax_teff.spines["top"].set_position(("axes", 1.14))
    ax_teff.set_xticks([teff_a_bp_rp(t) for t in TEFF_TICKS])
    ax_teff.set_xticklabels([str(t) for t in TEFF_TICKS], fontsize=8)
    ax_teff.set_xlabel(r"$T_{\mathrm{eff}}$  [K]", color=COLOR_TEXTO, fontsize=9)
    ax_teff.xaxis.set_label_coords(0.5, 1.22)
    ax_teff.tick_params(axis="x", colors=COLOR_TEXTO)
    ax_teff.spines["top"].set_color(COLOR_TEXTO)

    return ax_spec, ax_teff


def generar_diagrama(df, ruta_salida):
    """Genera y guarda el diagrama HR.

    Configuración:
        - Fondo negro y campo estelar (puntos coloreados por temperatura).
        - Eje inferior: color (G_BP - G_RP).
        - Eje superior: tipo espectral (O..M) y temperatura efectiva (Teff).
        - Eje izquierdo: luminosidad en unidades solares (L / L_sol).
        - Eje derecho: magnitud absoluta (M_G).
    """
    fig, ax = plt.subplots(figsize=(9, 8))

    # Fondo negro.
    fig.patch.set_facecolor(COLOR_FONDO)
    ax.set_facecolor(COLOR_FONDO)

    # Campo estelar del diagrama HR.
    dibujar_densidad(ax, df["BP_RP"], df["M_G"])
    dibujar_regiones(ax)

    # Límites de los ejes. El eje Y queda invertido (luminosas arriba).
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(MG_MAX, MG_MIN)

    # Eje izquierdo: luminosidades solares (escala logarítmica).
    pos_mg = luminosidad_a_mg(np.array(LUMINOSIDADES))
    ax.set_yticks(pos_mg)
    ax.set_yticklabels([formatear_luminosidad(l) for l in LUMINOSIDADES])
    ax.set_ylabel(r"Luminosidad  $L\,/\,L_\odot$")

    # Eje derecho: magnitud absoluta.
    ax_mag = ax.twinx()
    ax_mag.set_ylim(ax.get_ylim())
    ax_mag.set_yticks(np.arange(-10, 20, 5))
    ax_mag.set_ylabel(r"$M_G$  [mag]")

    # Eje inferior: índice de color (análogo aproximado de B-V, no el mismo
    # índice; BP_RP es la fotometría de Gaia).
    ax.set_xlabel(r"Color  $G_{BP} - G_{RP}$  [mag]")

    # Ejes superiores: tipo espectral y temperatura efectiva.
    ax_spec, ax_teff = dibujar_ejes_superiores(ax)
    ax_spec.set_xlabel("Tipo espectral", color=COLOR_TEXTO, fontsize=10)
    ax_spec.xaxis.set_label_coords(0.5, 1.11)

    # Estilo general (texto y espinas en blanco sobre fondo negro).
    fig.suptitle("Diagrama Hertzsprung-Russell — Gaia DR3", color=COLOR_TEXTO, fontsize=13, y=0.96)
    ax.grid(True, alpha=0.15, linewidth=0.5, color=COLOR_TEXTO)
    ax.tick_params(colors=COLOR_TEXTO)
    ax_mag.tick_params(colors=COLOR_TEXTO)
    ax.xaxis.label.set_color(COLOR_TEXTO)
    ax.yaxis.label.set_color(COLOR_TEXTO)
    ax_mag.yaxis.label.set_color(COLOR_TEXTO)
    for espina in ("bottom", "left"):
        ax.spines[espina].set_color(COLOR_TEXTO)
    ax_mag.spines["right"].set_color(COLOR_TEXTO)

    # Deja margen superior para el título y los dos ejes superiores.
    fig.subplots_adjust(left=0.14, right=0.87, bottom=0.11, top=0.76)

    ruta = Path(ruta_salida)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta, dpi=300, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    args = parse_args()

    print("Leyendo datos desde SQLite...")
    df = leer_datos(args.db)

    print("Calculando distancia, color y magnitud absoluta...")
    df = calcular_variables(df)

    mostrar_estadisticas(df)

    print("Generando el diagrama HR...")
    generar_diagrama(df, args.out)
    print(f"Diagrama guardado en '{args.out}'.")


if __name__ == "__main__":
    main()
