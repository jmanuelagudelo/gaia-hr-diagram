#!/usr/bin/env python3
"""
- Distancia:            d[pc] = 1000 / p[mas]
- Índice de color:      BP_RP = G_BP - G_RP
- Magnitud absoluta:    M_G = G + 5*log10(p) - 10   (con p en mas)

Genera un diagrama HR e identifica de forma aproximada la secuencia principal y la
región de gigantes rojas.
"""

import argparse
import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon, Rectangle

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

# --- Ejes de luminosidad y magnitud absoluta ---------------------------------
# Relación estándar entre luminosidad y magnitud absoluta, tomando como
# L / L_sol = 10^{-0.4 * (M - M_sol)}
MAG_ABS_SOL = 4.74
LUMINOSIDADES = [1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

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

# --- Límites y estilo del diagrama --------------------------------------------
X_MIN, X_MAX = -0.6, 4.0        # eje de color (BP_RP)
MG_MIN, MG_MAX = -6.0, 18.0     # magnitud absoluta (arriba y abajo)

COLOR_FONDO = "black"
COLOR_TEXTO = "white"

# Número de celdas (bins) del mapa de densidad del diagrama HR.
BINS_DENSIDAD = (600, 600)

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


def poligono_secuencia_principal():
    """Devuelve los vértices (x, y) del polígono de la banda de secuencia principal."""
    a = SECUENCIA_PRINCIPAL["a"]
    x = np.linspace(
        SECUENCIA_PRINCIPAL["bp_rp_min"], SECUENCIA_PRINCIPAL["bp_rp_max"], 100
    )
    y_brillante = a * x + SECUENCIA_PRINCIPAL["b_brillante"]
    y_debil = a * x + SECUENCIA_PRINCIPAL["b_debil"]
    xs = np.concatenate([x, x[::-1]])
    ys = np.concatenate([y_brillante, y_debil[::-1]])
    return xs, ys


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
    print("=== Regiones aproximadas del diagrama (cortes geométricos) ===")
    print(f"  Estrellas en la secuencia principal: {n_ms:,}")
    print(f"  Estrellas en la región de gigantes:  {n_gigantes:,}")
    print("  Nota: NO es una clasificación espectroscópica.")


def dibujar_regiones(ax):
    """Delimita (con contornos) las regiones de secuencia principal y gigantes rojas."""
    g = GIGANTES_ROJAS
    ax.add_patch(
        Rectangle(
            (g["bp_rp_min"], g["mg_min"]),
            g["bp_rp_max"] - g["bp_rp_min"],
            g["mg_max"] - g["mg_min"],
            facecolor="none",
            edgecolor="tab:red",
            linewidth=2.0,
            linestyle="--",
            alpha=0.9,
            label="Gigantes rojas (aprox.)",
        )
    )

    xs, ys = poligono_secuencia_principal()
    ax.add_patch(
        Polygon(
            np.column_stack([xs, ys]),
            closed=True,
            facecolor="none",
            edgecolor="tab:blue",
            linewidth=2.0,
            linestyle="--",
            alpha=0.9,
            label="Secuencia principal (aprox.)",
        )
    )


def dibujar_densidad(ax, bp_rp, mg):
    """Dibuja el diagrama HR como mapa de densidad (histograma 2D en escala log)."""
    hist, x_edges, y_edges = np.histogram2d(
        bp_rp,
        mg,
        bins=BINS_DENSIDAD,
        range=[[X_MIN, X_MAX], [MG_MIN, MG_MAX]],
    )
    norm = LogNorm(vmin=1.0, vmax=max(hist.max(), 1.0))
    ax.pcolormesh(
        x_edges,
        y_edges,
        hist.T,
        cmap="inferno",
        norm=norm,
        rasterized=True,
        shading="auto",
    )


def dibujar_ejes_superiores(ax):
    """Añade dos ejes superiores: tipo espectral (letras) y temperatura (Teff)."""
    # Eje de tipo espectral: letras coloreadas, pegado al borde superior.
    ax_spec = ax.twiny()
    ax_spec.set_xlim(ax.get_xlim())
    posiciones = [bp_rp for _, bp_rp, _ in CLASES_ESPECTRALES]
    etiquetas = [letra for letra, _, _ in CLASES_ESPECTRALES]
    ax_spec.set_xticks(posiciones)
    ax_spec.set_xticklabels(etiquetas, fontsize=12, fontweight="bold")
    for etiqueta, (_, _, color) in zip(ax_spec.get_xticklabels(), CLASES_ESPECTRALES):
        etiqueta.set_color(color)
    ax_spec.tick_params(axis="x", length=0)
    ax_spec.spines["top"].set_color(COLOR_TEXTO)

    # Eje de temperatura efectiva, desplazado por encima del de tipo espectral.
    ax_teff = ax.twiny()
    ax_teff.set_xlim(ax.get_xlim())
    ax_teff.spines["top"].set_position(("axes", 1.10))
    ax_teff.set_xticks([teff_a_bp_rp(t) for t in TEFF_TICKS])
    ax_teff.set_xticklabels([str(t) for t in TEFF_TICKS], fontsize=9)
    ax_teff.set_xlabel(r"$T_{\mathrm{eff}}$  [K]", color=COLOR_TEXTO, fontsize=10)
    ax_teff.tick_params(axis="x", colors=COLOR_TEXTO)
    ax_teff.spines["top"].set_color(COLOR_TEXTO)

    return ax_spec, ax_teff


def generar_diagrama(df, ruta_salida):
    """Genera y guarda el diagrama HR.

    Configuración:
        - Fondo negro y mapa de densidad (histograma 2D en escala log).
        - Eje inferior: color (G_BP - G_RP).
        - Eje superior: tipo espectral (O..M) y temperatura efectiva (Teff).
        - Eje izquierdo: luminosidad en unidades solares (L / L_sol).
        - Eje derecho: magnitud absoluta (M_G).
    """
    fig, ax = plt.subplots(figsize=(9, 8))

    # Fondo negro.
    fig.patch.set_facecolor(COLOR_FONDO)
    ax.set_facecolor(COLOR_FONDO)

    # Mapa de densidad del diagrama HR.
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
    ax_mag.set_yticks(np.arange(-5, 20, 5))
    ax_mag.set_ylabel(r"$M_G$  [mag]")

    # Eje inferior: índice de color.
    ax.set_xlabel(r"Color  $G_{BP} - G_{RP}$  [mag]")

    # Ejes superiores: tipo espectral y temperatura efectiva.
    dibujar_ejes_superiores(ax)

    # Estilo general (texto y espinas en blanco sobre fondo negro).
    fig.suptitle("Diagrama Hertzsprung-Russell — Gaia DR3", color=COLOR_TEXTO, fontsize=13)
    ax.grid(True, alpha=0.15, linewidth=0.5, color=COLOR_TEXTO)
    ax.tick_params(colors=COLOR_TEXTO)
    ax_mag.tick_params(colors=COLOR_TEXTO)
    ax.xaxis.label.set_color(COLOR_TEXTO)
    ax.yaxis.label.set_color(COLOR_TEXTO)
    ax_mag.yaxis.label.set_color(COLOR_TEXTO)
    for espina in ("bottom", "left"):
        ax.spines[espina].set_color(COLOR_TEXTO)
    ax_mag.spines["right"].set_color(COLOR_TEXTO)

    leyenda = ax.legend(
        loc="lower left",
        fontsize=9,
        facecolor=COLOR_FONDO,
        edgecolor=COLOR_TEXTO,
        framealpha=0.9,
    )
    for texto in leyenda.get_texts():
        texto.set_color(COLOR_TEXTO)

    # Deja margen superior para el título y los dos ejes superiores.
    fig.subplots_adjust(left=0.14, right=0.87, bottom=0.11, top=0.84)

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
