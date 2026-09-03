import argparse
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

Req_Cols = ["Plx", "e_Plx", "Gmag", "BPmag", "RPmag"]
col_names= {
    "Plx": "p",
    "e_Plx": "p_error",
    "Gmag": "phot_g_mean_mag",
    "BPmag": "phot_bp_mean_mag",
    "RPmag": "phot_rp_mean_mag",
}
SIGNAL_TO_NOISE_MIN = 10.0


def parse_args():
    parser = argparse.ArgumentParser(description="Construye la base de datos SQLite.")
    parser.add_argument("--csv", required=True, help="Ruta del CSV descargado.")
    parser.add_argument("--db", required=True, help="Ruta de salida de la base SQLite.")
    return parser.parse_args()


def leer_csv(ruta_csv):
    """Lee el CSV y devuelve un DataFrame con nombres de columna sin espacios."""
    ruta = Path(ruta_csv)
    if not ruta.exists():
        raise FileNotFoundError(f"No se encontró el archivo CSV: {ruta}")

    df = pd.read_csv(ruta, comment="#")
    df.columns = [col.strip() for col in df.columns]
    return df

def limpiar_datos(df):
    """Elimina filas con valores faltantes y aplica validaciones básicas."""
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=Req_Cols).copy()

    for columna in Req_Cols:
        df[columna] = pd.to_numeric(df[columna], errors="coerce")
    df = df.dropna(subset=Req_Cols)


    df = df[df["Plx"] > 0]
    df = df[(df["Plx"] / df["e_Plx"]) > SIGNAL_TO_NOISE_MIN]

    df = df.rename(columns=col_names)
    return df.reset_index(drop=True)


def crear_base_datos(df, ruta_db):
    """Crea la base de datos SQLite e inserta los datos en la tabla `estrellas`."""
    ruta = Path(ruta_db)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    if ruta.exists():
        ruta.unlink()

    with sqlite3.connect(ruta) as conexion:
        df.to_sql("estrellas", conexion, index=False, if_exists="replace")
    return ruta


def main():
    args = parse_args()

    print("Leyendo el CSV descargado...")
    df = leer_csv(args.csv)
    print(f"  Filas leídas del CSV: {len(df):,}")

    df = limpiar_datos(df)
    print(f"  Filas tras limpieza y validación: {len(df):,}")

    crear_base_datos(df, args.db)
    print(f"Base de datos creada en '{args.db}'.")
    print(f"Estrellas almacenadas en la tabla 'estrellas': {len(df):,}")


if __name__ == "__main__":
    main()
