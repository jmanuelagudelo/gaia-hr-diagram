#Diagrama HR de Gaia DR3
set -euo pipefail #mostrar errores 

# --- Rutas ----
DATOS="datos"
RESULTADOS="resultados"
CSV="${DATOS}/datos_gaia.csv"
DB="${DATOS}/datos_mision.db"
PNG="${RESULTADOS}/diagrama_hr.png"
ENDPOINT="https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

#muestra local d < ~50 pc señal/ruido de paralaje alto (Plx/e_Plx > 10).
QUERY='SELECT TOP 50000 Plx, e_Plx, Gmag, BPmag, RPmag FROM "I/355/gaiadr3" WHERE Plx > 20 AND Plx / e_Plx > 10 AND Gmag IS NOT NULL'


# --- Utilidades -----
#Funciones para mensajes de información y errores

info() { echo "[INFO] $*"; }
error() { echo "[ERROR] $*" >&2; }

info "Comprobando dependencias..."

if ! python3 -c "import pandas, numpy, matplotlib" >/dev/null 2>&1; then
    error "Faltan módulos de Python (pandas, numpy, matplotlib)."
    error "Instálalos con:  python3 -m pip install -r requirements.txt"
    exit 1
fi

mkdir -p "${DATOS}" "${RESULTADOS}"

# --- 2. Descargar datos (con indicador giratorio) ----
info "Descargando datos de Gaia DR3 desde VizieR/TAP..."

#Puesto que VizieR recibe la URL con + por cada espacio, usamos quote_plus
#para codificar la consulta. Además llamamos por POST.
QUERY_ENCODED=$(python3 -c 'import urllib.parse, sys; print(urllib.parse.quote_plus(sys.argv[1]))' "${QUERY}")

# Lanzamos wget en segundo plano y, mientras tanto, mostramos un indicador
wget -q -O "${CSV}" \
    --timeout=60 --tries=3 --retry-connrefused \
    --post-data="REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=${QUERY_ENCODED}" \
    "${ENDPOINT}" &
PID_WGET=$!

spinner='|/\-'
i=0
while kill -0 "$PID_WGET" 2>/dev/null; do
    char="${spinner:$((i % 4)):1}"
    printf "\r  Descargando... %s " "$char"
    i=$((i + 1))
    sleep 0.1
done

if wait "$PID_WGET"; then
    printf "\r  Descarga completada.        \n"
else
    printf "\n"
    error "La descarga falló. VizieR puede limitar temporalmente las peticiones; reintenta en unos minutos."
    exit 1
fi
info "Datos guardados en '${CSV}'."

# --- 3. Construir base de datos ---
python3 constructor_db.py --csv "${CSV}" --db "${DB}"

# --- 4. Visualización ---
info "Generando el diagrama Hertzsprung-Russell..."
python3 analisis_visual.py --db "${DB}" --out "${PNG}"

info "Pipeline completado. Resultado: ${PNG}"
if [[ -f "${PNG}" ]]; then
    if [[ "$OSTYPE" == "darwin"* ]]; then
        open "${PNG}"
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        xdg-open "${PNG}"
    else
        error "No sé cómo abrir automáticamente la imagen en este sistema."
    fi
else
    error "No se encontró la imagen: ${PNG}"
    exit 1
fi