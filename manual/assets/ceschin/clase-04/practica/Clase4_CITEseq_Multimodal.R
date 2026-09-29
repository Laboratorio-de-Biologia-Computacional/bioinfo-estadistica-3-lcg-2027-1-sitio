# =============================================================================
# Bioinformática y Estadística 3 — Módulo scRNA-seq + CITE-seq
# Clase 4 — CITE-seq y análisis multimodal
# Docente: Dr. Danilo Ceschin
#
# CÓMO USAR ESTE SCRIPT EN RSTUDIO
#   1. Abre `Curso_scRNAseq.Rproj` (doble clic). Eso fija el directorio de
#      trabajo; sin eso las rutas relativas `data/` y `results/` no funcionan.
#   2. Ejecuta línea por línea o bloque por bloque con Ctrl+Enter
#      (Cmd+Enter en macOS). No hace falta correr todo el archivo de una sola vez.
#   3. El índice de secciones se abre con Ctrl+Shift+O (Cmd+Shift+O), o con el
#      botón del margen inferior izquierdo del editor.
#   4. La sección 0 (instalación) solo hace falta la primera vez en cada
#      computadora. Las siguientes veces corre solo las líneas de `library()`.
#   5. Los gráficos van al panel Plots. Usa el botón "Zoom" para verlos en
#      grande, o `guardar_fig()` para exportarlos a `figures/`.
# =============================================================================

# Bioinformática y Estadística 3 — Licenciatura en Ciencias Genómicas --------
#
## Módulo scRNA-seq + CITE-seq · Clase 4 de 4 --------------------------------
### CITE-seq: análisis conjunto de ARN y proteínas de superficie -------------
#
# **Docente:** Dr. Danilo Ceschin · **Fecha:** jueves 24 de septiembre de 2026 ·
# **Duración:** 2 h
#
# --------------------------------------------------------------------------
#
# La clase pasada terminamos con un problema abierto: los marcadores canónicos de
# inmunología son **proteínas de superficie**, pero nosotros medimos **ARNm**. `CD4`
# casi no se detecta en scRNA-seq aunque la proteína esté ahí. Anotamos linfocitos T CD4
# sin poder ver CD4.
#
# **CITE-seq mide las dos cosas en la misma célula.**
#
# Hoy vamos a ver por qué los datos de proteína no se pueden tratar como los de ARN
# —tienen una estructura de ruido completamente distinta— y cómo se combinan las dos
# modalidades en un único análisis. La estrategia que vamos a usar (WNN) es la misma que
# se aplica para combinar ARN y accesibilidad de cromatina en datos Multiome.
#
# --------------------------------------------------------------------------
#
### Objetivos de la clase ----------------------------------------------------
#
# 1. Entender la química de CITE-seq: ADTs, controles isotipo, hashing, diseño de panel.
# 2. Reconocer por qué el ruido de fondo domina la señal de ADT y qué implica para la
#    normalización.
# 3. Comparar **CLR** y **dsb** sobre los mismos anticuerpos y ver en qué cambia la
#    interpretación.
# 4. Integrar ARN + proteína con **WNN** e interpretar los pesos de modalidad por
#    célula.
# 5. Contrastar la anotación multimodal con la anotación por ARN de la clase 3.
#
# --------------------------------------------------------------------------

## 0. Preparación del entorno ------------------------------------------------
#
# Este bloque prepara el entorno local: crea las carpetas del proyecto,
# instala lo que falte e importa los paquetes. La instalación completa puede
# tardar la primera vez; después arranca en segundos.

# ---------------------------------------------------------------------------
# 0.1  Carpetas de trabajo y rutas
#
# Todo lo que el script descarga o produce vive dentro del proyecto:
#   data/     matrices descargadas de 10x Genomics (se bajan una sola vez)
#   results/  checkpoints .rds y tablas de resultados
#   figures/  figuras exportadas con guardar_fig()
#
# Estas carpetas viven en el disco, dentro del proyecto: NO se borran al cerrar
# RStudio. El checkpoint que guarda cada clase queda disponible para la siguiente.
# ---------------------------------------------------------------------------
dir_datos <- "data"
dir_res   <- "results"
dir_fig   <- "figures"
for (d in c(dir_datos, dir_res, dir_fig)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Rutas usadas a lo largo de las cuatro clases
f_h5_filt <- file.path(dir_datos, "pbmc10k_filtered.h5")   #  21 MB
f_h5_raw  <- file.path(dir_datos, "pbmc10k_raw.h5")        # 154 MB (solo clase 4)
ck1 <- file.path(dir_res, "clase1_pbmc_reducido.rds")
ck2 <- file.path(dir_res, "clase2_pbmc_clusterizado.rds")
ck3 <- file.path(dir_res, "clase3_pbmc_anotado.rds")
ck4 <- file.path(dir_res, "clase4_pbmc_multimodal.rds")

# Tamaño de archivo sin depender del shell del sistema: funciona igual en
# macOS, Linux y Windows.
tam_archivo <- function(f) {
  if (!file.exists(f)) {
    cat("[falta]", f, "\n")
  } else {
    cat(sprintf("%s - %.1f MB\n", f, file.size(f) / 1024^2))
  }
  invisible(file.exists(f))
}

# Exportar una figura al directorio figures/
guardar_fig <- function(p, nombre, w = 12, h = 6, dpi = 150) {
  ggplot2::ggsave(file.path(dir_fig, nombre), p, width = w, height = h, dpi = dpi)
  cat("figura guardada en", file.path(dir_fig, nombre), "\n")
}

# ---------------------------------------------------------------------------
# 0.2  Paquetes
#
# Instalación condicional: solo se instala lo que falta. La primera vez puede
# tardar bastante (Seurat y sus dependencias compilan código C++); las
# siguientes veces este bloque termina en segundos.
#
# REQUISITOS DEL SISTEMA (una sola vez, fuera de R):
#   - macOS  : instalar las Command Line Tools (`xcode-select --install`) y
#              HDF5 (`brew install hdf5`), que es lo que necesita `hdf5r`
#              para leer los archivos .h5 de 10x.
#   - Linux  : `sudo apt install libhdf5-dev libcurl4-openssl-dev libssl-dev
#              libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev`
#   - Windows: instalar Rtools de la misma versión mayor que tu R.
#
# Si `hdf5r` no compila, la alternativa es bajar el .tar.gz de la matriz desde
# 10x y usar Read10X() en lugar de Read10X_h5().
# ---------------------------------------------------------------------------
instalar_si_falta <- function(pkgs, fuente = c("CRAN", "Bioconductor")) {
  fuente <- match.arg(fuente)
  faltan <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(faltan)) return(invisible(NULL))
  message("Instalando desde ", fuente, ": ", paste(faltan, collapse = ", "))
  if (fuente == "CRAN") {
    install.packages(faltan)
  } else {
    BiocManager::install(faltan, ask = FALSE, update = FALSE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes",     quietly = TRUE)) install.packages("remotes")

instalar_si_falta(c("Seurat", "SeuratObject", "hdf5r", "tidyverse", "patchwork",
                    "cowplot", "ggridges", "pheatmap", "clustree", "R.utils",
                    "viridis", "gridExtra", "ggrepel"), "CRAN")

# Nota sobre `scrapper`, `viridis` y `gridExtra`: SingleR los declara en
# Suggests, no en Imports, así que NO se instalan junto con SingleR. Hacen falta
# igual en la clase 3: `scrapper` para SingleR(clusters = ...) y los otros dos
# para plotScoreHeatmap() y plotDeltaDistribution(). Sin ellos el error aparece
# recién en medio de la práctica, así que se instalan acá.
instalar_si_falta(c("SingleCellExperiment", "scuttle", "scDblFinder",
                    "SingleR", "celldex", "glmGamPoi", "scrapper"), "Bioconductor")

# presto acelera ~100x el test de Wilcoxon de FindAllMarkers. No está en CRAN.
if (!requireNamespace("presto", quietly = TRUE)) {
  remotes::install_github("immunogenomics/presto", upgrade = "never")
}

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(tidyverse)
  library(patchwork); library(cowplot)
})

cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("R:", R.version.string, "\n")

# El tamaño de las figuras lo define el panel Plots de RStudio.
# Para ver una figura en grande: botón "Zoom". Para exportarla: guardar_fig().

set.seed(42)  # reproducibilidad: TODO lo que sigue depende de esto

# Paquetes específicos de esta clase
if (!requireNamespace("dsb", quietly = TRUE)) install.packages("dsb", quiet = TRUE)
# `fast.km` existe desde dsb 2.0.0. Si hay una versión vieja instalada, la llamada de
# la sección 4.2 fallaría con "unused argument", así que la actualizamos acá.
if (packageVersion("dsb") < "2.0.0") install.packages("dsb", quiet = TRUE)
suppressPackageStartupMessages({ library(dsb); library(ggridges) })
cat("dsb:", as.character(packageVersion("dsb")), "\n")

# --------------------------------------------------------------------------
## 1. Qué es CITE-seq --------------------------------------------------------
#
# **CITE-seq** (*Cellular Indexing of Transcriptomes and Epitopes by Sequencing*) usa
# anticuerpos conjugados a un oligonucleótido de ADN con código de barras: los **ADT**
# (*antibody-derived tags*).
#
# El procedimiento:
#
# 1. Se tiñen las células con el panel de anticuerpos-ADT (como en citometría).
# 2. Se lavan (los anticuerpos no unidos, en teoría, se van).
# 3. Se encapsulan en gotas junto con las perlas de captura.
# 4. El código de barras de la célula se agrega **tanto al ARNm como al oligo del ADT**.
# 5. Se secuencian dos librerías separadas y se reúnen por código de barras celular.
#
# Resultado: para cada célula, un vector de expresión génica **y** un vector de
# abundancia de proteínas de superficie.
#
# | | ARN | ADT (proteína) |
# |---|---|---|
# | Nº de features | ~20.000–30.000 | 10–300 (el panel lo eliges tú) |
# | Conteos por célula | miles | decenas de miles (mucho más profundo) |
# | Sparsity | 90–95 % ceros | baja: **casi todo tiene señal** |
# | Ruido dominante | dropout (falsos ceros) | **fondo** (falsos positivos) |
#
# Esa última fila es la clave de toda la clase. **En ARN el problema es no ver lo que
# está. En proteína el problema es ver lo que no está.**
#
### 1.1. Diseño experimental -------------------------------------------------
#
# - **Panel.** De 10 a 300 anticuerpos. La elección es irreversible: lo que no pusiste,
#   no lo puedes medir después. Panel chico y bien elegido > panel grande y genérico.
# - **Controles isotipo.** Anticuerpos que no se unen a nada específico, del mismo
#   isotipo y conjugados al mismo tipo de oligo-ADT (en citometría sería el mismo
#   fluoróforo; acá no hay fluoróforo). Su señal **es** el fondo inespecífico. Deberían estar en todo
#   panel; en la práctica muchos datasets públicos no los tienen.
# - **Cell hashing.** Anticuerpos contra proteínas ubicuas (β2-microglobulina, MHC-I)
#   con un código por muestra. Permiten multiplexar varias muestras en un mismo carril y
#   —beneficio principal— **identificar doublets entre muestras** de forma directa.
#   Seurat los maneja con `HTODemux`.
# - **Limitación.** Solo proteínas de **superficie**. Nada intracelular, nada de
#   fosfoproteínas.
#
### 1.2. El dataset ----------------------------------------------------------
#
# El mismo PBMC 10k de las clases anteriores, ahora con su modalidad `Antibody Capture`:
# **17 anticuerpos TotalSeq-B**. Las células son literalmente las mismas.

url_filt <- paste0("https://cf.10xgenomics.com/samples/cell-exp/3.0.0/",
                   "pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")
if (!file.exists(f_h5_filt))
  download.file(url_filt, f_h5_filt, mode = "wb", quiet = TRUE)

data <- Read10X_h5(f_h5_filt)
str(data, max.level = 1)

adt_raw <- data[["Antibody Capture"]]
cat("\nNombres tal como vienen del archivo:\n")
print(head(rownames(adt_raw), 4))

# ---------------------------------------------------------------------------
# ⚠️ Los nombres de los anticuerpos traen el sufijo del kit ("CD3_TotalSeqB") y
# un guion bajo. Eso trae dos problemas:
#
#   1. Seurat NO admite guiones bajos en los nombres de features: al crear el
#      assay los convierte en guiones, sin avisar. Si después comparamos esa
#      matriz con otra leída directamente del .h5 —que conserva el guion bajo—,
#      los nombres ya no coinciden. dsb aborta con:
#
#          rows of cell and background matrices have mis-matching names
#
#   2. "CD3_TotalSeqB" es ilegible en un gráfico. Queremos ver "CD3".
#
# Lo resolvemos una sola vez, al leer: acortamos los nombres y hacemos nosotros
# la conversión que Seurat haría por su cuenta. Todas las matrices de ADT del
# script pasan por acá, así que quedan con nombres idénticos.
# ---------------------------------------------------------------------------
nombres_adt <- function(x) gsub("_", "-", sub("_TotalSeqB$", "", x))

rownames(adt_raw) <- nombres_adt(rownames(adt_raw))
cat("\nPanel de", nrow(adt_raw), "anticuerpos:\n")
print(rownames(adt_raw))

# Controles isotipo del panel. Son anticuerpos que NO reconocen nada: lo que
# miden es unión inespecífica pura, o sea el ruido técnico de cada célula.
# dsb los usa en el paso II para estimar ese ruido, y es su configuración
# recomendada.
isotipos <- grep("[Ii]so|[Cc]trl|IgG", rownames(adt_raw), value = TRUE)
cat("\nControles isotipo detectados:", paste(isotipos, collapse = ", "), "\n")

# --------------------------------------------------------------------------
## 2. La observación central: las gotas vacías tienen señal de proteína ------
#
# Este es el hecho que organiza todo el análisis de ADT.
#
# Vamos a bajar la matriz **raw** (todos los códigos de barras, incluidos los cientos de
# miles de gotas que **no contienen ninguna célula**) y mirar cuánta señal de ADT tienen
# esas gotas vacías.
#
# Si la tinción y el lavado fueran perfectos, deberían tener cero.

# ~154 MB, alrededor de un minuto
url_raw <- paste0("https://cf.10xgenomics.com/samples/cell-exp/3.0.0/",
                  "pbmc_10k_protein_v3/pbmc_10k_protein_v3_raw_feature_bc_matrix.h5")
if (!file.exists(f_h5_raw))
  download.file(url_raw, f_h5_raw, mode = "wb", quiet = TRUE)

raw <- Read10X_h5(f_h5_raw)
rna_raw_all <- raw[["Gene Expression"]]
adt_raw_all <- raw[["Antibody Capture"]]
rownames(adt_raw_all) <- nombres_adt(rownames(adt_raw_all))   # mismos nombres que arriba

cat("Códigos de barras en la matriz raw:", ncol(rna_raw_all), "\n")
cat("Códigos de barras con célula (matriz filtrada):", ncol(adt_raw), "\n")
cat("→ el", round(100 * (1 - ncol(adt_raw) / ncol(rna_raw_all)), 1), "% son gotas vacías\n")

# Métricas por código de barras
qc <- data.frame(
  barcode   = colnames(rna_raw_all),
  rna_size  = log10(Matrix::colSums(rna_raw_all) + 1),
  prot_size = log10(Matrix::colSums(adt_raw_all) + 1),
  n_gene    = Matrix::colSums(rna_raw_all > 0)
)
qc$es_celula <- qc$barcode %in% colnames(adt_raw)

ggplot(qc[qc$rna_size > 0 | qc$prot_size > 0, ],
       aes(rna_size, prot_size, color = es_celula)) +
  geom_point(size = 0.25, alpha = 0.4) +
  scale_color_manual(values = c("grey65", "firebrick"),
                     labels = c("gota vacía", "célula")) +
  labs(title = "Librería de ARN vs librería de proteína, por código de barras",
       subtitle = "Las gotas vacías (gris) tienen cientos de conteos de ADT y casi nada de ARN",
       x = "log10(conteos de ARN)", y = "log10(conteos de ADT)", color = NULL) +
  theme_bw() + theme(legend.position = "top")

# > **Esa nube gris a la izquierda es el hallazgo.** Son gotas sin célula, con señal de
#   ARN prácticamente nula (mediana: 8 conteos), y sin embargo con cientos de conteos
#   de ADT.
# >
# > ¿De dónde salen? De anticuerpos que quedaron libres en la suspensión después del
#   lavado. Se reparten por todas las gotas.
# >
# > **La consecuencia práctica:** cuando ves 90 conteos de CD15 en un linfocito T, no
#   sabes si son 90 de fondo o 85 de fondo y 5 de señal real (en este dataset, una gota
#   vacía tiene en mediana 88 conteos de CD15; lo vas a ver en la tabla de abajo). Y como las gotas
#   vacías miden **exactamente ese fondo**, tenemos un estimador directo de él. Esa es
#   la idea de dsb.

# Cuantificación del fondo: primero en total, después anticuerpo por anticuerpo
vacias <- qc[!qc$es_celula & qc$rna_size < 2 & qc$prot_size > 1, ]
cat("Conteos de ADT — mediana en células     :",
    round(median(10^qc$prot_size[qc$es_celula])), "\n")
cat("Conteos de ADT — mediana en gotas vacías:",
    round(median(10^vacias$prot_size)), "\n\n")

# El total engaña. Lo que importa es el fondo de CADA anticuerpo comparado con lo
# que mide ese anticuerpo en una célula típica.
fondo_ab <- data.frame(
  vacias  = apply(adt_raw_all[, vacias$barcode], 1, median),
  celulas = apply(adt_raw_all[, qc$barcode[qc$es_celula]], 1, median)
)
fondo_ab$pct_fondo <- round(100 * fondo_ab$vacias / fondo_ab$celulas, 1)
fondo_ab[order(-fondo_ab$pct_fondo), ]

# > **Cómo leer estos números.** En total, una gota vacía tiene unos 260 conteos de ADT
#   contra ~5.800 de una célula: el fondo parece chico, un 4–5 %. Pero no se reparte
#   parejo, y la tabla por anticuerpo lo muestra.
# >
# > - Para **CD3, CD4 o CD45RA**, antígenos que tiene la mayoría de las células de un
#   PBMC, el fondo es el 1–7 % de la mediana: hay mucha señal real encima.
# > - Para **CD8a, CD14, CD16, CD19 o CD15**, que la mayoría de las células NO tiene,
#   la célula típica mide apenas un poco más que una gota vacía: entre el 25 % y el
#   85 % de esos conteos es fondo.
# > - Los **controles isotipo** quedan en ~45–50 %: no se unen a nada, así que lo que
#   miden en una célula es fondo más la unión inespecífica de esa célula.
# >
# > Para un anticuerpo cuyo antígeno la célula no tiene, el fondo es el 100 % de lo que
#   medimos. Por eso, en proteína, el problema no es el cero: es el fondo.

# --------------------------------------------------------------------------
## 3. El objeto multimodal en Seurat v5 --------------------------------------
#
# Un `Seurat` puede contener varios **assays**. Vamos a construir uno con `RNA` y `ADT`
# sobre el mismo conjunto de células.

rna <- data[["Gene Expression"]]
adt <- data[["Antibody Capture"]]
rownames(adt) <- nombres_adt(rownames(adt))   # mismos nombres que arriba

# ---------------------------------------------------------------------------
# Partimos de las MISMAS células de la clase 3
#
# Si repitiéramos aquí el control de calidad y la búsqueda de doublets, el
# resultado no sería idéntico al de la clase 1: scDblFinder tiene un componente
# aleatorio y descarta unas decenas de células distintas en cada corrida. La
# comparación con la anotación de la clase 3 (§6) quedaría con células de más y
# de menos. Por eso, si el checkpoint existe, tomamos sus células y traemos su
# anotación. Si no existe, repetimos el QC de la clase 1.
# ---------------------------------------------------------------------------
if (file.exists(ck3)) {
  ant <- readRDS(ck3)
  celulas <- colnames(ant)
  cols_c3 <- intersect(c("clusters", "tipo_celular", "confianza_anotacion"),
                       colnames(ant@meta.data))
  meta_c3 <- ant@meta.data[, cols_c3, drop = FALSE]
  colnames(meta_c3) <- paste0(cols_c3, "_clase3")
  rm(ant); invisible(gc())

  pbmc <- CreateSeuratObject(rna[, celulas], project = "pbmc10k_cite",
                             min.cells = 3, min.features = 200)
  pbmc <- AddMetaData(pbmc, meta_c3)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  cat("Células tomadas de la clase 3:", ncol(pbmc), "\n")
} else {
  suppressPackageStartupMessages(library(scDblFinder))
  pbmc <- CreateSeuratObject(rna, project = "pbmc10k_cite", min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

  # QC por MAD, igual que en la clase 1
  io <- function(x, k = 5, hi = FALSE) { m <- median(x); d <- mad(x)
    if (hi) x > m + k * d else x < m - k * d | x > m + k * d }
  pbmc <- pbmc[, !(io(log1p(pbmc$nCount_RNA)) | io(log1p(pbmc$nFeature_RNA)) |
                   io(pbmc$percent.mt, 3, TRUE) | pbmc$percent.mt > 20)]

  sce  <- scDblFinder(as.SingleCellExperiment(pbmc), verbose = FALSE)
  pbmc <- pbmc[, colData(sce)$scDblFinder.class == "singlet"]
  cat("Sin checkpoint de la clase 3: QC repetido,", ncol(pbmc), "células\n")
}

# Agregamos la modalidad de proteína SOBRE LAS MISMAS CÉLULAS
pbmc[["ADT"]] <- CreateAssay5Object(counts = adt[, colnames(pbmc)])

Assays(pbmc)
pbmc

# Ahora hay dos conjuntos de métricas por célula
head(pbmc@meta.data[, c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT")])

VlnPlot(pbmc, features = c("nCount_RNA", "nCount_ADT", "nFeature_ADT"),
        pt.size = 0, ncol = 3, log = TRUE) & NoLegend()

# > `nFeature_ADT` es casi siempre 17 de 17: **cada célula tiene señal de todos los
#   anticuerpos del panel**. En ARN esto sería impensable. Otra vez: aquí el problema no
#   es el cero, es el fondo.

# --------------------------------------------------------------------------
## 4. Normalización de ADT: CLR vs dsb ---------------------------------------
#
### 4.1. CLR — centered log-ratio --------------------------------------------
#
# $$ \text{CLR}(x_i) = \log\left(\frac{x_i}{g(\mathbf{x})}\right), \qquad g(\mathbf{x})
# = \left(\prod_j x_j\right)^{1/n} $$
#
# Es la transformación estándar para **datos composicionales**: divide cada conteo por
# la media geométrica y toma logaritmo.
#
# `margin = 2` normaliza **por célula** (a través de los anticuerpos); `margin = 1`
# normaliza **por anticuerpo** (a través de las células).
#
# ⚠️ Ojo con el default: `NormalizeData()` usa `margin = 1`. Para ADT la viñeta oficial
# de CITE-seq de Seurat pasa explícitamente `margin = 2`, que es lo que hacemos acá.
# Si te olvidas del argumento, estás normalizando por otra cosa sin darte cuenta.
#
# **Qué hace bien:** corrige diferencias de profundidad de la librería de ADT entre
# células.
# **Qué no hace:** nada respecto del fondo. Si CD15 tiene 85 conteos de fondo en un
# linfocito T, después de CLR sigue teniendo señal — solo que reescalada.
#
### 4.2. dsb — denoised and scaled by background -----------------------------
#
# `dsb` ataca dos fuentes de ruido por separado:
#
# 1. **Ruido de fondo específico por anticuerpo** → lo estima **a partir de las gotas
#    vacías**, que son una medición directa del fondo. Resta esa media, anticuerpo por
#    anticuerpo.
# 2. **Ruido técnico específico por célula** → algunas células capturan más anticuerpo
#    libre que otras. dsb lo estima con un componente principal del "nivel de fondo por
#    célula" y lo elimina por regresión. Para construirlo usa los **controles isotipo** si el panel
#    los tiene —es el caso de este dataset, que trae tres— y, si no, un modelo de dos
#    poblaciones sobre la distribución de cada célula.
#
# El resultado está en unidades interpretables: **un valor cercano a 0 significa
# "indistinguible del fondo"**, y valores positivos son señal real por encima del fondo.
#
# > **Versión actual: dsb 2.0.1** (noviembre 2025). Trae el argumento `fast.km = TRUE`,
#   unas 10 veces más rápido, y la función `ModelNegativeADTnorm()` para el caso
#   —frecuente— en que no se dispone de la matriz raw con gotas vacías.

# --- Definir la población de gotas vacías que dsb usará como estimador del fondo
qc$es_celula <- qc$barcode %in% colnames(pbmc)

vacias_bc <- qc$barcode[
  !qc$es_celula &
  qc$rna_size  > 0.5 & qc$rna_size  < 2.5 &   # algo de ARN, pero poco: no es célula
  qc$prot_size > 1.5 & qc$prot_size < 3.5     # señal de ADT clara
]
cat("Gotas vacías seleccionadas para estimar el fondo:", length(vacias_bc), "\n")

# dsb necesita MATRICES (células x anticuerpos, densas)
cells_adt <- as.matrix(GetAssayData(pbmc, assay = "ADT", layer = "counts"))
empty_adt <- as.matrix(adt_raw_all[, vacias_bc])

dim(cells_adt); dim(empty_adt)

# Comprobación antes de llamar a dsb: las dos matrices tienen que tener los
# MISMOS anticuerpos, en el MISMO orden. Es el error más común de este paso.
stopifnot(identical(rownames(cells_adt), rownames(empty_adt)))

# Este panel SÍ tiene controles isotipo, así que usamos la configuración
# recomendada por dsb: denoise.counts = TRUE con use.isotype.control = TRUE.
adt_dsb <- DSBNormalizeProtein(
  cell_protein_matrix      = cells_adt,
  empty_drop_matrix        = empty_adt,
  denoise.counts           = TRUE,   # paso II: ruido técnico por célula
  use.isotype.control      = TRUE,   # recomendado cuando hay isotipos
  isotype.control.name.vec = isotipos,
  fast.km                  = TRUE    # dsb >= 2.0.0: ~10x más rápido
)

dim(adt_dsb)
round(adt_dsb[1:5, 1:4], 2)

# Guardamos las DOS normalizaciones para poder compararlas
# (a) CLR, en el assay ADT
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR",
                      margin = 2, verbose = FALSE)

# (b) dsb, en un assay aparte
pbmc[["ADTdsb"]] <- CreateAssay5Object(data = adt_dsb)

Assays(pbmc)

### 4.3. La comparación ------------------------------------------------------
#
# Miremos los mismos anticuerpos bajo las tres representaciones. Lo que buscamos es
# **bimodalidad**: una población negativa y una positiva bien separadas, como en un
# histograma de citometría.

abs_ver <- c("CD3", "CD4", "CD8a", "CD19", "CD14", "CD16")
abs_ver <- abs_ver[abs_ver %in% rownames(pbmc[["ADT"]])]
if (length(abs_ver) < 3) abs_ver <- rownames(pbmc[["ADT"]])[1:6]
abs_ver

hacer_ridge <- function(mat, titulo) {
  df <- as.data.frame(t(as.matrix(mat[abs_ver, , drop = FALSE])))
  df <- tidyr::pivot_longer(df, everything(), names_to = "ab", values_to = "v")
  ggplot(df, aes(v, ab, fill = ab)) +
    ggridges::geom_density_ridges(alpha = .75, scale = 1.4, linewidth = .3) +
    labs(title = titulo, x = NULL, y = NULL) +
    theme_bw() + theme(legend.position = "none", plot.title = element_text(size = 11))
}

crudo <- GetAssayData(pbmc, assay = "ADT",    layer = "counts")[abs_ver, ]
clr   <- GetAssayData(pbmc, assay = "ADT",    layer = "data")[abs_ver, ]
dsbm  <- GetAssayData(pbmc, assay = "ADTdsb", layer = "data")[abs_ver, ]

hacer_ridge(log1p(crudo), "Conteos crudos (log1p)") |
  hacer_ridge(clr, "CLR") |
  hacer_ridge(dsbm, "dsb")

### Qué mirar en estos tres paneles ------------------------------------------
#
# - **Crudo:** todo apilado contra el extremo izquierdo. La población positiva apenas se
#   despega.
# - **CLR:** mejor separación, pero la población negativa sigue siendo una masa ancha y
#   su posición **cambia de anticuerpo a anticuerpo**. No hay un umbral común.
# - **dsb:** la población negativa se centra alrededor de **0** para todos los
#   anticuerpos, porque 0 significa "igual al fondo". La población positiva se despega
#   claramente.
#
# Esa propiedad —que el cero signifique lo mismo en todos los anticuerpos— es la que
# permite comparar entre marcadores y poner umbrales interpretables, igual que en
# citometría.

# Biaxiales tipo citometría: CD4 vs CD8, bajo CLR y bajo dsb
if (all(c("CD4", "CD8a") %in% rownames(pbmc[["ADT"]]))) {
  df <- data.frame(
    CD4_clr = clr["CD4", ],  CD8_clr = clr["CD8a", ],
    CD4_dsb = dsbm["CD4", ], CD8_dsb = dsbm["CD8a", ]
  )
  p1 <- ggplot(df, aes(CD4_clr, CD8_clr)) +
    geom_point(size = .25, alpha = .3) + labs(title = "CLR") + theme_bw()
  p2 <- ggplot(df, aes(CD4_dsb, CD8_dsb)) +
    geom_point(size = .25, alpha = .3) +
    geom_hline(yintercept = 0, color = "red", linetype = 2) +
    geom_vline(xintercept = 0, color = "red", linetype = 2) +
    labs(title = "dsb (las líneas rojas marcan el nivel del fondo)") + theme_bw()

  print(p1 | p2)
}

# > Con dsb los cuadrantes tienen sentido: arriba-izquierda son CD8+CD4−, abajo-derecha
#   son CD4+CD8−, y el cuadrante doble-negativo está donde tiene que estar. Es
#   exactamente la lectura de un dot-plot de citometría.
# >
# > Fíjate en la nube intermedia de CD4 (alrededor de 10 en dsb, entre el doble negativo
#   y los T CD4). No es ruido: son los monocitos, que expresan CD4 en superficie en
#   menor nivel que los linfocitos T. Con CLR esa población queda mezclada con los
#   extremos; con dsb aparece como un nivel propio.
#
### 4.4. Otras opciones que conviene conocer ---------------------------------
#
# | Método | Cuándo |
# |---|---|
# | **dsb** | Tienes la matriz raw con gotas vacías. Es la opción por defecto. |
# | **`dsb::ModelNegativeADTnorm()`** | No tienes la matriz raw (te dieron solo la filtrada). Modela la población negativa directamente. |
# | **ADTnorm** (*Nat Commun* 2025) | Vas a **integrar varios lotes o estudios**. Alinea los picos de las poblaciones negativa y positiva entre datasets, que es donde CLR y dsb fallan. |
# | **totalVI** (Python, scvi-tools) | Modelo generativo que trata ARN y proteína conjuntamente, con el fondo como parte del modelo. Potente y más costoso. |

# Alternativa sin gotas vacías (útil para datasets públicos que solo publican
# la matriz filtrada). No la usamos hoy, queda como referencia.
# adt_dsb2 <- ModelNegativeADTnorm(cell_protein_matrix = cells_adt,
#                                  denoise.counts = TRUE,
#                                  use.isotype.control = TRUE,
#                                  isotype.control.name.vec = isotipos,
#                                  fast.km = TRUE)
cat("Ver comentario.\n")

# --------------------------------------------------------------------------
## 5. Integración multimodal: WNN --------------------------------------------
#
# Tenemos dos vistas de las mismas células. ¿Cómo las combinamos?
#
# | Estrategia | Problema |
# |---|---|
# | Concatenar las matrices | 20.000 genes contra 17 proteínas: la proteína queda enterrada. |
# | Analizar cada una por separado | Dos UMAPs y dos clusterings que no se hablan. |
# | **Ponderar por célula (WNN)** | ✅ Lo que vamos a hacer. |
#
# **Weighted Nearest Neighbors** aprende, **para cada célula individualmente**, cuánto
# pesa cada modalidad. La intuición es sencilla: evalúa qué tan bien predice cada
# modalidad el perfil de esa célula comparado con sus vecinos, y le da más peso a la que
# resulta más informativa **para esa célula**.
#
# Esto refleja algo real: para separar un monocito de un linfocito T, el ARN alcanza y
# sobra. Para separar T CD4 naive de T CD4 de memoria, la proteína (CD45RA, CD45RO,
# CD127) es mucho más informativa.
#
# > Es la **misma función** que se usa para combinar ARN y accesibilidad de cromatina en
#   datos Multiome. Cambian las modalidades, no el método.

# --- Rama ARN
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc, verbose = FALSE) |>
        FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
        ScaleData(verbose = FALSE) |>
        RunPCA(npcs = 30, reduction.name = "pca", verbose = FALSE)

# --- Rama proteína (usamos dsb)
DefaultAssay(pbmc) <- "ADTdsb"
VariableFeatures(pbmc) <- rownames(pbmc[["ADTdsb"]])   # con 17 anticuerpos, todos
pbmc <- ScaleData(pbmc, verbose = FALSE) |>
        RunPCA(npcs = 15, reduction.name = "apca",
               reduction.key = "aPC_", verbose = FALSE)

Reductions(pbmc)

pbmc <- FindMultiModalNeighbors(
  pbmc,
  reduction.list = list("pca", "apca"),
  dims.list      = list(1:30,  1:15),
  modality.weight.name = "RNA.weight",
  verbose = FALSE
)

pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn",
                reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
                verbose = FALSE)

# Leiden (igraph), el mismo algoritmo que elegimos en la clase 2: garantiza
# comunidades conectadas, y en Seurat >= 5.3.1 no necesita Python.
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 4,
                     leiden_method = "igraph", resolution = 0.8, verbose = FALSE)
pbmc$wnn_clusters <- pbmc$seurat_clusters
cat("Clusters WNN:", nlevels(pbmc$wnn_clusters), "\n")

# UMAPs de cada modalidad por separado, para comparar
pbmc <- RunUMAP(pbmc, reduction = "pca",  dims = 1:30,
                reduction.name = "rna.umap", reduction.key = "rnaUMAP_", verbose = FALSE)
pbmc <- RunUMAP(pbmc, reduction = "apca", dims = 1:15,
                reduction.name = "adt.umap", reduction.key = "adtUMAP_", verbose = FALSE)

(DimPlot(pbmc, reduction = "rna.umap", group.by = "wnn_clusters", label = TRUE, label.size = 3) +
   ggtitle("UMAP solo ARN") + NoLegend()) |
(DimPlot(pbmc, reduction = "adt.umap", group.by = "wnn_clusters", label = TRUE, label.size = 3) +
   ggtitle("UMAP solo proteína") + NoLegend()) |
(DimPlot(pbmc, reduction = "wnn.umap", group.by = "wnn_clusters", label = TRUE, label.size = 3) +
   ggtitle("UMAP WNN (ambas)") + NoLegend())

### 5.1. Los pesos de modalidad son un resultado biológico -------------------
#
# `RNA.weight` es un número por célula, entre 0 y 1: cuánto pesó el ARN en la
# construcción del grafo para esa célula. `1 - RNA.weight` es el peso de la proteína.
#
# No es un parámetro técnico: **es interpretable**.

p1 <- VlnPlot(pbmc, "RNA.weight", group.by = "wnn_clusters", pt.size = 0, sort = TRUE) +
  NoLegend() + labs(title = "Peso de la modalidad ARN por cluster",
                    subtitle = "Valores bajos = la PROTEÍNA fue más informativa para esas células")
p2 <- FeaturePlot(pbmc, "RNA.weight", reduction = "wnn.umap") +
  ggtitle("Peso del ARN sobre el UMAP WNN")
p1 | p2

# > **Cómo se lee.** El violín está ordenado de mayor a menor. Los clusters de la
#   izquierda (peso de ARN alto) son aquellos donde
#   el ARN ya alcanzaba: en nuestros datos, las DC convencionales y plasmacitoides, los
#   monocitos CD16+ y los B de memoria, linajes con perfiles transcripcionales muy
#   distintivos. Hacia la derecha quedan las subpoblaciones de linfocitos T, que
#   transcripcionalmente se parecen mucho entre sí y se separan bien por superficie
#   (CD45RA/CD45RO, CD127, CD25).
# >
# > **Pero fíjate en los dos clusters con el peso de ARN más bajo de todos.** No son
#   poblaciones T: son un cluster de monocitos con señal alta de un control isotipo y
#   un cluster mixto que expresa marcadores de varios linajes a la vez. Los vamos a
#   identificar en §6. La lección: un peso de proteína alto quiere decir "la proteína
#   separó a estas células", no "aquí hay biología nueva". Un artefacto de la proteína
#   también separa células.

# Los marcadores de proteína sobre el UMAP WNN
DefaultAssay(pbmc) <- "ADTdsb"
mostrar <- intersect(c("CD3", "CD4", "CD8a", "CD19", "CD14", "CD16", "CD56", "CD45RA"),
                     rownames(pbmc[["ADTdsb"]]))

FeaturePlot(pbmc, features = mostrar, reduction = "wnn.umap",
            ncol = 4, cols = c("lightgrey", "darkgreen"), min.cutoff = 0) &
  theme(plot.title = element_text(size = 11)) & NoLegend()

# El mismo marcador en las dos modalidades, lado a lado.
# Este es el argumento de la sección 1 de la clase 3, hecho visible.
DefaultAssay(pbmc) <- "RNA"
p_rna <- FeaturePlot(pbmc, "CD4", reduction = "wnn.umap") +
  ggtitle("CD4 — ARNm") + theme(plot.title = element_text(size = 12))
DefaultAssay(pbmc) <- "ADTdsb"
p_adt <- FeaturePlot(pbmc, "CD4", reduction = "wnn.umap",
                     cols = c("lightgrey", "darkgreen"), min.cutoff = 0) +
  ggtitle("CD4 — proteína (dsb)") + theme(plot.title = element_text(size = 12))

p_rna | p_adt

# > **Este par de gráficos es el argumento entero de la clase.** El ARNm de `CD4` se ve
#   sobre todo en los monocitos (abajo a la izquierda), que lo expresan; en los propios
#   linfocitos T CD4 aparece en unos pocos puntos dispersos, indistinguibles del ruido.
#   La proteína CD4 marca de forma nítida toda la población T CD4, y más tenue a los
#   monocitos: el mismo patrón biológico, medido sin dropout.
# >
# > En la clase 3 anotamos "T CD4" sin poder ver CD4. Funcionó porque usamos `IL7R` como
#   sustituto. Aquí no hace falta el sustituto.

# --------------------------------------------------------------------------
## 6. Anotación asistida por proteína ----------------------------------------

DefaultAssay(pbmc) <- "ADTdsb"
# scale = FALSE: el color es el valor dsb medio, sin reescalar. Con el escalado por
# defecto (z-score por anticuerpo) perderíamos justamente lo que dsb nos da: que el
# 0 signifique "igual al fondo" en todos los anticuerpos.
DotPlot(pbmc, features = rownames(pbmc[["ADTdsb"]]), group.by = "wnn_clusters",
        scale = FALSE) +
  RotatedAxis() +
  labs(title = "Perfil de superficie de cada cluster WNN",
       subtitle = "Color = valor dsb medio (0 = fondo) · Tamaño = % de células con dsb > 0")

# Marcadores de proteína por cluster (17 features: es instantáneo)
Idents(pbmc) <- "wnn_clusters"
mk_adt <- FindAllMarkers(pbmc, assay = "ADTdsb", only.pos = TRUE,
                         logfc.threshold = 0.3, verbose = FALSE)

mk_adt |> group_by(cluster) |> slice_max(avg_log2FC, n = 3) |>
  select(cluster, gene, avg_log2FC) |> as.data.frame()

# > **Tres clusters que no son un tipo celular nuevo.** Antes de etiquetar, busca en la
#   tabla y en el DotPlot los clusters cuyo marcador principal es raro:
# >
# > - **Un cluster de monocitos CD14+ cuyo marcador principal es `IgG2b-control`.** Un
#   control isotipo no reconoce nada: esas células pegan anticuerpo de forma
#   inespecífica (los monocitos tienen receptores Fc). En ARN son monocitos idénticos
#   a los demás; la proteína los separó por un artefacto.
# > - **Un cluster con CD3, CD4, CD14, CD15 y CD16 a la vez**, y en ARN `PPBP` (plaqueta)
#   y `LYZ` (monocito). Ninguna célula tiene todo eso: son agregados y doublets
#   heterotípicos que scDblFinder no atrapó. Sus células venían de casi todos los tipos
#   de la clase 3, repartidas.
# > - **Un cluster chico donde se juntan las DC plasmacitoides y los plasmablastos.** El
#   ARN los separaba sin problema; en el panel no hay ningún anticuerpo que los
#   distinga (faltan CD123, CD38, CD27), y WNN los fusionó. Solo ves lo que pusiste en
#   el panel.

# ---------------------------------------------------------------------------
# Edita ESTE VECTOR, igual que en la clase 3, pero ahora con evidencia de
# ARN **y** de proteína.
#
# Los nombres son los identificadores de cluster tal como los devolvió Leiden
# (empiezan en 1). Si tu corrida dio otro número de clusters, ajusta el vector a
# TU tabla: el código de abajo avisa si sobra o falta alguno.
# ---------------------------------------------------------------------------
DefaultAssay(pbmc) <- "RNA"
etiquetas_wnn <- c(
  "1"  = "Monocito CD14+",                # CD14; LYZ
  "2"  = "T CD4 naive",                   # CD3, CD4, CD45RA
  "3"  = "NK",                            # CD16, CD45RA; GNLY, GZMB
  "4"  = "T CD4 memoria",                 # CD127, CD45RO
  "5"  = "T regulador / activado",        # CD25 alto, PD-1, TIGIT; FOXP3 en ARN
  "6"  = "MAIT",                          # CD8a, CD3, CD127; SLC4A10
  "7"  = "T CD8 memoria",                 # CD8a, TIGIT, CD45RO; GZMK
  "8"  = "T CD8 naive",                   # CD8a, CD45RA
  "9"  = "B naive",                       # CD19, CD45RA
  "10" = "Monocito CD14+ (isotipo alto)", # IgG2b-control: unión inespecífica
  "11" = "B memoria",                     # CD19 alto; IGHG1 en ARN
  "12" = "Agregados / doublets",          # CD4 + CD15 + CD16 a la vez; PPBP
  "13" = "NK CD56bright",                 # CD56 alto, CD16 bajo
  "14" = "Monocito CD16+",                # CD16 y algo de CD14: intermedio
  "15" = "DC convencional",               # sin ADT propio; CD1C en ARN
  "16" = "T CD4 (naive/memoria)",         # cluster chico de transición
  "17" = "T no convencionales",           # sin ADT propio; verificar con ARN
  "18" = "DC plasmacitoide",              # sin ADT propio; CLEC4C en ARN
  "19" = "Plasmablasto"                   # sin ADT propio; MZB1, IGHA1 en ARN
)

sobran <- setdiff(names(etiquetas_wnn), levels(pbmc$wnn_clusters))
if (length(sobran)) {
  cat("Etiquetas para clusters que NO existen:", paste(sobran, collapse = ", "),
      "— revisa el vector\n")
  etiquetas_wnn <- etiquetas_wnn[setdiff(names(etiquetas_wnn), sobran)]
}
faltan <- setdiff(levels(pbmc$wnn_clusters), names(etiquetas_wnn))
if (length(faltan)) etiquetas_wnn[faltan] <- "Sin asignar"

# Nombramos por célula antes de asignar: al indexar `etiquetas_wnn` con los
# clusters, el vector resultante queda nombrado por id de cluster, y Seurat
# rechaza la metadata con "No cell overlap between new meta data and Seurat
# object". `factor()` tampoco borra esos nombres.
pbmc$tipo_wnn <- setNames(factor(etiquetas_wnn[as.character(pbmc$wnn_clusters)],
                                 levels = unique(etiquetas_wnn)),
                          colnames(pbmc))

DimPlot(pbmc, reduction = "wnn.umap", group.by = "tipo_wnn",
        label = TRUE, repel = TRUE, label.size = 3.5) +
  ggtitle("Anotación multimodal (WNN)") + NoLegend()

### 6.1. Contraste con la anotación de la clase 3 ----------------------------
#
# Las células son las mismas de la clase 3 y su etiqueta ya está en la metadata
# (`tipo_celular_clase3`), así que la comparación es directa. Lo interesante no es la
# concordancia global, sino **dónde** difieren.

if ("tipo_celular_clase3" %in% colnames(pbmc@meta.data)) {
  tab <- table(RNA_clase3 = pbmc$tipo_celular_clase3, WNN_clase4 = pbmc$tipo_wnn)
  suppressPackageStartupMessages(library(pheatmap))
  pheatmap(log1p(tab), cluster_rows = FALSE, cluster_cols = FALSE,
           display_numbers = tab, number_format = "%.0f", fontsize_number = 8,
           angle_col = 90,
           main = "Anotación solo ARN (clase 3, filas) vs multimodal (clase 4, columnas)")
} else {
  cat("Falta el checkpoint de la clase 3 (", ck3, "). Corre la clase 3 primero.\n")
}

# > **Qué mirar en la tabla.**
# >
# > - **"T CD4 memoria" de la clase 3 se abre en tres:** memoria convencional, T
#   reguladores (CD25 alto, CD127 bajo) y una población PD-1+ TIGIT+. El ARN los tenía
#   juntos; la proteína los separa. Es el caso para el que WNN existe.
# > - **Los monocitos CD14+ se parten en dos**, pero por un motivo técnico: el cluster
#   de isotipo alto. Aquí el desacuerdo lo resuelve el ARN, no la proteína.
# > - **Unas 160 células de todos los tipos terminan en "Agregados / doublets".** La
#   proteína detectó lo que scDblFinder no vio.

### 6.2. ¿Qué dice la proteína de la anotación dudosa de la clase 3? ---------
#
# En la clase 3 el cluster "T CD8 naive" terminó con confianza **baja**: Blueprint por
# cluster decía CD4, pero el 79 % de sus células, Monaco y los marcadores decían CD8.
# Hoy tenemos la medición que faltaba: CD4 y CD8a como proteína, en esas mismas células.

if ("tipo_celular_clase3" %in% colnames(pbmc@meta.data)) {
  tipos_t <- c("T CD4 naive", "T CD4 memoria", "T CD8 naive", "T CD8 memoria", "MAIT")
  df_t <- data.frame(
    CD4  = dsbm["CD4", ],
    CD8a = dsbm["CD8a", ],
    tipo = pbmc$tipo_celular_clase3,
    conf = if ("confianza_anotacion_clase3" %in% colnames(pbmc@meta.data))
             pbmc$confianza_anotacion_clase3 else NA
  )
  df_t <- df_t[df_t$tipo %in% tipos_t, ]
  df_t$tipo <- factor(df_t$tipo, levels = tipos_t)

  # Una célula es "CD8 por proteína" si su señal de CD8a supera a la de CD4
  resumen_t <- df_t |>
    group_by(tipo) |>
    summarise(n = n(),
              conf_clase3  = first(as.character(conf)),
              pct_CD8_prot = round(100 * mean(CD8a > CD4), 1))
  print(as.data.frame(resumen_t), row.names = FALSE)

  ggplot(df_t, aes(CD4, CD8a)) +
    geom_point(size = 0.3, alpha = 0.4, colour = "#0B2E33") +
    geom_hline(yintercept = 0, colour = "firebrick", linetype = 2) +
    geom_vline(xintercept = 0, colour = "firebrick", linetype = 2) +
    facet_wrap(~ tipo, nrow = 1) +
    labs(title = "Proteína CD4 vs CD8a (dsb), según la etiqueta por ARN de la clase 3",
         subtitle = "La clase 3 declaró confianza BAJA para \"T CD8 naive\": la proteína decide",
         x = "CD4 (dsb)", y = "CD8a (dsb)") +
    theme_bw() + theme(strip.text = element_text(size = 11))
}

# > **La proteína confirma la etiqueta.** El 96 % de las células que la clase 3 llamó
#   "T CD8 naive" es CD8a+ CD4− por proteína. Blueprint por cluster estaba equivocado;
#   la mayoría por célula, Monaco y los marcadores tenían razón. "Confianza baja" no
#   significaba "incorrecta": significaba "hace falta otra evidencia". Hoy la tenemos.
# >
# > Mira también "T CD8 memoria": ~87 % CD8a+, y un grupo doble negativo (CD4− CD8−)
#   abajo a la izquierda. Son células T que no son CD4 ni CD8 —probablemente γδ—, lo
#   que Monaco ya sugería en la clase 3 ("Non-Vd2 gd T cells") y dejó ese cluster con
#   confianza media.

saveRDS(pbmc, ck4)
tam_archivo(ck4)
sessionInfo()

# --------------------------------------------------------------------------
## 7. Síntesis del módulo ----------------------------------------------------
#
# **De esta clase:**
#
# 1. En ADT el ruido dominante es el **fondo**, no el dropout. Es la inversión exacta
#    del problema del ARN.
# 2. **Las gotas vacías miden el fondo directamente.** Esa observación es la base de
#    dsb.
# 3. **CLR** corrige profundidad; **dsb** corrige fondo y devuelve valores donde el 0
#    significa algo. **ADTnorm** resuelve la integración entre lotes; **totalVI** modela
#    todo conjuntamente.
# 4. **WNN** pondera las modalidades por célula, y esos pesos son un resultado biológico
#    interpretable.
# 5. La proteína recupera marcadores —CD4 es el caso de manual— que el ARNm simplemente
#    no ve, y resuelve las anotaciones dudosas del ARN (T CD8 naive, T reguladores).
# 6. **La proteína también tiene artefactos propios**: unión inespecífica (isotipos
#    altos) y agregados. Un peso de proteína alto no garantiza biología.
#
# **De las cuatro clases juntas:**
#
# | Clase | Idea que hay que llevarse |
# |---|---|
# | 1 | Normalizar y reducir dimensiones **determina** la estructura que después vas a interpretar. UMAP es un mapa, no una medición. |
# | 2 | El clustering siempre devuelve clusters. Los p-valores post-clustering no son p-valores. Para condiciones: pseudobulk. |
# | 3 | La anotación es una hipótesis. Se sostiene con evidencia convergente y se reporta con su nivel de confianza. |
# | 4 | Cada modalidad tiene su propia estructura de ruido. Integrar no es concatenar. |
#
# En el análisis de cromatina y en los flujos multimodales van a reencontrar
# `FindTransferAnchors` (clase 3), `FindMultiModalNeighbors` (esta clase) y toda la
# lógica de que **cada modalidad exige su propio preprocesamiento antes de poder
# combinarla**.
