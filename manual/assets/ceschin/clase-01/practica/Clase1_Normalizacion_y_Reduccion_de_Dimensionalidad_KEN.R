# =============================================================================
# Bioinformática y Estadística 3 — Módulo scRNA-seq + CITE-seq
# Clase 1 — Normalización, escalamiento y reducción de dimensionalidad
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
#
# SI TRABAJAS EN EL CLÚSTER `ken` (variante de esta copia)
#   El entorno ya está instalado: aquí no se instala nada, se verifica.
#   La sesión tiene que arrancar con estos módulos, y en este orden:
#
#       module load shared hdf5/1.14.0 gcc/14.2.0 r/4.6.0
#
#   `hdf5` no es opcional ni es sólo para instalar: sin él, `hdf5r` figura como
#   instalado y aun así falla al cargar, con `libhdf5_hl.so.310: cannot open
#   shared object file`. Y sin `hdf5r` no hay `Read10X_h5()`, que es la primera
#   instrucción de la sección 1.
#
#   El proyecto va en `/mnt/data/bioinfo3/<tu_usuario>/`, NO en tu `$HOME`: los
#   cuatro objetos `.rds` de las cuatro clases pesan cientos de MB cada uno y la
#   cuota del `$HOME` es de 5 GB.
# =============================================================================

# Bioinformática y Estadística 3 — Licenciatura en Ciencias Genómicas --------
#
## Módulo scRNA-seq + CITE-seq · Clase 1 de 4 --------------------------------
### Normalización, escalamiento y reducción de dimensionalidad ---------------
#
# **Docente:** Dr. Danilo Ceschin · **Fecha:** martes 8 de septiembre de 2026 ·
# **Duración:** 2 h
#
# --------------------------------------------------------------------------
#
# Este módulo es autosuficiente: trae su propio conjunto de datos y hace su propio
# control de calidad. El script descarga la matriz de PBMC 10k, la filtra en la sección 2
# y sigue desde ahí. Al final de ese paso tenemos una matriz de conteos limpia... y
# todavía no podemos comparar nada.
#
# El problema es simple de enunciar: **dos células con el mismo estado biológico pueden
# tener conteos totales muy distintos**, porque la captura de ARNm y la profundidad de
# secuenciación varían de gota en gota. Si comparamos conteos crudos, lo que vamos a
# encontrar es cuánto se secuenció cada célula, no qué le pasa a cada célula.
#
# Todo lo que hagamos hoy —normalizar, elegir genes, escalar, reducir dimensiones—
# **determina la estructura que vamos a interpretar en las próximas tres clases**. No
# son pasos de plomería: son decisiones analíticas.
#
# --------------------------------------------------------------------------
#
### Objetivos de la clase ----------------------------------------------------
#
# 1. Entender qué corrige y qué **no** corrige cada método de normalización.
# 2. Elegir genes variables (HVGs) con criterio y ver el efecto de esa elección.
# 3. Ejecutar PCA y justificar cuántas componentes conservar.
# 4. Ejecutar t-SNE y UMAP, y —sobre todo— saber **qué no se puede leer** en esas
#    figuras.
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
# ¿Estamos en el clúster o en una computadora personal? Lmod define LMOD_CMD, y
# el nombre del nodo empieza por `ken`. Con cualquiera de las dos basta.
EN_CLUSTER <- nzchar(Sys.getenv("LMOD_CMD")) ||
  grepl("^ken", Sys.info()[["nodename"]], ignore.case = TRUE)

if (EN_CLUSTER) {
  message("Entorno detectado: clúster. No se instala nada; sólo se verifica.")
}

instalar_si_falta <- function(pkgs, fuente = c("CRAN", "Bioconductor")) {
  fuente <- match.arg(fuente)
  faltan <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(faltan)) return(invisible(NULL))

  # En el clúster la biblioteca es compartida y de sólo lectura para el grupo.
  # Intentar instalar ahí no repara nada: abre el diálogo de biblioteca personal
  # y, con todo el grupo a la vez, se lleva por delante la primera media hora de
  # clase. Preferimos detenernos con un mensaje que diga exactamente qué pedir.
  if (EN_CLUSTER) {
    stop("Faltan paquetes en el entorno del clúster: ",
         paste(faltan, collapse = ", "),
         "\n  Esto NO se resuelve desde aquí: la biblioteca del módulo es",
         " compartida.\n  Avisa a quien coordina el curso y sigue con la clase;",
         " no intentes instalarlos.",
         call. = FALSE)
  }

  message("Instalando desde ", fuente, ": ", paste(faltan, collapse = ", "))
  if (fuente == "CRAN") {
    install.packages(faltan)
  } else {
    BiocManager::install(faltan, ask = FALSE, update = FALSE)
  }
}

if (!EN_CLUSTER) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("remotes",     quietly = TRUE)) install.packages("remotes")
}

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
if (!EN_CLUSTER && !requireNamespace("presto", quietly = TRUE)) {
  remotes::install_github("immunogenomics/presto", upgrade = "never")
}

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject)
  # `tidyverse` sólo adjunta a otros ocho paquetes; no aporta código propio. En
  # el clúster están los ocho —ggplot2 4.0.3, dplyr 1.2.1, las versiones que fija
  # envs/renv.lock— y falta únicamente el meta-paquete. Este guion usa dos.
  library(ggplot2); library(dplyr)
  library(patchwork); library(cowplot)
})

cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("R:", R.version.string, "\n")

# El tamaño de las figuras lo define el panel Plots de RStudio.
# Para ver una figura en grande: botón "Zoom". Para exportarla: guardar_fig().

set.seed(42)  # reproducibilidad: TODO lo que sigue depende de esto

# --------------------------------------------------------------------------
## 1. Los datos --------------------------------------------------------------
#
# Vamos a usar el mismo conjunto de datos durante **las cuatro clases**: PBMC (células
# mononucleares de sangre periférica) de un donante sano, procesadas por 10x Genomics
# con el kit *Cell Surface Protein*.
#
# Esto significa que el archivo contiene **dos modalidades medidas en las mismas
# células**:
#
# | Modalidad | Qué es | Cuándo la usamos |
# |---|---|---|
# | `Gene Expression` | ~33.000 genes, conteos de ARNm | clases 1, 2 y 3 |
# | `Antibody Capture` | 17 anticuerpos con oligo (ADT) | clase 4 |
#
# Hoy usamos solo la primera. En la clase 4 vamos a agregar la proteína **sobre estas
# mismas células**, y vamos a poder preguntarnos si la proteína confirma o corrige lo
# que dijo el ARN.
#
# > Fuente: [10x Genomics — 10k PBMCs with
#   TotalSeq-B](https://www.10xgenomics.com/datasets), matriz filtrada por Cell Ranger
#   (21 MB).

# Descarga de la matriz filtrada (~21 MB, unos 15 segundos)
url_filt <- paste0("https://cf.10xgenomics.com/samples/cell-exp/3.0.0/",
                   "pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")

if (!file.exists(f_h5_filt)) {
  download.file(url_filt, f_h5_filt, mode = "wb", quiet = TRUE)
}
tam_archivo(f_h5_filt)

# Read10X_h5 devuelve una LISTA cuando el archivo tiene más de un tipo de feature
data <- Read10X_h5(f_h5_filt)
str(data, max.level = 1)

### 1.1. La matriz dispersa --------------------------------------------------
#
# Fíjate en la clase del objeto: `dgCMatrix`. Es una **matriz dispersa** (sparse): en
# vez de guardar todos los ceros, guarda solo las posiciones y valores de las entradas
# no nulas.
#
# En scRNA-seq entre el 90 % y el 95 % de las entradas son cero (la mayoría de los genes
# no se detecta en la mayoría de las células). Guardar eso como matriz densa sería
# absurdo.

rna <- data[["Gene Expression"]]

cat("Clase:", class(rna), "\n")
cat("Dimensiones:", nrow(rna), "genes x", ncol(rna), "células\n")
cat("Proporción de ceros:", round(100 * (1 - length(rna@x) / (nrow(rna) * ncol(rna))), 2), "%\n\n")

# Cuánto ocupa cada representación
cat("Como matriz dispersa:", format(object.size(rna), units = "MB"), "\n")
cat("Como matriz densa   :",
    format(object.size(numeric(1)) * nrow(rna) * ncol(rna), units = "MB"), "(estimado)\n")

# Un vistazo a las esquinas de la matriz
rna[c("CD3E", "CD14", "MS4A1", "NKG7"), 1:6]

# --------------------------------------------------------------------------
## 2. Repaso ejecutable de control de calidad --------------------------------
#
# > Este bloque hace el control de calidad **de este dataset**. No hereda ningún objeto
#   de ningún lado: en unos 3 minutos tenemos la matriz filtrada desde la que arranca
#   todo el módulo.
# >
# > Usamos **umbrales adaptativos por MAD**, no umbrales fijos del tipo `percent.mt < 10`.
#   La idea es que el umbral se derive de la distribución de *estos* datos, en lugar de
#   ser una constante que uno copia de un tutorial.

pbmc <- CreateSeuratObject(counts = rna, project = "pbmc10k",
                           min.cells = 3, min.features = 200)

# Porcentaje de lecturas mitocondriales y ribosomales
pbmc[["percent.mt"]]   <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.ribo"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

pbmc

### 2.1. Umbrales por MAD ----------------------------------------------------
#
# La **desviación absoluta mediana** (MAD) es un estimador robusto de dispersión: a
# diferencia de la desviación estándar, no se deja arrastrar por los valores extremos,
# que es justo lo que estamos tratando de detectar.
#
# Una célula es *outlier* si está a más de `k` MADs de la mediana. La convención actual
# ([single-cell best
# practices](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html))
# es:
#
# - **5 MADs** para conteos y número de genes (permisivo, para no perder poblaciones
#   raras)
# - **3 MADs** para el porcentaje mitocondrial, **más un tope absoluto** de seguridad
#
# Los valores de `k` siguen siendo una decisión nuestra. La diferencia es que ahora es
# una decisión **explícita y sobre la escala de estos datos**, no un número mágico.

is_outlier <- function(x, nmads = 5, type = c("both", "lower", "higher")) {
  type <- match.arg(type)
  med <- median(x); mad_x <- mad(x)
  lo <- med - nmads * mad_x
  hi <- med + nmads * mad_x
  switch(type, both = x < lo | x > hi, lower = x < lo, higher = x > hi)
}

md <- pbmc@meta.data

# Trabajamos en escala log para los conteos: sus distribuciones son muy asimétricas
out_counts   <- is_outlier(log1p(md$nCount_RNA),   5)
out_features <- is_outlier(log1p(md$nFeature_RNA), 5)
out_mito     <- is_outlier(md$percent.mt, 3, "higher") | md$percent.mt > 20

pbmc$qc_outlier <- out_counts | out_features | out_mito

cat("Umbrales derivados de los datos:\n")
cat("  nCount_RNA   :", round(expm1(median(log1p(md$nCount_RNA)) - 5 * mad(log1p(md$nCount_RNA)))), "-",
                        round(expm1(median(log1p(md$nCount_RNA)) + 5 * mad(log1p(md$nCount_RNA)))), "\n")
cat("  nFeature_RNA :", round(expm1(median(log1p(md$nFeature_RNA)) - 5 * mad(log1p(md$nFeature_RNA)))), "-",
                        round(expm1(median(log1p(md$nFeature_RNA)) + 5 * mad(log1p(md$nFeature_RNA)))), "\n")
cat("  percent.mt   : <", round(median(md$percent.mt) + 3 * mad(md$percent.mt), 2), "% (tope duro 20 %)\n\n")
print(table(Descartadas = pbmc$qc_outlier))


VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "qc_outlier", ncol = 3, pt.size = 0) &
  labs(x = "¿Descartada por QC?")

# Doublets con scDblFinder (sigue siendo el método recomendado en 2026)
suppressPackageStartupMessages(library(scDblFinder))

sce <- tryCatch(
  as.SingleCellExperiment(pbmc),
  error = function(e) {
    message("as.SingleCellExperiment() no pudo usar la capa `data` ",
            "(el objeto aún no está normalizado). Se construye el ",
            "SingleCellExperiment desde los conteos, que es lo único que ",
            "scDblFinder necesita.")
    SingleCellExperiment::SingleCellExperiment(
      list(counts = SeuratObject::LayerData(pbmc, assay = "RNA", layer = "counts")))
  }
)
sce <- scDblFinder(sce, verbose = FALSE)

pbmc$doublet_class <- colData(sce)$scDblFinder.class
pbmc$doublet_score <- colData(sce)$scDblFinder.score
table(pbmc$doublet_class)

# Filtrado final
pbmc <- subset(pbmc, subset = qc_outlier == FALSE & doublet_class == "singlet")
cat("Objeto final:", nrow(pbmc), "genes x", ncol(pbmc), "células\n")
pbmc

# --------------------------------------------------------------------------
## 3. Normalización ----------------------------------------------------------
#
### El problema --------------------------------------------------------------
#
# Comparemos dos células cualquiera:


FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA", pt.size = 0.3) +
  labs(title = "Cada punto es una célula",
       subtitle = "El rango de conteos totales abarca más de un orden de magnitud") +
  NoLegend()

q <- quantile(pbmc$nCount_RNA, c(0.05, 0.5, 0.95))
cat("Conteos totales por célula — percentil 5:", q[1],
    "| mediana:", q[2], "| percentil 95:", q[3], "\n")
cat("La célula del percentil 95 tiene", round(q[3] / q[1], 1),
    "veces más lecturas que la del percentil 5.\n")
cat("¿Es", round(q[3] / q[1], 1), "veces más activa transcripcionalmente? Casi seguro que no.\n")

### 3.1. LogNormalize: el caballito de batalla -------------------------------
#
# $$ \tilde{x}_{gc} = \log\left(1 + \frac{x_{gc}}{\sum_g x_{gc}} \times 10^4\right) $$
#
# Tres pasos:
#
# 1. **Dividir por el total de la célula** → cada célula suma 1 (fracciones, no
#    conteos).
# 2. **Multiplicar por 10.000** (*counts per 10K*). El factor es arbitrario; solo evita
#    trabajar con números diminutos.
# 3. **`log(1 + x)`** → comprime la cola derecha y hace que las diferencias sean
#    multiplicativas, no aditivas.
#
# **Qué supone este método:** que el ARN total por célula es una molestia técnica.
# **Cuándo eso es falso:** cuando el contenido de ARN es la biología que te interesa (un
# plasmablasto tiene genuinamente mucho más ARN que un linfocito en reposo; una célula
# en mitosis, también).
#
# **El pseudoconteo (`+1`)** no es inocente: introduce un sesgo dependiente de la
# profundidad, especialmente en genes de baja expresión. Es el precio de poder tomar
# logaritmo de un cero.

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize",
                      scale.factor = 1e4, verbose = FALSE)

# En Seurat v5 los datos viven en LAYERS, no en slots.
# Esta es la API pública correcta (el viejo pbmc@assays$RNA@counts ya NO funciona
# sobre un objeto Assay5):
Layers(pbmc[["RNA"]])

# Comparación cruda vs normalizada, para el mismo gen y las mismas células
crudo <- pbmc[["RNA"]]$counts["CD3E", 1:8]
norm  <- pbmc[["RNA"]]$data["CD3E", 1:8]
tot   <- pbmc$nCount_RNA[1:8]

data.frame(conteo_crudo = as.numeric(crudo),
           total_celula = as.numeric(tot),
           lognorm      = round(as.numeric(norm), 3))

### 3.2. SCTransform v2: modelar en vez de dividir ---------------------------
#
# `LogNormalize` asume que un solo factor de escala por célula alcanza. `SCTransform`
# toma otro camino: ajusta una **regresión binomial negativa regularizada** gen por gen,
# modelando explícitamente la relación media–varianza, y devuelve **residuos de
# Pearson** como valores "normalizados".
#
# Ventajas: no necesita pseudoconteo, no sobre-corrige los genes de baja expresión, y
# hace en un paso lo que `LogNormalize` + `FindVariableFeatures` + `ScaleData` hacen en
# tres.
#
# Costo: más lento (con `glmGamPoi` es manejable) y los valores resultantes son
# residuos, no expresión — hay que tener cuidado al interpretarlos directamente.
#
# > **Detalle importante:** `SCTransform` crea un **assay nuevo** llamado `"SCT"`. El
#   assay `"RNA"` sigue ahí. Podemos comparar los dos caminos sobre el mismo objeto.

# ~2-3 minutos con glmGamPoi
pbmc <- SCTransform(pbmc, vst.flavor = "v2", method = "glmGamPoi",
                    vars.to.regress = NULL, verbose = FALSE)

Assays(pbmc)

cat("Assay activo por defecto:", DefaultAssay(pbmc), "\n")
cat("Genes en el assay SCT:", nrow(pbmc[["SCT"]]), "(SCTransform recorta a los informativos)\n")
cat("Genes en el assay RNA:", nrow(pbmc[["RNA"]]), "\n")

# --------------------------------------------------------------------------
## 4. Genes variables (HVGs) -------------------------------------------------
#
# El archivo trae ~33.500 genes, pero `CreateSeuratObject(min.cells = 3)` ya descartó los
# que no se detectan en al menos 3 células: quedan ~20.000. De esos, la mayoría no aporta
# nada para distinguir tipos celulares: o no se expresan, o se expresan igual en todas
# las células (genes *housekeeping*).
#
# Nos quedamos con los genes cuya varianza **excede lo que se esperaría dado su nivel
# medio de expresión**. Ese "dado su nivel medio" es la clave: en datos de conteo, los
# genes muy expresados son automáticamente más variables. El método `vst` ajusta esa
# relación media–varianza y busca los genes que se salen de la curva.

DefaultAssay(pbmc) <- "RNA"
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",
                             nfeatures = 2000, verbose = FALSE)

top15 <- head(VariableFeatures(pbmc), 15)
cat("Top 15 genes variables:\n"); print(top15)


p <- VariableFeaturePlot(pbmc)
LabelPoints(p, points = top15, repel = TRUE, xnudge = 0, ynudge = 0) +
  labs(title = "Relación media–varianza",
       subtitle = "En rojo, los 2000 genes seleccionados") +
  theme(legend.position = "bottom")

### 4.1. ¿Cuántos genes? — una decisión, no una constante --------------------
#
# `nfeatures = 2000` es el valor por defecto de Seurat y aparece en prácticamente todos
# los tutoriales. Veamos qué tan estable es esa elección.

tmp <- FindVariableFeatures(pbmc, nfeatures = 500,  verbose = FALSE); hvg500  <- VariableFeatures(tmp)
tmp <- FindVariableFeatures(pbmc, nfeatures = 1000, verbose = FALSE); hvg1000 <- VariableFeatures(tmp)
tmp <- FindVariableFeatures(pbmc, nfeatures = 3000, verbose = FALSE); hvg3000 <- VariableFeatures(tmp)
hvg2000 <- VariableFeatures(pbmc)

cat("Solapamiento con la selección de 2000 genes:\n")
cat("  500  :", length(intersect(hvg500,  hvg2000)), "/ 500\n")
cat("  1000 :", length(intersect(hvg1000, hvg2000)), "/ 1000\n")
cat("  3000 :", length(intersect(hvg3000, hvg2000)), "/ 2000 de los nuestros están incluidos\n\n")

# ¿Marcadores conocidos quedan afuera con 500?
marcadores <- c("CD3E", "CD8A", "MS4A1", "CD14", "FCGR3A", "NKG7", "PPBP", "FCER1A", "IL7R")
data.frame(marcador = marcadores,
           en_500  = marcadores %in% hvg500,
           en_2000 = marcadores %in% hvg2000)

# > **Para discutir:** `CD3E` es *el* marcador de linfocitos T y en muchas selecciones
#   acotadas no aparece entre los genes variables. ¿Es un problema?
# >
# > No necesariamente: los HVGs se usan para **construir el espacio**, no para anotar.
#   Un gen puede ser excelente marcador (alto y consistente en un tipo celular) sin ser
#   de los más variables del dataset. Pero conviene tenerlo presente cuando el
#   clustering "no separa" una población que uno espera ver.

# --------------------------------------------------------------------------
## 5. Escalamiento -----------------------------------------------------------
#
# `ScaleData` centra cada gen en media 0 y lo escala a varianza 1 (z-score por gen).
#
# **Por qué hace falta:** PCA busca direcciones de máxima varianza. Sin escalar, los
# genes muy expresados dominan las primeras componentes solo por su magnitud, no por su
# capacidad de discriminar.
#
# **Qué se pierde:** después de escalar, un gen expresado en el 1 % de las células pesa
# lo mismo que uno expresado en el 60 %.
#
# Por defecto Seurat escala **solo los HVGs**, que es lo que se usa para PCA. Escalar
# los ~20.000 genes solo hace falta si vas a hacer un heatmap con genes que no son HVGs
# — y en un dataset de 10k células consume varios GB de RAM.

pbmc <- ScaleData(pbmc, verbose = FALSE)   # solo HVGs
Layers(pbmc[["RNA"]])

### 5.1. Regresar covariables: cuidado ---------------------------------------
#
# ---------------------------------------------------------------------------
# Nota de vocabulario: en inglés esto se llama *regress out*, y a veces se traduce
# como "regresar la covariable". Es un calco confuso —acá "regresar" no es volver ni
# devolver—, así que en este curso decimos **eliminar por regresión** o **corregir por**.
#
# La operación, gen por gen: se ajusta una regresión de la expresión contra la
# covariable y se conserva SOLO EL RESIDUO, o sea la parte del gen que la covariable
# no explica.
# ---------------------------------------------------------------------------
# `vars.to.regress` permite **eliminar por regresión** el efecto de una covariable (por
# ejemplo `percent.mt`
# o la fase del ciclo celular) antes de PCA.
#
# Es tentador y a veces necesario. Pero es una operación agresiva: si la covariable
# **está correlacionada con la biología** —y el porcentaje mitocondrial lo está: hay
# tipos celulares genuinamente más mitocondriales, como los cardiomiocitos— corregir por ella
# borra señal real.
#
# Regla práctica: **primero mira si la covariable está estructurando el embedding**. Si
# no lo está, no la corrijas.

# NO ejecutamos esto en clase (tarda), queda como referencia:
# pbmc <- ScaleData(pbmc, vars.to.regress = c("percent.mt"), verbose = FALSE)
cat("Ver más abajo el diagnóstico visual antes de decidir si hace falta.\n")

# --------------------------------------------------------------------------
## 6. PCA --------------------------------------------------------------------
#
# PCA es el paso más importante de todos y el que menos atención suele recibir. **Todo
# lo que sigue —vecinos, clusters, UMAP, integración— se calcula sobre el espacio de
# PCA, no sobre los genes.**
#
# Cada componente principal es una combinación lineal de genes. Las primeras capturan
# los ejes de variación dominantes, que en un PBMC típicamente son "linfoide vs
# mieloide", "T vs B", etc.

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 50, verbose = FALSE)

# Qué genes definen las primeras componentes
print(pbmc[["pca"]], dims = 1:4, nfeatures = 8)

DimHeatmap(pbmc, dims = 1:6, cells = 500, balanced = TRUE, ncol = 3)

### 6.1. ¿Cuántas componentes conservar? -------------------------------------
#
# No hay respuesta correcta. Hay criterios, y conviene mirar más de uno.
#
# **a) Codo (`ElbowPlot`).** Desde Seurat 5.5 se puede graficar directamente el
# **porcentaje de varianza explicada**, que es más informativo que la desviación
# estándar cruda.


p1 <- ElbowPlot(pbmc, ndims = 50) + labs(title = "Desviación estándar por PC")

# Varianza explicada, calculada a mano (transparente y siempre funciona)
ev  <- Stdev(pbmc, reduction = "pca")^2
tot <- pbmc[["pca"]]@misc$total.variance
if (is.null(tot)) tot <- sum(apply(GetAssayData(pbmc, layer = "scale.data"), 1, var))
vpc <- 100 * ev / tot
df  <- data.frame(PC = 1:length(vpc), var = vpc, acum = cumsum(vpc))

p2 <- ggplot(df, aes(PC, acum)) +
  geom_line(color = "steelblue") + geom_point(size = 1) +
  geom_hline(yintercept = df$acum[30], linetype = 3, color = "red") +
  labs(title = "Varianza acumulada", y = "% acumulado") + theme_bw()

p1 | p2

cat("Varianza explicada por las primeras PCs:\n")
print(round(head(df, 10), 2))
cat("\nCon 10 PCs:", round(df$acum[10], 1), "% | 20 PCs:", round(df$acum[20], 1),
    "% | 30 PCs:", round(df$acum[30], 1), "%\n")

# **Un dato incómodo:** con 30 componentes solemos explicar apenas un 20–30 % de la
# varianza total. Eso no significa que estemos tirando el 70 % de la biología — buena
# parte de esa varianza restante es ruido de muestreo (shot noise) de los conteos.
#
# **b) Criterio práctico.** En la mayoría de los datasets de PBMC, entre 20 y 30 PCs
# funciona bien y el resultado es robusto a esa elección. **Errar por exceso es más
# seguro que por defecto:** las PCs de más aportan ruido, pero las PCs de menos borran
# poblaciones enteras.
#
# **c) Diagnóstico técnico:** ¿alguna PC está capturando algo que no queremos?


FeaturePlot(pbmc, reduction = "pca",
            features = c("percent.mt", "nCount_RNA"), ncol = 2) &
  theme(plot.title = element_text(size = 11))

# ¿Las primeras PCs correlacionan con covariables técnicas?
emb <- Embeddings(pbmc, "pca")[, 1:10]
cor_tec <- data.frame(
  PC          = paste0("PC", 1:10),
  r_nCount    = round(apply(emb, 2, cor, pbmc$nCount_RNA), 2),
  r_percentMT = round(apply(emb, 2, cor, pbmc$percent.mt), 2)
)
cor_tec

# > **Cómo leer esta tabla:** una correlación |r| > 0.5 entre una PC de las primeras y
#   `nCount_RNA` o `percent.mt` es una señal de alarma. Significa que esa componente
#   está codificando profundidad de secuenciación o estrés celular, no biología, y que
#   los clusters que salgan de ahí van a reflejar calidad técnica.
# >
# > Si eso pasa: revisar el QC, considerar `vars.to.regress`, o usar `SCTransform` (que
#   modela la profundidad explícitamente).

# --------------------------------------------------------------------------
## 7. t-SNE y UMAP -----------------------------------------------------------
#
# PCA es lineal: fácil de interpretar, pero malo para visualizar estructura no lineal en
# 2D. Para eso usamos métodos de *embedding* no lineal.
#
# | | t-SNE | UMAP |
# |---|---|---|
# | Preserva | estructura local | local, algo de global |
# | Velocidad | lento | rápido |
# | Parámetro clave | `perplexity` | `n.neighbors`, `min.dist` |
# | Reproducible | solo con misma semilla | solo con misma semilla |
#
# Ambos se calculan **sobre el espacio de PCA**, no sobre los genes: la elección de
# `dims` de la sección anterior se propaga aquí.

pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- RunTSNE(pbmc, dims = 1:30, verbose = FALSE, check_duplicates = FALSE)

# ---------------------------------------------------------------------------
# Un color para poder mirar los dibujos
#
# Todavía no clusterizamos —eso es la clase que viene— así que un `DimPlot` sin
# más sale de un solo color, y así no se puede comparar nada entre proyecciones.
# Necesitamos una etiqueta por célula, y la vamos a construir del modo más
# transparente posible: el promedio de unos pocos marcadores canónicos por linaje,
# y a cada célula le ponemos el linaje con el puntaje más alto.
#
# ⚠️ Esto NO es anotación celular (eso es la clase 3) ni clustering (clase 2). Es
#    una etiqueta gruesa y provisoria, solo para poder ver. Pero fíjate en una cosa
#    que importa para el argumento de esta clase: la etiqueta viene de marcadores
#    conocidos DE ANTEMANO, no de haber clusterizado estos datos. Por eso es legítimo
#    usarla para juzgar las proyecciones — no hay circularidad.

marcadores_linaje <- list(
  "Linfocitos T / NK" = c("CD3E", "CD3D", "IL7R", "TRAC", "NKG7", "GNLY"),
  "Linfocitos B"      = c("MS4A1", "CD79A", "CD79B", "BANK1"),
  "Mieloides"         = c("LYZ", "S100A9", "CD14", "FCN1", "FCGR3A")
)

datos_norm <- GetAssayData(pbmc, assay = "RNA", layer = "data")
puntajes <- sapply(marcadores_linaje, function(genes) {
  genes <- intersect(genes, rownames(datos_norm))
  Matrix::colMeans(datos_norm[genes, , drop = FALSE])
})

# Si ningún linaje supera un puntaje mínimo, la célula queda "sin asignar" en vez
# de forzarle una etiqueta. Preferimos un gris honesto a un color inventado.
mejor <- apply(puntajes, 1, max)
pbmc$linaje <- factor(
  ifelse(mejor < 0.05, "sin asignar", colnames(puntajes)[max.col(puntajes)]),
  levels = c(names(marcadores_linaje), "sin asignar")
)

# La misma paleta que se usa en las diapositivas, para que lo proyectado y lo que
# te sale en pantalla sean la misma figura.
COL_LINAJE <- c("Linfocitos T / NK" = "#2a78d6",
                "Linfocitos B"      = "#1baf7a",
                "Mieloides"         = "#eb6834",
                "sin asignar"       = "#C9CFD6")

print(table(pbmc$linaje))
# ---------------------------------------------------------------------------

(DimPlot(pbmc, reduction = "umap", group.by = "linaje", cols = COL_LINAJE) +
   ggtitle("UMAP (30 PCs)") + NoLegend()) |
(DimPlot(pbmc, reduction = "tsne", group.by = "linaje", cols = COL_LINAJE) +
   ggtitle("t-SNE (30 PCs)"))

### 7.1. ⚠️ Qué NO se puede leer en un UMAP ----------------------------------
#
# Esta es la parte más importante de la clase. Un UMAP es un dibujo, y la mayoría de las
# cosas que uno quiere leer en él no están ahí:
#
# | Uno tiende a leer... | ¿Es válido? |
# |---|---|
# | "Estas dos islas están lejos → son muy distintas" | ❌ **No.** Las distancias entre grupos no tienen significado cuantitativo. |
# | "Esta isla es más grande → hay más células o más heterogeneidad" | ❌ **No.** UMAP no preserva densidad ni tamaño. |
# | "Este grupo está en el medio → es intermedio / de transición" | ❌ **No** por sí solo. Hay que verificarlo en el espacio de expresión. |
# | "Estas células están juntas → se parecen entre sí" | ✅ **Sí**, esto es lo que UMAP sí preserva (estructura local). |
# | "Hay un puente continuo entre A y B" | ⚠️ **Sugerente**, no probatorio. Puede ser un artefacto de doublets o de ARN ambiental (clase 3). |
#
# Demostrémoslo cambiando solo los hiperparámetros, sin tocar los datos.
#
# Esta es la figura de los seis UMAPs de la presentación teórica. Aquí la
# reproduces tú:
# los mismos datos, el mismo espacio de PCA, y una grilla de 2 x 3. Lo único que
# cambia entre paneles son dos números:
#
# - **filas**: cuántas componentes principales entran al UMAP (10 y 30);
# - **columnas**: cuántos vecinos considera el algoritmo (`n.neighbors` = 5, 30 y 100).
#
# `min.dist` queda fijo en 0.3 (el valor por defecto) para que la comparación
# aísle un factor por eje. La semilla también es fija: lo que veas entre paneles
# es el efecto del hiperparámetro, no del azar.
#
# Ojo con el tiempo: son seis UMAPs y cada uno tarda unos segundos. Si tienes
# cache activado en el .qmd, se calcula una sola vez.

grilla_pcs <- c(10, 30)          # filas
grilla_nn  <- c(5, 30, 100)      # columnas

plots <- list()
for (d in grilla_pcs) {
  for (nn in grilla_nn) {
    u <- RunUMAP(pbmc, dims = 1:d, n.neighbors = nn, min.dist = 0.3,
                 reduction.name = "tmpumap", seed.use = 42, verbose = FALSE)
    plots[[length(plots) + 1]] <-
      DimPlot(u, reduction = "tmpumap", group.by = "linaje", cols = COL_LINAJE,
              pt.size = 0.2) +
      ggtitle(paste0(d, " PCs  ·  n.neighbors = ", nn)) + NoLegend() +
      theme(plot.title = element_text(size = 10),
            axis.title = element_text(size = 8))
  }
}

# Las seis proyecciones, en el mismo orden que en la diapositiva. El panel de
# abajo al centro (30 PCs, n.neighbors = 30) es exactamente el UMAP que
# calculamos en la sección 7: sirve de ancla para comparar los otros cinco.
wrap_plots(plots, nrow = 2, ncol = 3)

# El color es la etiqueta de linaje que construimos arriba, la misma en los seis
# paneles y la misma que está proyectada. Eso es lo que hace legible el argumento:
# **las tres poblaciones aparecen en los seis**. Lo que cambia es cuánto se
# separan, qué forma tienen y cuánto se fragmenta por dentro el compartimento
# T/NK — que es la parte del dataset con estructura continua, y por eso la más
# sensible al parámetro.

# > **Ninguna de estas seis figuras es la verdadera.** Todas son proyecciones legítimas
#   de los mismos datos.
# >
# > Por eso, en un paper, un UMAP sin la información de cuántas PCs y qué
#   hiperparámetros se usaron es una figura **irreproducible**. Y por eso el UMAP se usa
#   para *comunicar* un resultado, nunca para *derivarlo*: las conclusiones
#   cuantitativas se sacan del espacio de PCA y de los conteos, no del dibujo.

# Diagnóstico final: ¿hay estructura técnica dominando el embedding?

FeaturePlot(pbmc, features = c("percent.mt", "nCount_RNA", "doublet_score"),
            reduction = "umap", ncol = 3) &
  theme(plot.title = element_text(size = 11))

# Y un vistazo a marcadores canónicos: ya se adivinan las poblaciones

FeaturePlot(pbmc, reduction = "umap", ncol = 3,
            features = c("CD3E", "CD8A", "IL7R", "CD14", "FCGR3A", "MS4A1",
                         "NKG7", "FCER1A", "PPBP")) &
  theme(plot.title = element_text(size = 11)) & NoLegend()

# --------------------------------------------------------------------------
## 8. Guardar el checkpoint --------------------------------------------------
#
# Guardamos el objeto para retomarlo en la clase 2. Queda en `results/` dentro del
# proyecto, así que la clase 2 lo va a encontrar sola: no hay que subir ni bajar nada.
# El archivo pesa varios cientos de MB — no lo pongas bajo control de versiones.

saveRDS(pbmc, ck1)
tam_archivo(ck1)

# Registro de la sesión: esto va SIEMPRE en el material suplementario de un paper
sessionInfo()

# --------------------------------------------------------------------------
## 9. Síntesis ---------------------------------------------------------------
#
# 1. **Normalizar** corrige diferencias de profundidad, no diferencias de contenido
#    celular de ARN. `LogNormalize` es simple y robusto; `SCTransform` v2 modela la
#    relación media–varianza y evita el pseudoconteo.
# 2. **Los HVGs** definen el espacio, no la anotación. Un marcador excelente puede no
#    ser un HVG.
# 3. **Escalar** iguala el peso de los genes; eliminar covariables por regresión es potente y
#    peligroso a la vez.
# 4. **PCA es el paso decisivo**: todo lo posterior vive en ese espacio. Errar por
#    exceso de PCs es más seguro que por defecto.
# 5. **UMAP y t-SNE son mapas, no mediciones.** Preservan vecindad local y nada más.
#
# En la próxima clase vamos a cortar este espacio en clusters — y a ver por qué los
# p-valores que Seurat devuelve para los marcadores de esos clusters **no son p-valores
# válidos**.
