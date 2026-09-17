# =============================================================================
# Bioinformática y Estadística 3 — Módulo scRNA-seq + CITE-seq
# Clase 2 — Clustering y expresión diferencial
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
## Módulo scRNA-seq + CITE-seq · Clase 2 de 4 --------------------------------
### Clustering y expresión diferencial: qué se puede afirmar y qué no --------
#
# **Docente:** Dr. Danilo Ceschin · **Fecha:** jueves 10 de septiembre de 2026 ·
# **Duración:** 2 h
#
# --------------------------------------------------------------------------
#
# En la clase 1 construimos un espacio de baja dimensión. Hoy lo cortamos en grupos y le
# buscamos genes característicos a cada grupo.
#
# Esta clase tiene dos mitades muy distintas. La primera es mecánica: cómo se construye
# el grafo, qué hace el algoritmo de comunidades, qué significa la resolución. La
# segunda es estadística, y es la razón de que esta materia se llame *Bioinformática **y
# Estadística***:
#
# > Cuando definimos clusters con ciertos datos y después testeamos diferencias entre
#   esos clusters con **los mismos datos**, los p-valores que obtenemos no son
#   p-valores. Este procedimiento —el que aparece en casi todos los tutoriales de
#   scRNA-seq— produce falsos positivos con una tasa que no controlamos.
#
# Vamos a demostrarlo empíricamente y a ver qué se hace en su lugar.
#
# --------------------------------------------------------------------------
#
### Objetivos de la clase ----------------------------------------------------
#
# 1. Construir el grafo KNN/SNN y aplicar algoritmos de detección de comunidades
#    (Louvain, Leiden).
# 2. Entender la resolución como parámetro sin valor "correcto" y evaluar la estabilidad
#    del clustering.
# 3. **Reconocer el problema de doble inmersión (*double dipping*)** y sus
#    consecuencias.
# 4. Aplicar el flujo correcto para comparar condiciones: **pseudobulk** con réplicas
#    biológicas.
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
## 1. Retomamos el objeto de la clase 1 --------------------------------------
#
# Dos opciones: si corriste la clase 1, el checkpoint ya está en `results/` y el
# primer bloque lo carga solo. Si no, el segundo bloque regenera todo en ~4 minutos.

# OPCIÓN A — ya tienes el checkpoint de la clase 1
if (file.exists(ck1)) {
  pbmc <- readRDS(ck1)
  DefaultAssay(pbmc) <- "RNA"
  cat("Objeto cargado:", ncol(pbmc), "células\n")
} else {
  cat("No encontrado. Ejecuta el bloque de abajo (OPCIÓN B).\n")
}

# OPCIÓN B — regenerar desde cero (~4 min)
if (!exists("pbmc")) {
  suppressPackageStartupMessages(library(scDblFinder))

  url_filt <- paste0("https://cf.10xgenomics.com/samples/cell-exp/3.0.0/",
                     "pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")
  if (!file.exists(f_h5_filt))
    download.file(url_filt, f_h5_filt, mode = "wb", quiet = TRUE)

  data <- Read10X_h5(f_h5_filt)
  pbmc <- CreateSeuratObject(data[["Gene Expression"]], project = "pbmc10k",
                             min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

  is_outlier <- function(x, k = 5, type = "both") {
    med <- median(x); m <- mad(x)
    if (type == "higher") x > med + k * m else x < med - k * m | x > med + k * m
  }
  keep <- !(is_outlier(log1p(pbmc$nCount_RNA)) | is_outlier(log1p(pbmc$nFeature_RNA)) |
            is_outlier(pbmc$percent.mt, 3, "higher") | pbmc$percent.mt > 20)
  pbmc <- pbmc[, keep]

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
  pbmc <- pbmc[, colData(sce)$scDblFinder.class == "singlet"]

  pbmc <- NormalizeData(pbmc, verbose = FALSE) |>
          FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
          ScaleData(verbose = FALSE) |>
          RunPCA(npcs = 50, verbose = FALSE) |>
          RunUMAP(dims = 1:30, verbose = FALSE)
  cat("Objeto regenerado:", ncol(pbmc), "células\n")
}
pbmc

# --------------------------------------------------------------------------
## 2. Clustering basado en grafos --------------------------------------------
#
### 2.1. Por qué no usamos k-means -------------------------------------------
#
# k-means, clustering jerárquico y compañía suponen clusters aproximadamente esféricos y
# de tamaño similar. Las poblaciones celulares no cumplen ninguna de las dos cosas: hay
# tipos con 4.000 células y tipos con 40, y muchos tienen forma alargada porque
# contienen un gradiente de estados.
#
# El clustering por grafos no hace supuestos sobre la forma. Son tres pasos:
#
# 1. **Grafo KNN** — cada célula se conecta con sus `k` vecinos más cercanos **en el
#    espacio de PCA**.
# 2. **Grafo SNN** (*shared nearest neighbors*) — el peso de la arista entre dos células
#    pasa a ser cuántos vecinos comparten. Esto suaviza el ruido: dos células conectadas
#    por casualidad rara vez comparten muchos vecinos.
# 3. **Detección de comunidades** — se busca la partición que maximiza la
#    **modularidad**: muchas aristas dentro de los grupos, pocas entre grupos.
#
# > Nota que la elección de `dims` de la clase 1 vuelve a aparecer aquí, y esta vez
#   determina de verdad quién es vecino de quién.

pbmc <- FindNeighbors(pbmc, dims = 1:30, k.param = 20, verbose = FALSE)
Graphs(pbmc)

### 2.2. Louvain vs Leiden ---------------------------------------------------
#
# `FindClusters` ofrece cuatro algoritmos:
#
# | `algorithm` | Método |
# |---|---|
# | 1 | Louvain (por defecto) |
# | 2 | Louvain con refinamiento multinivel |
# | 3 | SLM |
# | 4 | **Leiden** |
#
# **Leiden corrige un defecto real de Louvain:** Louvain puede producir comunidades
# *internamente desconectadas* — grupos que el algoritmo declara como una unidad pero
# que en el grafo no están conectados entre sí. Leiden garantiza que eso no pase.
#
# Históricamente Leiden en Seurat requería Python (`reticulate` + `leidenalg`), lo cual
# era una fuente constante de problemas de instalación. **Desde Seurat 5.3.1 se puede
# usar la implementación de `igraph`, sin Python:**
#
# ```r
# FindClusters(obj, algorithm = 4, leiden_method = "igraph")
# ```

# `cluster.name` guarda cada partición en su propia columna. Sin ese argumento las dos
# corridas escribirían en `RNA_snn_res.0.5` y la segunda pisaría a la primera — un error
# silencioso: el barrido de la sección siguiente quedaría mezclando algoritmos.
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1,
                     cluster.name = "louvain", verbose = FALSE)

pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 4,
                     leiden_method = "igraph",
                     cluster.name = "leiden", verbose = FALSE)

cat("Louvain:", nlevels(pbmc$louvain), "clusters\n")
cat("Leiden :", nlevels(pbmc$leiden),  "clusters\n")

(DimPlot(pbmc, group.by = "louvain", label = TRUE) + ggtitle("Louvain (res 0.5)") + NoLegend()) |
(DimPlot(pbmc, group.by = "leiden",  label = TRUE) + ggtitle("Leiden (res 0.5)")  + NoLegend())

# ¿En qué se diferencian? Tabla de contingencia
tab <- table(Louvain = pbmc$louvain, Leiden = pbmc$leiden)
tab

# Índice de Rand ajustado: 1 = particiones idénticas, 0 = coincidencia al azar
ari <- function(x, y) {
  tb <- table(x, y)
  n  <- sum(tb)
  ch2 <- function(k) k * (k - 1) / 2
  s_ij <- sum(ch2(tb)); s_i <- sum(ch2(rowSums(tb))); s_j <- sum(ch2(colSums(tb)))
  exp_i <- s_i * s_j / ch2(n)
  (s_ij - exp_i) / ((s_i + s_j) / 2 - exp_i)
}
cat("ARI Louvain vs Leiden:", round(ari(pbmc$louvain, pbmc$leiden), 3), "\n")

# ---------------------------------------------------------------------------
# La comparación con k-means, que es la figura proyectada en la clase.
# Para que se vea el problema hay que mirar una partición GRUESA: pedimos pocos
# grupos y comparamos qué hace cada método con la población más chica.
km_k <- 5
pbmc <- FindClusters(pbmc, resolution = 0.05, algorithm = 4,
                     leiden_method = "igraph",
                     cluster.name = "grafo_grueso", verbose = FALSE)

set.seed(42)
pbmc$kmeans <- factor(kmeans(Embeddings(pbmc, "pca")[, 1:30],
                             centers = km_k, nstart = 10)$cluster)

t_grafo <- sort(table(pbmc$grafo_grueso), decreasing = TRUE)
t_km    <- sort(table(pbmc$kmeans),       decreasing = TRUE)
cat("Grupo más chico  — grafo :", min(t_grafo), "células\n")
cat("Grupo más chico  — k-means:", min(t_km),   "células\n")
cat("Razón mayor/menor — grafo :", round(max(t_grafo)/min(t_grafo)), "x\n")
cat("Razón mayor/menor — k-means:", round(max(t_km)/min(t_km)),     "x\n")
cat("ARI grafo vs k-means:", round(ari(pbmc$grafo_grueso, pbmc$kmeans), 3), "\n")

(DimPlot(pbmc, group.by = "kmeans") + ggtitle(paste0("k-means (k = ", km_k, ")")) + NoLegend()) |
(DimPlot(pbmc, group.by = "grafo_grueso") + ggtitle("Grafo SNN + Leiden") + NoLegend())

# > Mirá el bloque grande: k-means lo parte en tajadas donde no hay ninguna
#   separación, porque busca celdas convexas de tamaño parecido. Y mirá la
#   población más chica: el grafo le da un cluster propio, k-means la absorbe.
#   En un experimento real esa población podría ser justo la que buscás.
#
# > Advertencia honesta: con muchos grupos y poblaciones bien separadas, k-means
#   no lo hace tan mal. El problema aparece donde importa — poblaciones chicas y
#   nubes alargadas.
# ---------------------------------------------------------------------------

# --------------------------------------------------------------------------
## 3. La resolución ----------------------------------------------------------
#
# La resolución controla la granularidad: valores bajos dan pocos clusters grandes,
# valores altos dan muchos clusters chicos.
#
# **No existe la resolución correcta.** La resolución define la pregunta que estás
# haciendo. "¿Cuántos tipos celulares hay?" no tiene una respuesta única: depende de si
# quieres distinguir *linfocito T* de *monocito*, o *T CD4 naive* de *T CD4 de memoria
# central*.
#
# Lo que sí se puede hacer es elegir con criterio y **reportar la elección**.

# Dos decisiones deliberadas:
#  - 0.5 está EN el barrido. La resolución que vamos a usar el resto del módulo tiene que
#    ser un nodo del árbol que estamos por mirar; si no, la justificamos con un gráfico
#    que no la contiene.
#  - El barrido corre con Leiden, el mismo algoritmo que recomendamos arriba. Mezclar
#    Louvain en el árbol y Leiden en la resolución de trabajo sería incoherente.
# Son nueve corridas: tarda alrededor de un minuto.
resoluciones <- c(0.05, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1.2, 2.0)
pbmc <- FindClusters(pbmc, resolution = resoluciones, algorithm = 4,
                     leiden_method = "igraph", verbose = FALSE)

n_clusters <- sapply(paste0("RNA_snn_res.", resoluciones),
                     function(x) nlevels(factor(pbmc[[x]][, 1])))
data.frame(resolucion = resoluciones, n_clusters = as.integer(n_clusters))

suppressPackageStartupMessages(library(clustree))

clustree(pbmc, prefix = "RNA_snn_res.") +
  labs(title = "Árbol de clusters a través de las resoluciones")

### Cómo se lee un clustree --------------------------------------------------
#
# - Cada fila es una resolución; cada nodo, un cluster.
# - Las **flechas** muestran de dónde vienen las células cuando aumenta la resolución.
# - Un cluster que se divide limpiamente en dos y esos dos se mantienen estables hacia
#   abajo → **subestructura real**.
# - Un nodo que recibe flechas de **varios** clusters de la fila anterior (mezcla de
#   colores en las aristas) → el algoritmo está *reorganizando* células, no *refinando*
#   grupos. Eso es señal de sobre-clusterización.
#
# Regla práctica: quédate en la resolución más alta **antes** de que aparezcan cruces
# masivos.

wrap_plots(lapply(c(0.1, 0.5, 1.2), function(r) {
  col <- paste0("RNA_snn_res.", r)
  DimPlot(pbmc, group.by = col, label = TRUE, pt.size = 0.15) +
    ggtitle(paste("resolución", r)) + NoLegend() +
    theme(plot.title = element_text(size = 11))
}), ncol = 3)

# Fijamos la resolución de trabajo para el resto del módulo
Idents(pbmc) <- "RNA_snn_res.0.5"
pbmc$clusters <- Idents(pbmc)
cat("Resolución de trabajo: 0.5 →", nlevels(pbmc$clusters), "clusters\n")
table(pbmc$clusters)

# --------------------------------------------------------------------------
## 4. ⚠️ El problema estadístico: doble inmersión (double dipping) -----------
#
# Este es el núcleo de la clase.
#
## 4.1. Qué hicimos, exactamente ---------------------------------------------
#
# 1. Usamos la expresión de ~2.000 genes para construir el grafo.
# 2. Buscamos la partición que **maximiza la separación** entre grupos.
# 3. Ahora vamos a testear si esos genes están diferencialmente expresados **entre esos
#    mismos grupos**.
#
# El paso 2 eligió la partición que hace máxima la diferencia. El paso 3 pregunta si la
# diferencia es significativa. **Es circular.** La hipótesis nula que el test de
# Wilcoxon evalúa ("estas dos muestras vienen de la misma distribución") ya fue violada
# por construcción, no por la biología.
#
## 4.2. Demostración ---------------------------------------------------------
#
# Simulemos una matriz de **ruido puro**: conteos de Poisson con una tasa fija por gen,
# igual para todas las células, sin ninguna estructura de grupos. Por construcción, **no hay ningún gen diferencialmente
# expresado, porque no hay grupos**. Después clusterizamos y corremos el test.

set.seed(2026)
n_celulas <- 2000     # con 800 el efecto existe pero queda al borde de lo visible
n_genes   <- 1500

# Conteos de Poisson con una tasa fija por gen (igual para TODAS las células).
# No hay grupos, no hay genes diferenciales, no hay estructura: solo ruido de
# muestreo, que es exactamente el ruido que tienen los datos reales.
tasas <- rlnorm(n_genes, meanlog = 0.5, sdlog = 1)
ruido <- matrix(rpois(n_genes * n_celulas, lambda = tasas),
                nrow = n_genes, ncol = n_celulas)
rownames(ruido) <- paste0("gen", seq_len(n_genes))
colnames(ruido) <- paste0("cel", seq_len(n_celulas))

# Le aplicamos EXACTAMENTE el mismo pipeline que a los datos reales
sim <- CreateSeuratObject(counts = ruido)
sim <- NormalizeData(sim, verbose = FALSE)
sim <- FindVariableFeatures(sim, nfeatures = 1000, verbose = FALSE)
sim <- ScaleData(sim, verbose = FALSE)
sim <- RunPCA(sim, npcs = 20, verbose = FALSE)
sim <- FindNeighbors(sim, dims = 1:20, verbose = FALSE)
sim <- FindClusters(sim, resolution = 0.5, verbose = FALSE)

# La demostración necesita AL MENOS DOS grupos para poder testear. Sobre ruido puro,
# cuántos aparecen depende de la versión de Seurat y del azar; si sale uno solo,
# subimos la resolución hasta que haya dos. Que haga falta subirla no debilita el
# argumento: lo refuerza. El clustering entrega grupos si se los pedís.
#
# ⚠️ El paso importa. Buscamos el corte MÁS GRUESO que el algoritmo produce: dos
#    grupos grandes. Si subimos de a saltos (0.5 → 1 → 2) pasamos de largo esa
#    zona y aterrizamos en seis o siete grupos chicos; el contraste entre los dos
#    mayores queda entonces con pocas células de cada lado y el efecto se diluye.
#    Por eso subimos de a 0.1, desde abajo.
# Primero, el dato que importa para el argumento: cuántos grupos aparecen a la
# resolución por defecto, la que uno usaría sin pensarlo.
cat("Clusters a resolución 0.5, sobre datos SIN ninguna estructura:",
    nlevels(sim$seurat_clusters), "\n")

# Para el TEST buscamos otra cosa: el corte más grueso, dos grupos grandes.
res_sim <- 0.1
sim <- FindClusters(sim, resolution = res_sim, verbose = FALSE)
while (nlevels(sim$seurat_clusters) < 2 && res_sim < 2) {
  res_sim <- round(res_sim + 0.1, 2)
  sim <- FindClusters(sim, resolution = res_sim, verbose = FALSE)
}
cat("Resolución más baja que parte el ruido en dos grupos:", res_sim, "\n")

cat("Grupos que vamos a comparar:", nlevels(sim$seurat_clusters), "\n")
print(table(sim$seurat_clusters))

sim <- RunUMAP(sim, dims = 1:20, verbose = FALSE)
DimPlot(sim, label = TRUE) +
  ggtitle("UMAP de conteos de Poisson sin ninguna estructura") +
  labs(subtitle = "Estos 'grupos' no existen: todas las células vienen de la MISMA distribución") +
  NoLegend()

# > Detengámonos aquí un momento. Ese UMAP tiene grupos separados, con forma, con
#   estructura. Y **los datos son ruido puro**.
# >
# > El clustering siempre devuelve clusters. Nunca dice "no hay grupos".
#
# Ahora la parte que importa: ¿cuántos "marcadores significativos" encuentra el test?

# OJO con `return.thresh`: su valor por defecto es 0.01 y FindAllMarkers descarta
# de entrada todo lo que tenga p-valor CRUDO mayor. Si lo dejamos así, estaríamos
# contando sobre una tabla ya filtrada —y si ningún gen pasa ese filtro, la función
# devuelve un data.frame vacío SIN COLUMNAS, que rompe cualquier indexado posterior.
# Lo ponemos en 1 para recibir todo y hacer nosotros el conteo.
marcadores_ruido <- FindAllMarkers(sim, only.pos = TRUE, logfc.threshold = 0,
                                   min.pct = 0, return.thresh = 1, verbose = FALSE)

cat("Genes testeados en total:", nrow(marcadores_ruido), "\n")

if (nrow(marcadores_ruido) == 0) {
  cat("FindAllMarkers no devolvió ninguna fila. Revisa cuántos clusters hay:\n")
  print(table(sim$seurat_clusters))
} else {
  cat("Genes con p_val_adj < 0.05 :", sum(marcadores_ruido$p_val_adj < 0.05), "\n")
  cat("Genes con p_val_adj < 0.01 :", sum(marcadores_ruido$p_val_adj < 0.01), "\n")
  cat("Genes con p_val_adj < 1e-10:", sum(marcadores_ruido$p_val_adj < 1e-10), "\n\n")
  cat("Recordatorio: la respuesta correcta es CERO.\n\n")

  cat("p-valor ajustado más chico encontrado:",
      format(min(marcadores_ruido$p_val_adj), digits = 3), "\n\n")

  print(head(marcadores_ruido[order(marcadores_ruido$p_val_adj),
                              c("cluster", "gene", "avg_log2FC", "p_val_adj")], 10))
}

# ---------------------------------------------------------------------------
# La comparación decisiva
#
# ⚠️ Para que esta comparación pruebe algo, los dos paneles tienen que diferir en UNA
# sola cosa. Si de un lado usáramos `FindAllMarkers` (uno-contra-el-resto, con
# `only.pos = TRUE`, que además trunca la distribución de p-valores) y del otro un
# contraste pareado sin filtro de signo, estaríamos cambiando tres cosas a la vez:
# el origen de los grupos, la estructura del contraste y el filtro. La conclusión
# no sería atribuible al double dipping.
#
# Por eso los dos paneles usan `FindMarkers`, contraste pareado, sin filtros.
# Lo único que cambia es DE DÓNDE salieron las dos etiquetas.
# ---------------------------------------------------------------------------

# (a) Dos grupos DERIVADOS de estos mismos datos: los dos clusters más grandes
Idents(sim) <- "seurat_clusters"
dos <- names(sort(table(sim$seurat_clusters), decreasing = TRUE))[1:2]
cat("Comparamos los clusters", dos[1], "y", dos[2], "\n")
pv_clusters <- FindMarkers(sim, ident.1 = dos[1], ident.2 = dos[2],
                           logfc.threshold = 0, min.pct = 0, verbose = FALSE)

# (b) Dos grupos asignados AL AZAR, sin mirar los datos
sim$azar <- factor(sample(c("A", "B"), ncol(sim), replace = TRUE))
Idents(sim) <- "azar"
pv_azar <- FindMarkers(sim, ident.1 = "A", ident.2 = "B",
                       logfc.threshold = 0, min.pct = 0, verbose = FALSE)

cat("Genes con p < 0.05 — clusters:", sum(pv_clusters$p_val < 0.05),
    "| al azar:", sum(pv_azar$p_val < 0.05), "\n")
cat("(con ~", nrow(pv_azar), "genes testeados, al azar se esperan ~",
    round(0.05 * nrow(pv_azar)), "por puro azar)\n\n")


p1 <- ggplot(pv_clusters, aes(p_val)) +
  geom_histogram(bins = 50, fill = "firebrick", alpha = .8) +
  labs(title = paste0("Cluster ", dos[1], " vs cluster ", dos[2], " (grupos derivados de los datos)"),
       subtitle = "Bajo H0 verdadera debería ser UNIFORME", x = "p-valor") + theme_bw()

p2 <- ggplot(pv_azar, aes(p_val)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = .8) +
  labs(title = "Grupo A vs grupo B (etiquetas al azar)",
       subtitle = "Aquí sí: uniforme, como debe ser", x = "p-valor") + theme_bw()

p1 | p2

### 4.3. Qué muestra esta comparación ----------------------------------------
#
# - **Izquierda:** los dos grupos salieron del clustering de estos mismos datos → los
#   p-valores se apilan contra cero. Falsos positivos masivos.
# - **Derecha:** los dos grupos se sortearon sin mirar los datos → distribución uniforme,
#   exactamente lo que predice la teoría bajo la hipótesis nula.
#
# Mismo dataset, mismo test, mismo tipo de contraste, mismo número de células por grupo.
# **La única diferencia es si las etiquetas se derivaron o no de los datos que después
# se testean.** Esa diferencia sola produce los dos histogramas.
#
### 4.4. Entonces, ¿qué hacemos? ---------------------------------------------
#
# | Situación | Qué hacer |
# |---|---|
# | Quiero **ordenar** genes candidatos a marcador de un cluster | `FindAllMarkers` está bien. Usa el ranking, **no reportes los p-valores como evidencia inferencial**. |
# | Quiero **testear** si un cluster es un grupo real | `ClusterDE` (control sintético) o *count splitting* (`countsplit`): partir los conteos en dos mitades independientes, clusterizar con una y testear con la otra. |
# | Quiero comparar **condiciones** (enfermo vs sano) dentro de un tipo celular | **Pseudobulk**. Ver la sección 6. |
# | Quiero validar un marcador | Datos independientes: otro dataset, citometría, IHC. |
#
# > **Regla operativa de este curso:** en scRNA-seq, un p-valor de `FindAllMarkers` es
#   una **medida de ranking**, no una probabilidad de error de tipo I. Si en tu tesis
#   escribes "*el gen X está significativamente sobreexpresado en el cluster 3 (p =
#   2×10⁻³⁰⁰)*", esa frase no significa lo que parece.

# --------------------------------------------------------------------------
## 5. Marcadores de cluster, usados correctamente ----------------------------
#
# Con la advertencia anterior en mente, buscamos genes característicos de cada cluster.
# `presto` (instalado en el setup) hace el test de Wilcoxon unas 100 veces más rápido;
# Seurat lo detecta solo.
#
# Los parámetros que importan más que el p-valor:
#
# - `only.pos = TRUE` — solo genes **enriquecidos** en el cluster (los marcadores útiles
#   son positivos).
# - `min.pct` — el gen debe detectarse en al menos ese % de células. Filtra genes que
#   "distinguen" por dropout.
# - `logfc.threshold` — magnitud mínima del efecto. **Este es el criterio que de verdad
#   importa.**

Idents(pbmc) <- "clusters"

marcadores <- FindAllMarkers(pbmc, only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.5,
                             verbose = FALSE)

top5 <- marcadores |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 5) |>
  ungroup()

print(as.data.frame(top5[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2")]), row.names = FALSE)

DoHeatmap(pbmc, features = top5$gene, size = 3) + NoLegend()

# El DotPlot suele ser más informativo que el heatmap:
# el TAMAÑO del punto = % de células que expresan; el COLOR = nivel medio.
genes_dot <- unique(marcadores |> group_by(cluster) |> slice_max(avg_log2FC, n = 3) |> pull(gene))

DotPlot(pbmc, features = genes_dot) +
  RotatedAxis() +
  labs(title = "Marcadores por cluster",
       subtitle = "Tamaño = % de células que lo expresan · Color = expresión media")

# Volcano de un contraste puntual entre dos clusters
c_a <- levels(pbmc$clusters)[1]; c_b <- levels(pbmc$clusters)[2]
de  <- FindMarkers(pbmc, ident.1 = c_a, ident.2 = c_b, verbose = FALSE)
de$gene <- rownames(de)

ggplot(de, aes(avg_log2FC, -log10(p_val_adj + 1e-300))) +
  geom_point(aes(color = abs(avg_log2FC) > 1 & p_val_adj < 0.05), alpha = .5, size = 1) +
  scale_color_manual(values = c("grey70", "firebrick"), guide = "none") +
  ggrepel::geom_text_repel(
    data = head(de[order(-abs(de$avg_log2FC)), ], 12),
    aes(label = gene), size = 3, max.overlaps = 20) +
  labs(title = paste("Cluster", c_a, "vs cluster", c_b),
       subtitle = "Ojo: el eje Y NO es una medida válida de evidencia (ver sección 4)",
       x = "log2 Fold Change", y = "-log10(p ajustado)") +
  theme_bw()

# --------------------------------------------------------------------------
## 6. Comparar condiciones: pseudobulk ---------------------------------------
#
# Hasta aquí comparamos clusters entre sí. El otro tipo de pregunta —el que aparece en
# casi todos los proyectos reales— es:
#
# > *"Dentro de los monocitos, ¿qué genes cambian entre pacientes y controles?"*
#
# Aquí hay un segundo problema estadístico, **independiente del double dipping**.
#
## 6.1. Pseudorreplicación ---------------------------------------------------
#
# Si tienes 3 pacientes y 3 controles, y de cada uno recuperas 2.000 monocitos, un test
# célula a célula ve **6.000 vs 6.000 observaciones**. Pero las 2.000 células de un
# paciente **no son 2.000 réplicas independientes**: comparten genotipo, tratamiento,
# hora de extracción, lote de procesamiento. Son *pseudorréplicas*.
#
# El tamaño muestral real es **3 vs 3**.
#
# Al ignorarlo, los p-valores se vuelven astronómicamente pequeños y aparecen cientos de
# "genes diferenciales" que son variación entre individuos. La documentación oficial de
# Seurat muestra un ejemplo donde **1.649 genes** son significativos en el análisis
# single-cell y no lo son en pseudobulk.
#
## 6.2. La solución ----------------------------------------------------------
#
# **Sumar los conteos de todas las células del mismo tipo celular dentro de cada
# muestra** → una muestra, un perfil. Con eso se recupera un diseño de RNA-seq bulk
# clásico, y se aplican métodos con 15 años de validación: DESeq2, edgeR, limma-voom.
#
# `AggregateExpression()` hace la agregación.

# Nuestro dataset es de UN donante, así que no tiene réplicas biológicas.
# Simulamos un diseño para practicar la mecánica: 6 "muestras", 2 condiciones.
# ⚠️ Esto es un ejercicio de sintaxis, NO un experimento: las condiciones son
#    arbitrarias, así que NO debería salir nada diferencial. Ese es el punto.

set.seed(7)
pbmc$muestra   <- sample(paste0("S", 1:6), ncol(pbmc), replace = TRUE)
pbmc$condicion <- ifelse(pbmc$muestra %in% c("S1", "S2", "S3"), "control", "tratado")

table(pbmc$muestra, pbmc$condicion)[, ]

# Usamos identificadores de cluster que empiecen con letra ("C0", "C1", ...)
# para que los nombres de columna del pseudobulk sean fáciles de partir.
pbmc$cl_id <- paste0("C", as.character(pbmc$clusters))

# Agregación: SUMAMOS conteos por (cluster x muestra)
pseudo <- AggregateExpression(pbmc, assays = "RNA",
                              group.by = c("cl_id", "muestra"),
                              return.seurat = FALSE)$RNA

cat("Matriz pseudobulk:", nrow(pseudo), "genes x", ncol(pseudo), "perfiles\n")
cat("(", nlevels(pbmc$clusters), "clusters x 6 muestras )\n\n")
cat("Formato de los nombres de columna:\n"); print(head(colnames(pseudo)))
pseudo[1:4, 1:5]

# Separamos cluster y muestra a partir de los nombres de columna
partes <- do.call(rbind, strsplit(colnames(pseudo), "_", fixed = TRUE))
info <- data.frame(col = colnames(pseudo),
                   cluster = partes[, 1],
                   muestra = partes[, 2])
head(info)

# Tomamos un cluster y corremos DESeq2 a nivel de MUESTRA
cl   <- info$cluster[1]
sub  <- info[info$cluster == cl, ]
mat  <- round(as.matrix(pseudo[, sub$col]))
mat  <- mat[rowSums(mat) > 10, ]

coldata <- data.frame(
  row.names = sub$col,
  condicion = factor(ifelse(sub$muestra %in% c("S1", "S2", "S3"), "control", "tratado"),
                     levels = c("control", "tratado"))
)
print(coldata)

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
suppressPackageStartupMessages(library(DESeq2))

dds <- DESeqDataSetFromMatrix(mat, coldata, ~ condicion)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds)

cat("Cluster", cl, "— pseudobulk (n = 3 vs 3):\n")
cat("  Genes con padj < 0.05:", sum(res$padj < 0.05, na.rm = TRUE), "\n")

# El mismo contraste, pero célula a célula (lo que NO hay que hacer)
cl_orig <- sub("^C", "", cl)          # volvemos al identificador original
Idents(pbmc) <- "clusters"

sc_res <- FindMarkers(pbmc, subset.ident = cl_orig, group.by = "condicion",
                      ident.1 = "tratado", ident.2 = "control",
                      logfc.threshold = 0, verbose = FALSE)

cat("MISMO contraste, test célula a célula:\n")
cat("  Genes con p_val_adj < 0.05:", sum(sc_res$p_val_adj < 0.05), "\n\n")
cat("Las condiciones fueron asignadas AL AZAR: la respuesta correcta es 0.\n")
cat("Cada gen que aparezca aquí es un falso positivo.\n")

# > **La comparación de estos dos bloques es el resumen de la clase.** Mismo dato, misma
#   pregunta, mismo software. La diferencia está solo en qué se considera una unidad
#   experimental independiente.
# >
# > En un experimento real con efecto verdadero, el pseudobulk seguiría encontrando los
#   genes que realmente cambian — pero sin arrastrar cientos de falsos positivos. No es
#   un método más conservador porque sí: es el método que responde la pregunta que uno
#   cree estar haciendo.
#
## 6.3. Por qué esto importa: una simulación ---------------------------------
#
# El bloque anterior enseña la **mecánica** de `AggregateExpression` + DESeq2, pero no
# puede demostrar el problema: nuestro dataset es de un solo donante, así que al sortear
# las "muestras" no hay ninguna variabilidad entre ellas. Sin variabilidad entre
# individuos, el test célula a célula tampoco falla — y entonces la comparación no
# prueba nada.
#
# Simulemos el caso real: **6 donantes**, 3 control y 3 tratado, 400 células cada uno,
# y **ningún efecto de condición**. Lo único que varía es de donante a donante, como en
# cualquier experimento con personas. La respuesta correcta sigue siendo cero.

set.seed(11)
n_g <- 1500; n_m <- 6; n_c <- 400
cond_m <- rep(c("control", "tratado"), each = 3)     # la condición es de la MUESTRA

simular <- function(sigma) {
  tasa_base <- rlnorm(n_g, meanlog = 0.5, sdlog = 1)
  # variabilidad biológica entre donantes; ningún término de condición
  tasa <- tasa_base * exp(matrix(rnorm(n_g * n_m, 0, sigma), n_g, n_m))
  M <- matrix(rpois(n_g * n_m * n_c, lambda = tasa[, rep(seq_len(n_m), each = n_c)]),
              nrow = n_g)
  list(M = M, muestra = rep(seq_len(n_m), each = n_c))
}

# Wilcoxon vectorizado (aproximación normal con corrección por empates): mismo test
# que usa Seurat por defecto, pero para las 1.500 filas de una sola vez.
wilcox_filas <- function(M, g) {
  n1 <- sum(g); n2 <- sum(!g); n <- n1 + n2
  R <- t(apply(M, 1, rank))
  W <- rowSums(R[, g, drop = FALSE])
  emp <- apply(M, 1, function(x) { tt <- table(x); sum(tt^3 - tt) })
  sdW <- sqrt(n1 * n2 / 12 * ((n + 1) - emp / (n * (n - 1))))
  2 * pnorm(-abs((W - n1 * (n + 1) / 2) / sdW))
}

falsos <- function(sigma) {
  sim <- simular(sigma)
  M <- sim$M; muestra <- sim$muestra
  cond_cel <- cond_m[muestra]

  # (a) cada célula como una réplica
  tot <- colSums(M)
  Mn  <- log1p(t(t(M) / pmax(tot, 1) * 1e4))
  p_cel <- wilcox_filas(Mn, cond_cel == "control")

  # (b) pseudobulk: sumar los conteos de cada muestra -> 6 perfiles -> 3 vs 3
  pb    <- sapply(seq_len(n_m), function(k) rowSums(M[, muestra == k, drop = FALSE]))
  lcpm  <- log2(t(t(pb) / colSums(pb)) * 1e6 + 1)
  p_pb  <- apply(lcpm, 1, function(x)
                 t.test(x[cond_m == "control"], x[cond_m == "tratado"])$p.value)

  c(celula     = sum(p.adjust(p_cel, "BH") < 0.05),
    pseudobulk = sum(p.adjust(p_pb,  "BH") < 0.05))
}

cat("Falsos positivos (p ajustado < 0.05) — la respuesta correcta es 0:\n\n")
tabla <- t(sapply(c(0, 0.1, 0.3), falsos))
rownames(tabla) <- paste("sigma =", c(0, 0.1, 0.3))
print(tabla)

# > Leé la tabla de arriba hacia abajo.
# >
# > Con `sigma = 0` —donantes idénticos— el test célula a célula **no falla**. O sea:
#   el problema no es tener muchas células.
# >
# > En cuanto los donantes empiezan a diferir entre sí, que es lo que pasa siempre, el
#   test célula a célula declara diferencial a buena parte del transcriptoma. El
#   pseudobulk no encuentra nada, en las tres filas, que es lo correcto.
# >
# > **El mecanismo:** las células de un mismo donante no son independientes. Al
#   ignorarlo, la variabilidad ENTRE INDIVIDUOS entra en la cuenta como si fuera
#   variabilidad ENTRE CÉLULAS, y el test la lee como señal.
#
# > Ojo con la confusión frecuente: la pseudorreplicación y el *double dipping* de la
#   sección 4 son problemas **independientes**. Uno es usar los datos dos veces; el otro
#   es contar mal las unidades experimentales. Se pueden tener los dos a la vez.

### Requisitos para pseudobulk -----------------------------------------------
#
# - **Al menos 3 réplicas biológicas por condición** (individuos distintos, no pocillos
#   distintos).
# - Suficientes células por muestra y tipo celular (~20 como mínimo grosero; con menos,
#   el perfil agregado es inestable).
# - Cuidado con la composición: si una condición tiene muchas más células de un tipo,
#   eso es en sí mismo un resultado (abundancia diferencial, que se analiza aparte con
#   métodos como `miloR` o `scCODA`).

saveRDS(pbmc, ck2)
tam_archivo(ck2)

# --------------------------------------------------------------------------
## 7. Síntesis ---------------------------------------------------------------
#
# 1. El clustering por grafos (KNN → SNN → modularidad) no supone forma ni tamaño de los
#    clusters. **Leiden** es preferible a Louvain y desde Seurat 5.3.1 corre en R puro.
# 2. **La resolución define la pregunta.** Se elige con `clustree` y criterios de
#    estabilidad, y se **reporta**.
# 3. **El clustering siempre devuelve clusters**, incluso sobre ruido puro.
# 4. **Los p-valores post-clustering no son p-valores.** Sirven para rankear candidatos.
#    Para inferencia hay `ClusterDE` y *count splitting*.
# 5. **Para comparar condiciones: pseudobulk.** Las células de una muestra no son
#    réplicas independientes.
