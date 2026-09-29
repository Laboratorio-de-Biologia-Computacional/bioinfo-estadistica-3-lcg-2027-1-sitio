# =============================================================================
# Bioinformática y Estadística 3 — Módulo scRNA-seq + CITE-seq
# Clase 3 — Anotación celular
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
## Módulo scRNA-seq + CITE-seq · Clase 3 de 4 --------------------------------
### Anotación celular: marcadores, referencias y ontologías ------------------
#
# **Docente:** Dr. Danilo Ceschin · **Fecha:** martes 22 de septiembre de 2026 ·
# **Duración:** 2 h
#
# --------------------------------------------------------------------------
#
# Tenemos clusters numerados. Hoy les ponemos nombre.
#
# Parece el paso trivial del pipeline y es el más difícil de todos, por una razón
# conceptual: decir *"el cluster 4 son monocitos CD14+"* **no es un resultado, es una
# hipótesis**. Es una afirmación sobre identidad biológica derivada de una medición de
# ARNm, mediada por decisiones de preprocesamiento que tomamos nosotros en las clases 1
# y 2.
#
# Todo el resto del análisis va a apoyarse en esos nombres. Si están mal, todo lo que
# sigue está mal, y —a diferencia de un error de código— nada va a fallar ruidosamente.
#
# --------------------------------------------------------------------------
#
### Objetivos de la clase ----------------------------------------------------
#
# 1. Usar marcadores canónicos y conocer sus límites en scRNA-seq.
# 2. Anotar automáticamente con **SingleR + celldex** (API nueva `fetchReference`), por
#    célula y por cluster.
# 3. Mapear contra una referencia curada (**reference mapping**) y comparar los tres
#    métodos.
# 4. Diagnosticar **RNA ambiental** y otros confundidores que producen anotaciones
#    falsas.
# 5. Reportar la anotación con un **nivel de confianza explícito**.
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

# --------------------------------------------------------------------------
## 1. Retomamos el objeto clusterizado ---------------------------------------

if (file.exists(ck2)) {
  pbmc <- readRDS(ck2)
  DefaultAssay(pbmc) <- "RNA"; Idents(pbmc) <- "clusters"
  cat("Cargado:", ncol(pbmc), "células,", nlevels(pbmc$clusters), "clusters\n")
} else {
  cat("No encontrado. Ejecuta el bloque de regeneración de abajo.\n")
}

# Regeneración desde cero si hace falta (~5 min).
# Reproduce los pasos de las clases 1 y 2 que importan para hoy, con los MISMOS
# parámetros: 2000 HVGs, 30 PCs, k = 20 y clustering con Leiden (igraph) a
# resolución 0.5 — el algoritmo que elegimos en la clase 2 porque corrige los
# clusters internamente desconectados de Louvain. Si aquí usáramos otro, quien
# regenera obtendría clusters distintos de quien cargó el checkpoint, y las
# etiquetas de la sección 7 no coincidirían.
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
  io <- function(x, k = 5, hi = FALSE) { m <- median(x); d <- mad(x)
    if (hi) x > m + k * d else x < m - k * d | x > m + k * d }
  pbmc <- pbmc[, !(io(log1p(pbmc$nCount_RNA)) | io(log1p(pbmc$nFeature_RNA)) |
                   io(pbmc$percent.mt, 3, TRUE) | pbmc$percent.mt > 20)]
  sce <- scDblFinder(as.SingleCellExperiment(pbmc), verbose = FALSE)
  pbmc$doublet_score <- colData(sce)$scDblFinder.score   # lo usa el diagnóstico de §5.2
  pbmc <- pbmc[, colData(sce)$scDblFinder.class == "singlet"]
  pbmc <- NormalizeData(pbmc, verbose = FALSE) |>
          FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
          ScaleData(verbose = FALSE) |>
          RunPCA(npcs = 50, verbose = FALSE) |>
          RunUMAP(dims = 1:30, verbose = FALSE) |>
          FindNeighbors(dims = 1:30, k.param = 20, verbose = FALSE) |>
          FindClusters(resolution = 0.5, algorithm = 4, leiden_method = "igraph",
                       verbose = FALSE)
  pbmc$clusters <- Idents(pbmc)
  cat("Regenerado:", ncol(pbmc), "células,", nlevels(pbmc$clusters), "clusters\n")
}
DimPlot(pbmc, label = TRUE) + NoLegend() + ggtitle("Punto de partida")

# --------------------------------------------------------------------------
## 2. Marcadores canónicos ---------------------------------------------------
#
### 2.1. De dónde vienen -----------------------------------------------------
#
# Los marcadores clásicos de inmunología (CD3, CD4, CD8, CD19, CD14...) vienen de
# **citometría de flujo**: son proteínas de superficie, elegidas porque son detectables
# con anticuerpos y porque discriminan bien.
#
# Cuando los trasladamos a scRNA-seq medimos otra cosa: **el ARN mensajero que codifica
# esa proteína**. Y ahí aparecen tres problemas:
#
# | Problema | Qué pasa |
# |---|---|
# | **Dropout** | Un gen expresado puede tener 0 conteos por azar de muestreo. `CD4` es notoriamente difícil de detectar en scRNA-seq, aunque la proteína esté ahí. |
# | **ARNm ≠ proteína** | La correlación entre ambos es moderada. `CD8A` puede transcribirse sin que la proteína esté en superficie, y viceversa. |
# | **Expresión compartida** | Casi ningún gen es exclusivo de un tipo celular. `NKG7` marca NK **y** T CD8 citotóxicos. |
#
# Este es, dicho sea de paso, uno de los argumentos más fuertes a favor de CITE-seq —
# que vemos el jueves.
#
### 2.2. Panel canónico de PBMC ----------------------------------------------

marcadores_canonicos <- list(
  "T (pan)"          = c("CD3D", "CD3E", "TRAC"),
  "T CD4"            = c("IL7R", "CD4", "CCR7"),
  "T CD8"            = c("CD8A", "CD8B", "GZMK"),
  "NK"               = c("GNLY", "NKG7", "KLRD1", "NCAM1"),
  "B"                = c("MS4A1", "CD79A", "CD79B"),
  "Plasmablasto"     = c("MZB1", "JCHAIN", "XBP1"),
  "Monocito CD14"    = c("CD14", "LYZ", "S100A8", "S100A9"),
  "Monocito CD16"    = c("FCGR3A", "MS4A7", "CDKN1C"),
  "DC convencional"  = c("FCER1A", "CST3", "CLEC10A"),
  "DC plasmacitoide" = c("LILRA4", "IL3RA", "GZMB"),
  "Plaqueta"         = c("PPBP", "PF4"),
  "Proliferación"    = c("MKI67", "TOP2A")
)

# Nos quedamos solo con los genes presentes en el objeto
marcadores_canonicos <- lapply(marcadores_canonicos, function(g) g[g %in% rownames(pbmc)])
str(marcadores_canonicos)


DotPlot(pbmc, features = marcadores_canonicos, cluster.idents = TRUE) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8),
        strip.text  = element_text(size = 7, angle = 60)) +
  labs(title = "Marcadores canónicos por cluster",
       subtitle = "Tamaño = % de células que lo expresan · Color = expresión media escalada")

# **Cómo se lee este gráfico.** Busca, para cada cluster (fila), qué bloque de
# marcadores se enciende. Un cluster con `CD3D`+`CD8A`+`GZMK` altos es T CD8. Un cluster
# con `CD14`+`LYZ`+`S100A8` es monocito clásico.
#
### 2.3. ¿Tipo celular o estado celular? -------------------------------------
#
# Antes de ponerle nombre a un cluster hay que descartar que lo que lo define no sea un
# **estado** en vez de un **tipo**. Un tipo celular es una identidad estable: un monocito
# no se convierte en linfocito. Un estado es transitorio: una célula estresada, activada
# o en división sigue siendo del mismo tipo.
#
# La señal de alarma es concreta: si lo que separa a un cluster son genes de **choque
# térmico**, de **respuesta temprana** o de **ciclo celular**, lo más probable es que el
# clustering esté capturando cómo se procesó la muestra, no qué célula es.
#
# Los genes de respuesta temprana (`FOS`, `JUN`, `EGR1`) son especialmente traicioneros:
# se inducen en minutos durante la disociación, así que su expresión mide **el
# protocolo**, no la biología del donante.

programas_estado <- list(
  "Estrés / choque térmico" = c("HSPA1A", "HSPA1B", "HSPB1", "DNAJB1"),
  "Respuesta temprana"      = c("FOS", "FOSB", "JUN", "JUNB", "EGR1"),
  "Ciclo celular"           = c("MKI67", "TOP2A", "PCNA", "CCNB1")
)
programas_estado <- lapply(programas_estado, function(g) g[g %in% rownames(pbmc)])


DotPlot(pbmc, features = programas_estado, cluster.idents = TRUE) +
  RotatedAxis() +
  labs(title = "Programas de ESTADO (no de tipo celular)",
       subtitle = "Si un cluster se enciende acá y no en el panel canónico, sospecha")

# Un score por célula resume cada programa mejor que gen por gen: promedia el bloque y
# lo corrige contra genes de expresión comparable, así no confunde "programa activo" con
# "célula con más conteos".
pbmc <- AddModuleScore(pbmc, features = programas_estado["Estrés / choque térmico"],
                       name = "score_estres", seed = 42)

VlnPlot(pbmc, "score_estres1", group.by = "clusters", pt.size = 0) +
  NoLegend() + labs(title = "Score de estrés por cluster", x = "cluster")

# > **Cómo se interpreta.** En este dataset —PBMC de donante sano, procesadas en fresco—
#   ningún cluster debería destacarse. Si mañana analizas un tejido sólido disociado en
#   caliente, este mismo gráfico va a mostrar un cluster entero encendido: ese cluster no
#   es un tipo celular, es el daño del protocolo. La decisión ahí no es anotarlo, es
#   volver al diseño experimental.
#
# > **En nuestros datos**, el DotPlot de estado muestra dos cosas distintas. Los genes de
#   respuesta temprana (`FOS`, `JUNB`) aparecen en TODOS los clusters, algo más en los
#   mieloides: es el piso que deja cualquier protocolo, y como es parejo no separa a
#   nadie. El ciclo celular, en cambio, se enciende en un solo cluster, el de los
#   plasmablastos. Ahí el estado (proliferación) es biología real superpuesta a un tipo
#   celular: los plasmablastos son células que se están dividiendo. Estado no siempre
#   quiere decir artefacto; quiere decir que no alcanza para definir un tipo.
#
# Fíjate también en lo que **no** es limpio: es habitual que un cluster tenga señal
# intermedia de dos bloques. Eso puede ser (a) un estado de transición real, (b)
# doublets remanentes, o (c) RNA ambiental. Volvemos a esto en la sección 5.

FeaturePlot(pbmc, ncol = 4, order = TRUE,
            features = c("CD3E", "IL7R", "CD8A", "NKG7",
                         "MS4A1", "CD14", "FCGR3A", "FCER1A",
                         "LILRA4", "PPBP", "MZB1", "MKI67")) &
  theme(plot.title = element_text(size = 11)) & NoLegend()

# --------------------------------------------------------------------------
## 3. Anotación automática: SingleR + celldex --------------------------------
#
# `SingleR` correlaciona el perfil de cada célula con perfiles de referencia de tipos
# celulares purificados (bulk RNA-seq), usando solo los genes que mejor discriminan cada
# par de tipos.
#
### ⚠️ Cambio de API respecto del material de 2024 ---------------------------
#
# La forma antigua, que todavía circula en tutoriales y material de cursos previos:
#
# ```r
# ref <- celldex::BlueprintEncodeData()          # ← forma antigua
# ```
#
# `celldex` migró de `ExperimentHub` a **gypsum**, y ahora las referencias se piden
# **con nombre y versión**:
#
# ```r
# ref <- celldex::fetchReference("blueprint_encode", "2024-02-26")
# ```
#
# El cambio no es cosmético: las referencias ahora están **versionadas** (puedes
# reproducir exactamente un análisis de hace dos años) y traen las etiquetas mapeadas a
# **Cell Ontology**, lo que las hace interoperables entre estudios.
#
# | Nombre | Contenido |
# |---|---|
# | `blueprint_encode` | Células inmunes y estromales humanas (Blueprint + ENCODE) |
# | `hpca` | Human Primary Cell Atlas — cobertura amplia de tejidos |
# | `monaco_immune` | Poblaciones inmunes humanas, la más fina para PBMC |
# | `dice` | Database of Immune Cell Expression |
# | `novershtern_hematopoietic` | Hematopoyesis humana, poblaciones progenitoras |
# | `immgen`, `mouse_rnaseq` | Ratón |

suppressPackageStartupMessages({ library(SingleR); library(celldex); library(SingleCellExperiment) })

# API NUEVA. La primera vez descarga y cachea (~1-2 min).
ref_bp    <- celldex::fetchReference("blueprint_encode", "2024-02-26")
ref_monaco <- celldex::fetchReference("monaco_immune",   "2024-02-26")

cat("Blueprint/ENCODE:", ncol(ref_bp), "perfiles,",
    length(unique(ref_bp$label.main)), "etiquetas principales\n")
cat("Monaco:", ncol(ref_monaco), "perfiles,",
    length(unique(ref_monaco$label.fine)), "etiquetas finas\n\n")
print(sort(unique(ref_bp$label.main)))

### 3.1. Por célula vs por cluster -------------------------------------------
#
# SingleR puede trabajar de dos maneras, y la diferencia importa:
#
# - **Por célula** (`clusters = NULL`): anota cada célula de forma independiente.
#   Ventaja: detecta poblaciones que el clustering fusionó. Desventaja: más ruidoso,
#   porque cada célula individual tiene poca información.
# - **Por cluster** (`clusters = pbmc$clusters`): promedia el perfil del cluster y anota
#   el promedio. Ventaja: mucho más estable. Desventaja: hereda todos los errores del
#   clustering.
#
# **Haz las dos.** Donde coinciden, tienes confianza alta. Donde no, hay algo que mirar.

# La matriz dispersa se pasa DIRECTAMENTE.
# ⚠️ NO usar as.matrix(): en 7000 células son varios GB y R se queda sin memoria.
expr <- GetAssayData(pbmc, assay = "RNA", layer = "data")
cat("Clase de la matriz:", class(expr)[1], "— dispersa, como debe ser\n")

# ---------------------------------------------------------------------------
# Un detalle de implementación que conviene conocer
#
# Cuando le pasas `clusters =`, SingleR no hace nada mágico: suma los perfiles de
# todas las células de cada cluster y anota ese perfil agregado. Esa suma la
# delega en el paquete `scrapper`, que SingleR declara en Suggests y por lo tanto
# NO se instala automáticamente con él. Si falta, el error es:
#
#     Error in loadNamespace(x) : there is no package called 'scrapper'
#
# La sección 0 ya lo instala. La función de abajo agrega además un plan B: si en
# alguna computadora `scrapper` no compila, hace la agregación a mano y la clase
# sigue igual. De paso deja a la vista qué significa realmente `clusters =`.
# ---------------------------------------------------------------------------
anotar_por_cluster <- function(expr, ref, labels, clusters) {
  if (requireNamespace("scrapper", quietly = TRUE)) {
    return(SingleR(test = expr, ref = ref, labels = labels, clusters = clusters))
  }
  message("scrapper no disponible: se agregan los clusters a mano.")
  cl  <- factor(clusters)
  ind <- Matrix::sparse.model.matrix(~ 0 + cl)   # células x clusters
  colnames(ind) <- levels(cl)
  SingleR(test = expr %*% ind, ref = ref, labels = labels)
}

# (a) por célula
anot_celula <- SingleR(test = expr, ref = ref_bp, labels = ref_bp$label.main)

# (b) por cluster
anot_cluster <- anotar_por_cluster(expr, ref_bp, ref_bp$label.main, pbmc$clusters)

pbmc$singler_celula <- anot_celula$pruned.labels
table(pbmc$singler_celula, useNA = "ifany")

# Etiqueta por cluster, propagada a las células
#
# ⚠️ Trampa de R que vale la pena conocer: al indexar un vector con nombres, el
# resultado hereda los nombres del ÍNDICE, no los del vector original. Acá
# `mapa[as.character(pbmc$clusters)]` devuelve los valores correctos pero
# nombrados "0", "1", "2"... en lugar de códigos de barras. Seurat, al recibir
# un vector con nombres, intenta emparejarlos con sus células, no encuentra
# ninguno y aborta con:
#
#     Error: No cell overlap between new meta data and Seurat object
#
# La solución es nombrar el vector por célula antes de asignarlo. `unname()`
# también funciona (Seurat cae en emparejamiento posicional), pero nombrar
# explícitamente deja claro que el orden es el de las células del objeto.
mapa <- setNames(anot_cluster$pruned.labels, rownames(anot_cluster))
pbmc$singler_cluster <- setNames(mapa[as.character(pbmc$clusters)], colnames(pbmc))

data.frame(cluster         = rownames(anot_cluster),
           etiqueta        = anot_cluster$labels,
           etiqueta_pruned = anot_cluster$pruned.labels,
           n_celulas       = as.integer(table(pbmc$clusters)[rownames(anot_cluster)]))

(DimPlot(pbmc, group.by = "singler_celula", label = TRUE, repel = TRUE, label.size = 3) +
   ggtitle("SingleR por célula") + NoLegend()) |
(DimPlot(pbmc, group.by = "singler_cluster", label = TRUE, repel = TRUE, label.size = 3) +
   ggtitle("SingleR por cluster") + NoLegend())

### 3.2. Diagnóstico: no confiar en la etiqueta, mirar el score --------------
#
# SingleR **siempre** devuelve una etiqueta: asigna la más parecida, incluso si el
# parecido es malo. Hay cuatro diagnósticos que hay que mirar sí o sí:
#
# 1. **`plotScoreHeatmap`** — los scores de todas las etiquetas para cada célula. Una
#    célula bien asignada tiene un score claramente alto para una etiqueta. Una célula
#    ambigua tiene scores parecidos para varias.
# 2. **`plotDeltaDistribution`** — el delta es la diferencia entre el mejor score y la
#    mediana de los demás. Delta chico = asignación poco confiable.
# 3. **`pruned.labels`** — SingleR pone `NA` a las asignaciones de baja confianza.
#    **Esos `NA` son información, no un error a arreglar.**
# 4. **La tabla cluster × etiqueta** — un cluster coherente concentra sus células en
#    una sola etiqueta. Uno repartido entre dos es una alerta.

plotScoreHeatmap(anot_celula, show.pruned = TRUE, max.labels = 20)

plotDeltaDistribution(anot_celula, ncol = 7, size = 0.4)

n_na <- sum(is.na(anot_celula$pruned.labels))
cat("Células descartadas por baja confianza (pruned):", n_na,
    sprintf("(%.1f %%)\n", 100 * n_na / ncol(pbmc)))

# ¿Dónde están esas células en el UMAP? Casi nunca están repartidas al azar.
pbmc$sin_asignar <- is.na(pbmc$singler_celula)
DimPlot(pbmc, group.by = "sin_asignar", cols = c("grey85", "firebrick"), order = TRUE) +
  ggtitle("Células que SingleR no pudo asignar con confianza")

# Tabla de contingencia: cluster x etiqueta. Es EL diagnóstico de coherencia.
tab <- table(cluster = pbmc$clusters, etiqueta = pbmc$singler_celula)

suppressPackageStartupMessages(library(pheatmap))
pheatmap(log1p(tab), cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = tab, number_format = "%.0f", fontsize_number = 7,
         main = "Cluster (filas) vs etiqueta SingleR por célula (columnas)")

# > **Cómo se lee.** Un cluster **coherente** concentra casi todas sus células en una
#   sola columna. Un cluster repartido entre dos etiquetas es una alerta: puede ser
#   sobre-clusterización, una población de transición, doublets, o que la referencia no
#   contiene ese tipo celular.
# >
# > Cuando la referencia no tiene el tipo celular que estás mirando, SingleR igual
#   asigna la etiqueta más cercana. Ese es el modo de falla más peligroso, porque el
#   resultado se ve perfectamente razonable.

# Referencia más fina para las poblaciones T, que es donde Blueprint se queda corto
anot_fina <- anotar_por_cluster(expr, ref_monaco, ref_monaco$label.fine, pbmc$clusters)

cl_ord <- rownames(anot_fina)
data.frame(cluster        = cl_ord,
           blueprint_main = anot_cluster[cl_ord, "pruned.labels"],
           monaco_fine    = anot_fina[cl_ord, "pruned.labels"])

# --------------------------------------------------------------------------
## 4. Reference mapping ------------------------------------------------------
#
# Un tercer enfoque, distinto de los dos anteriores: en lugar de correlacionar contra
# perfiles bulk, **proyectar nuestras células dentro del espacio de una referencia
# single-cell ya anotada** y transferir las etiquetas de los vecinos.
#
# Es la lógica de **Azimuth** (referencia PBMC multimodal del laboratorio Satija, con
# más de 160.000 células anotadas en tres niveles jerárquicos). Por debajo son las
# funciones `FindTransferAnchors` + `TransferData` de Seurat — exactamente las mismas
# que se usan para transferir etiquetas entre scATAC-seq y scRNA-seq. Conviene ver la
# mecánica hoy.
#
# **Ventajas:** la referencia es single-cell (no bulk), tiene jerarquía de niveles, y
# devuelve un score de predicción por célula.
# **Desventajas:** solo sirve si existe una referencia del tejido que estás estudiando,
# y hereda todos los sesgos de esa referencia.

# Azimuth descarga una referencia grande. Si el tiempo aprieta, salta a la sección 5.
instalar_azimuth <- FALSE   # ← pon TRUE para ejecutarlo

if (instalar_azimuth) {
  if (!requireNamespace("Azimuth", quietly = TRUE))
    remotes::install_github("satijalab/azimuth", upgrade = "never", quiet = TRUE)
  suppressPackageStartupMessages(library(Azimuth))
  pbmc <- RunAzimuth(pbmc, reference = "pbmcref")
  head(pbmc@meta.data[, grep("predicted", colnames(pbmc@meta.data))])
} else {
  cat("Azimuth desactivado. Abajo hacemos la misma idea 'a mano' con FindTransferAnchors.\n")
}

### 4.1. La mecánica, a mano -------------------------------------------------
#
# Para ver qué hace Azimuth por debajo, hagamos transferencia de etiquetas usando como
# "referencia" una mitad de los datos, con las etiquetas **por célula** de SingleR
# (§3.1), y prediciendo la otra mitad. Es un ejercicio artificial, pero muestra las
# funciones que importan y permite **medir** qué tan bien funciona el método, porque
# conocemos la respuesta.
#
# ¿Por qué las etiquetas por célula y no las por cluster? Las etiquetas por cluster son,
# por construcción, una función del grafo de vecinos: transferirlas es casi trivial
# (concordancia ~99 %, score 1 en todas las células) y el ejercicio no enseña nada. Las
# etiquetas por célula tienen la incertidumbre real de la anotación, y el score de
# predicción la hace visible.

set.seed(11)
idx <- sample(colnames(pbmc), floor(ncol(pbmc) / 2))
ref_obj   <- subset(pbmc, cells = idx)
query_obj <- subset(pbmc, cells = setdiff(colnames(pbmc), idx))

ref_obj <- ref_obj |> NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(verbose = FALSE) |> ScaleData(verbose = FALSE) |>
  RunPCA(npcs = 30, verbose = FALSE)
query_obj <- NormalizeData(query_obj, verbose = FALSE)

anchors <- FindTransferAnchors(reference = ref_obj, query = query_obj,
                               dims = 1:30, reduction = "pcaproject", verbose = FALSE)

# Las etiquetas de referencia no pueden tener NA: TransferData las trataría como
# una categoría inválida. Los clusters que SingleR dejó sin asignar (pruned) pasan
# a ser una etiqueta explícita, "Sin asignar", que también se transfiere.
# (El aviso "Different features in new layer data..." que aparece al preparar
# la referencia es inofensivo: ScaleData recalcula sobre los HVGs de esta mitad.)
etiq_ref <- as.character(ref_obj$singler_celula)
etiq_ref[is.na(etiq_ref)] <- "Sin asignar"

pred <- TransferData(anchorset = anchors, refdata = etiq_ref,
                     dims = 1:30, verbose = FALSE)

query_obj <- AddMetaData(query_obj, pred)
head(pred[, 1:3])

# ⚠️ Cómo se llaman las columnas del resultado: depende de cómo llamaste a
# `TransferData`, y es una fuente clásica de confusión.
#
#   TransferData(anchorset, refdata = <vector>)          ← lo que hicimos acá
#       devuelve un data.frame con:
#         predicted.id                 etiqueta predicha
#         prediction.score.<etiqueta>  una columna por etiqueta posible
#         prediction.score.max         el score de la etiqueta ganadora
#
#   TransferData(anchorset, refdata = list(id = <vector>), query = obj)
#       devuelve el OBJETO, con la metadata predicted.id y predicted.id.score.
#       Es la convención que usan Azimuth y MapQuery, donde vas a ver nombres
#       como `predicted.celltype.l2.score`.
#
# Los dos caminos dan lo mismo; solo cambia el nombre de la columna del score.
cat("Columnas devueltas:", paste(colnames(pred), collapse = ", "), "\n\n")

etiq_query <- as.character(query_obj$singler_celula)
etiq_query[is.na(etiq_query)] <- "Sin asignar"
concordancia <- mean(query_obj$predicted.id == etiq_query)
cat("Concordancia con la anotación original:", round(100 * concordancia, 1), "%\n\n")

p1 <- ggplot(query_obj@meta.data, aes(prediction.score.max)) +
  geom_histogram(bins = 40, fill = "steelblue") +
  labs(title = "Score de predicción por célula", x = "score") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- ggplot(query_obj@meta.data, aes(predicted.id, prediction.score.max)) +
  geom_boxplot(outlier.size = .3) + coord_flip() +
  labs(title = "Score por etiqueta predicha", x = NULL) + theme_bw()
p1 | p2

# ¿Dónde están las células con score bajo? La mitad "query" conserva el UMAP
# del objeto completo, así que se puede pintar directamente.
FeaturePlot(query_obj, "prediction.score.max", reduction = "umap", pt.size = 0.4) +
  scale_colour_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                         oob = scales::squish) +
  labs(title = "Score de predicción sobre el UMAP",
       subtitle = "Claro = asignación poco confiable · oscuro = score cercano a 1")

# > **Lo importante no es la concordancia global**, sino que el score identifica *dónde*
#   el método es confiable. Las etiquetas con score bajo y disperso son justamente las
#   poblaciones que hay que revisar a mano.
#
### 4.2. El panorama completo de métodos (2026) ------------------------------
#
# | Familia | Ejemplos | Qué asume | Cuándo usarlo |
# |---|---|---|---|
# | Correlación con bulk | **SingleR** + celldex | La referencia contiene el tipo celular | Primera pasada, siempre |
# | Reference mapping | **Azimuth**, `TransferData` | Existe referencia single-cell del tejido | Tejidos bien caracterizados (PBMC, médula, pulmón) |
# | Clasificador entrenado | **CellTypist** (Python) | El modelo vio tipos similares | Inmunología; modelos preentrenados de bajo y alto nivel |
# | Consenso | **popV** | Varios métodos, se mide el acuerdo | Cuando la confianza importa más que la velocidad |
# | Modelos fundacionales | **scGPT**, **Geneformer** | Preentrenamiento a gran escala transfiere | Ver abajo ⚠️ |
#
# **Sobre los modelos fundacionales.** Es el tema de moda y hay que tener una opinión
# informada. Los benchmarks independientes publicados en 2025–2026 (por ejemplo el
# framework BioLLM, en *Patterns*) muestran de forma consistente que en anotación de
# tipo celular **zero-shot** estos modelos **no superan de manera confiable** a métodos
# clásicos que son órdenes de magnitud más baratos. Con fine-tuning sobre datos del
# dominio mejoran, pero entonces ya no son "zero-shot" y compiten en igualdad de
# condiciones con un clasificador entrenado.
#
# Conclusión para la práctica: es un área de investigación muy activa y vale la pena
# seguirla, pero en 2026 **no reemplaza** a SingleR ni a reference mapping para el
# trabajo de todos los días.

# --------------------------------------------------------------------------
## 5. Confundidores que producen anotaciones falsas --------------------------
#
# Tres artefactos hacen que una anotación se vea perfecta y sea incorrecta.
#
### 5.1. RNA ambiental (la "sopa") -------------------------------------------
#
# Cuando las células se lisan durante la preparación, su ARNm queda libre en la
# suspensión y **se reparte en todas las gotas**. Cada célula recibe entonces un poco
# del transcriptoma de sus vecinas muertas.
#
# Consecuencia directa sobre la anotación: en un PBMC, los monocitos expresan `LYZ` a
# niveles altísimos. La sopa lleva `LYZ` a todas las gotas, y de pronto los linfocitos T
# "expresan" `LYZ`. Uno concluye que hay una población T atípica. No la hay.
#
# Diagnóstico rápido, sin instalar nada: **buscar marcadores exclusivos de un tipo
# celular en los clusters donde no deberían estar.**

genes_sopa <- c("LYZ", "S100A8", "HBB", "PPBP")
genes_sopa <- genes_sopa[genes_sopa %in% rownames(pbmc)]

# Cuantificación: % de células de cada tipo con conteo > 0 de cada gen.
# Un violín no sirve para esto: fuera de su tipo celular estos genes valen 0 en
# casi todas las células y el violín se ve como una línea. Lo que importa es
# CUÁNTAS células tienen al menos un conteo donde no deberían.
sopa_tab <- sapply(genes_sopa, function(g) {
  ct <- GetAssayData(pbmc, layer = "counts")[g, ]
  round(100 * tapply(ct > 0, pbmc$singler_cluster, mean), 1)
})
sopa_tab

sopa_df <- as.data.frame(as.table(sopa_tab))
colnames(sopa_df) <- c("tipo", "gen", "pct")

ggplot(sopa_df, aes(gen, tipo, fill = pct)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.1f %%", pct)), size = 4) +
  scale_fill_gradient(low = "#F4F8F8", high = "#0F7173", limits = c(0, 100)) +
  labs(title = "% de células con al menos un conteo del gen",
       subtitle = "Filas: etiqueta de SingleR por cluster. La sopa sube el piso de TODAS las filas por igual",
       x = NULL, y = NULL, fill = "%") +
  theme_minimal(base_size = 13) + theme(panel.grid = element_blank())

# > **Cómo leerlo.** Mira la columna de `LYZ` en las filas que NO son monocitos. Ese
#   porcentaje es el piso que pone la sopa: en datasets muy contaminados supera el 30 %.
#   En este dataset ronda el 6–10 %, parejo en todos los linfocitos. Hay sopa, pero
#   poca: PBMC frescas, pocas células lisadas.
# >
# > Compáralo con `PPBP` (plaquetas): casi 0 % en linfocitos y ~3 % en monocitos. Ese
#   patrón **no** es sopa, porque la sopa llega a todas las gotas por igual. Una señal
#   restringida a un solo tipo celular sugiere otra cosa: aquí, agregados
#   plaqueta–monocito, que son biológicamente reales en sangre. Mismo gen "fuera de
#   lugar", dos explicaciones distintas: por eso se mira el patrón, no un número suelto.
# >
# > **Herramientas de corrección:** `SoupX` (estima la fracción de contaminación a
#   partir de las gotas vacías), `DecontX` (modelo de mezcla, no necesita gotas vacías)
#   y `CellBender` (red neuronal, Python, el más agresivo).
# >
# > **Advertencia importante**, respaldada por los benchmarks de descontaminación
#   publicados en 2026: la corrección **no es gratis**. Algunos métodos generan conteos
#   artificiales —agregan señal donde no había— lo cual es peor que el problema
#   original. Recomendación práctica: **diagnosticar siempre, corregir solo si el
#   diagnóstico lo justifica, y comparar los resultados con y sin corrección.**
#
### 5.2. Doublets remanentes -------------------------------------------------
#
# `scDblFinder` no atrapa todo, sobre todo los doublets homotípicos (dos células del
# mismo tipo, que son casi indetectables). Un doublet T+monocito aparece como una célula
# que expresa `CD3E` **y** `CD14`, y es fácil bautizarla como "población híbrida novel".
#
# Regla: antes de reportar una población que coexpresa marcadores de dos linajes,
# **descarta que sean doublets**.

# ¿Hay clusters enriquecidos en score de doublet o con perfil intermedio?
if ("doublet_score" %in% colnames(pbmc@meta.data)) {
  print(VlnPlot(pbmc, "doublet_score", group.by = "clusters", pt.size = 0) + NoLegend())
}

# Coexpresión de linajes mutuamente excluyentes
co <- FetchData(pbmc, vars = c("CD3E", "CD14", "MS4A1", "clusters"), layer = "data")
co$doble_T_mono <- co$CD3E > 0.5 & co$CD14 > 0.5
cat("Células que coexpresan CD3E y CD14:", sum(co$doble_T_mono),
    sprintf("(%.2f %%)\n", 100 * mean(co$doble_T_mono)))
round(100 * tapply(co$doble_T_mono, co$clusters, mean), 2)

# > **Cómo leerlo.** En este dataset coexpresa menos del 1 % de las células, y el
#   porcentaje más alto está en el cluster de monocitos CD14+. Es lo esperable de
#   unos pocos doublets T+monocito remanentes (son los dos tipos más abundantes, así
#   que son la combinación más probable) y no alcanza para formar un cluster propio.
#   Si un cluster ENTERO coexpresara marcadores de dos linajes, la primera hipótesis
#   sería doublets, no una población nueva.

### 5.3. Efectos de lote -----------------------------------------------------
#
# Si el dataset tiene varias muestras y no se integró, los clusters pueden separar por
# **muestra** en vez de por tipo celular. Anotar en esa condición produce "tipos
# celulares" que son en realidad lotes.
#
# Nuestro dataset es de un solo donante, así que no aplica. Pero es exactamente el tema
# de la corrección de efectos de lote, y reaparece en el análisis multimodal con Bridge
# Integration.
#
# **La verificación mínima siempre es la misma:** grafica el UMAP coloreado por
# muestra/lote antes de anotar nada.

# --------------------------------------------------------------------------
## 6. Ponerle el nombre correcto: ontologías ---------------------------------
#
# `"Mono"`, `"monocyte"`, `"CD14+ Mono"`, `"Monocito clásico"` — si cada laboratorio usa
# su vocabulario, los datos no se pueden comparar entre estudios.
#
# La **[Cell Ontology (CL)](https://obophenotype.github.io/cell-ontology/)** es un
# vocabulario controlado, jerárquico y con identificadores estables:
#
# - `CL:0000860` — *classical monocyte*
# - `CL:0000624` — *CD4-positive, alpha-beta T cell*
#
# Las referencias de `celldex` ya traen este mapeo, y **CZ CELLxGENE** exige anotaciones
# en CL para cualquier dataset que se deposite. Usar CL desde el principio hace que tu
# anotación sea reutilizable.
#
# Recursos prácticos:
#
# - **[CellGuide (CZ CELLxGENE)](https://cellxgene.cziscience.com/cellguide)** — para
#   cada tipo celular: definición, marcadores canónicos, marcadores calculados sobre
#   millones de células, y la posición en la jerarquía CL.
# - **[CZ CELLxGENE Discover](https://cellxgene.cziscience.com/)** — buscar la expresión
#   de un gen a través de decenas de millones de células anotadas.
# - **[Human Cell Atlas](https://data.humancellatlas.org/)** — datos de referencia por
#   tejido.

# Las referencias de celldex traen las etiquetas de Cell Ontology. Monaco, que es
# la que usamos para el detalle fino, trae un identificador CL por etiqueta:
onto_monaco <- unique(data.frame(
  main = ref_monaco$label.main,
  fine = ref_monaco$label.fine,
  ont  = ref_monaco$label.ont
))
print(onto_monaco, row.names = FALSE)

# Fíjate en que dos etiquetas distintas pueden compartir identificador (Vd2 y
# non-Vd2 gd T cells son ambas CL:0000798, gamma-delta T cell): la ontología tiene
# su propio nivel de detalle, que no siempre coincide con el de la referencia.

# --------------------------------------------------------------------------
## 7. La anotación final -----------------------------------------------------
#
# Combinamos las tres fuentes de evidencia y —esto es lo que distingue una anotación
# seria de una automática— **declaramos el nivel de confianza**.
#
# La regla: la etiqueta la pones tú, no el software. SingleR y Azimuth son evidencia;
# los marcadores canónicos son evidencia; el criterio biológico es tuyo.

# Tabla de trabajo: qué dice cada método para cada cluster.
# Indexamos POR NOMBRE, nunca por posición: el orden de las filas de SingleR
# no tiene por qué coincidir con el de levels().
cl <- levels(pbmc$clusters)

# Para cada cluster, la etiqueta mayoritaria cuando se anota célula por célula,
# y qué % de las células del cluster la comparte (coherencia).
mayoritaria <- sapply(cl, function(k) {
  x <- pbmc$singler_celula[pbmc$clusters == k]
  names(which.max(table(x)))
})
coherencia <- sapply(cl, function(k) {
  x <- pbmc$singler_celula[pbmc$clusters == k]
  round(100 * mean(!is.na(x) & x == mayoritaria[k]), 1)
})

# Monaco da etiquetas finas; su etiqueta principal (label.main) sirve para
# compararla con Blueprint al mismo nivel de detalle.
fina_a_main <- setNames(ref_monaco$label.main, ref_monaco$label.fine)
fina_a_cl   <- setNames(ref_monaco$label.ont,  ref_monaco$label.fine)

resumen <- data.frame(
  cluster       = cl,
  n_celulas     = as.integer(table(pbmc$clusters)[cl]),
  bp_cluster    = anot_cluster[cl, "pruned.labels"],   # Blueprint, por cluster
  bp_celula     = mayoritaria,                         # Blueprint, mayoría por célula
  coherencia    = coherencia,                          # % de células con esa mayoría
  monaco_fina   = anot_fina[cl, "pruned.labels"],      # Monaco, por cluster
  stringsAsFactors = FALSE
)
resumen$monaco_main <- unname(fina_a_main[resumen$monaco_fina])

# Top 3 marcadores propios de cada cluster (mismo procedimiento que la clase 2)
top_mk <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.4,
                         logfc.threshold = 1, verbose = FALSE) |>
  group_by(cluster) |> slice_max(avg_log2FC, n = 3) |>
  summarise(marcadores = paste(gene, collapse = ", ")) |>
  mutate(cluster = as.character(cluster))

resumen <- left_join(resumen, top_mk, by = "cluster")
print(resumen[, c("cluster", "n_celulas", "bp_cluster", "bp_celula", "coherencia",
                  "monaco_fina", "marcadores")], row.names = FALSE)

# > **Antes de seguir, lee la tabla fila por fila.** Tres casos que conviene encontrar:
# >
# > - Un cluster donde **Blueprint por cluster y Blueprint por célula no coinciden**.
#   Es exactamente la situación de §3.1: el promedio del cluster dice una cosa y la
#   mayoría de sus células, otra. ¿A quién le crees? Mira los marcadores.
# > - Un cluster donde **Blueprint dice "Monocytes" y Monaco dice "dendritic cells"**.
#   Blueprint tiene una etiqueta "DC", pero con pocos perfiles; ante la duda, asigna
#   el vecino más poblado. Es el modo de falla de la referencia incompleta.
# > - Clusters chicos donde **Blueprint devuelve `NA`**. No es un error: la referencia
#   no tiene plasmablastos ni DC plasmacitoides, y SingleR lo dice en vez de inventar.

# ---------------------------------------------------------------------------
# Nivel de confianza, con una regla explícita
#
# La confianza NO la decide el software, pero sí se puede calcular con una regla
# que cualquiera pueda auditar. Usamos tres comprobaciones independientes, cada
# una vale un punto:
#
#   1. acuerdo_bp   Blueprint por cluster == Blueprint mayoritario por célula
#   2. coherente    al menos el 80 % de las células comparte esa mayoría
#   3. acuerdo_ref  Blueprint y Monaco coinciden en el tipo principal
#                   (comparamos los nombres sin mayúsculas, guiones ni espacios:
#                   "CD4+ T-cells" y "CD4+ T cells" son lo mismo)
#
#   3 puntos = alta · 2 = media · 0-1 = baja
#
# "Baja" no quiere decir "incorrecta": quiere decir que los métodos automáticos
# no convergen, y que la etiqueta descansa en tu criterio (marcadores). Esas son
# las que hay que defender —y las que el jueves vamos a poder contrastar con la
# proteína.
# ---------------------------------------------------------------------------
normalizar <- function(x) gsub("[^a-z0-9]", "", tolower(x))

resumen <- resumen |>
  mutate(
    acuerdo_bp  = !is.na(bp_cluster) & bp_cluster == bp_celula,
    coherente   = coherencia >= 80,
    acuerdo_ref = !is.na(bp_cluster) & !is.na(monaco_main) &
                  normalizar(bp_cluster) == normalizar(monaco_main),
    puntos      = acuerdo_bp + coherente + acuerdo_ref,
    confianza   = factor(c("baja", "baja", "media", "alta")[puntos + 1],
                         levels = c("alta", "media", "baja"))
  )

print(resumen[, c("cluster", "bp_cluster", "bp_celula", "monaco_main",
                  "acuerdo_bp", "coherente", "acuerdo_ref", "confianza")],
      row.names = FALSE)

# > **Cómo leer la columna de confianza.** No todas las "bajas" significan lo mismo:
# >
# > - **El cluster que Blueprint llama CD4 y sus células llaman CD8** falla las tres
#   comprobaciones. Es un conflicto real entre métodos: el promedio del cluster se
#   parece a CD4, pero el 79 % de sus células, Monaco y los marcadores (`CD8B`,
#   `LINC02446`, `NELL2`) dicen T CD8 naive. Es la etiqueta que más hay que
#   defender, y la primera que vamos a contrastar el jueves con la proteína CD8.
# > - **Plasmablastos y DC plasmacitoides** son "baja" por otro motivo: Blueprint
#   no tiene esos tipos celulares y devuelve `NA`. Monaco y los marcadores coinciden
#   sin ambigüedad. La regla penaliza el hueco de la referencia, y está bien que lo
#   haga: la etiqueta se apoya en menos fuentes independientes.
# > - **MAIT y T CD8 memoria** quedan en "media" porque Monaco los ubica en "T cells"
#   (T no convencionales: MAIT, γδ) y Blueprint en "CD8+ T-cells". Parece una
#   diferencia de vocabulario, pero esconde biología real.

# ---------------------------------------------------------------------------
# Edita ESTE VECTOR. Es el paso donde la anotación pasa a ser tuya.
# Mira la tabla de arriba, el DotPlot de marcadores canónicos y los diagnósticos,
# y escribe la etiqueta que tú defiendes para cada cluster.
#
# Los nombres del vector son los identificadores de cluster tal como los
# devolvió Leiden en la clase 2 (empiezan en 1). Si tu clustering salió
# distinto, ajusta el vector a TU tabla: el código de abajo avisa si sobra o
# falta alguno.
# ---------------------------------------------------------------------------
etiquetas <- c(
  "1"  = "Monocito CD14+",       # Blueprint + Monaco + CD14/LYZ/S100A8
  "2"  = "T CD4 memoria",        # Monaco Th1/Th17; AQP3, TNFRSF25; sin CCR7
  "3"  = "T CD4 naive",          # CCR7, LEF1; Monaco naive CD4
  "4"  = "NK",                   # GNLY, FGFBP2: NK CD56dim
  "5"  = "T CD8 memoria",        # CD8A/B, GZMK; Monaco ve γδ (TRGC2): mezcla
  "6"  = "T CD8 naive",          # Blueprint por cluster dice CD4; las células y Monaco, CD8
  "7"  = "MAIT",                 # SLC4A10, KLRB1; Monaco MAIT
  "8"  = "B naive",              # TCL1A, IGHD, FCER2
  "9"  = "B memoria",            # IGHG1, TNFRSF13B; Monaco "Exhausted B" es engañoso
  "10" = "NK CD56bright",        # KLRC1, XCL2
  "11" = "Monocito CD16+",       # CDKN1C, HES4; FCGR3A en el DotPlot
  "12" = "DC convencional",      # FCER1A, CD1C: cDC2. Blueprint dice monocito
  "13" = "Plasmablasto",         # IGHA1/2, MZB1; ciclo celular (CDC20)
  "14" = "DC plasmacitoide"      # CLEC4C, LILRA4; Blueprint no tiene pDC
)

# Control: ¿el vector cubre exactamente los clusters del objeto?
sobran    <- setdiff(names(etiquetas), levels(pbmc$clusters))
faltantes <- setdiff(levels(pbmc$clusters), names(etiquetas))
if (length(sobran)) {
  cat("Etiquetas para clusters que NO existen en el objeto:",
      paste(sobran, collapse = ", "), "— revisa el vector\n")
  etiquetas <- etiquetas[setdiff(names(etiquetas), sobran)]
}
# Cualquier cluster sin etiqueta queda como "Sin asignar" — y eso es un
# resultado legítimo, no una falla.
if (length(faltantes)) {
  etiquetas[faltantes] <- "Sin asignar"
  cat("Clusters sin etiqueta definida:", paste(faltantes, collapse = ", "), "\n")
}

# ---------------------------------------------------------------------------
# Identificador de Cell Ontology para cada etiqueta (ver §6).
# Los de las poblaciones que Monaco tiene salen de su columna label.ont; el resto
# se buscó en CellGuide / OLS. Verifica cualquiera en
# https://www.ebi.ac.uk/ols4/ontologies/cl
# ---------------------------------------------------------------------------
cl_ontologia <- c(
  "Monocito CD14+"   = "CL:0000860",   # classical monocyte
  "T CD4 memoria"    = "CL:0000897",   # CD4-positive, alpha-beta memory T cell
  "T CD4 naive"      = "CL:0000895",   # naive thymus-derived CD4-positive, alpha-beta T cell
  "NK"               = "CL:0000939",   # CD16-positive, CD56-dim natural killer cell, human
  "T CD8 memoria"    = "CL:0000913",   # effector memory CD8-positive, alpha-beta T cell
  "T CD8 naive"      = "CL:0000900",   # naive thymus-derived CD8-positive, alpha-beta T cell
  "MAIT"             = "CL:0000940",   # mucosal-associated invariant T cell
  "B naive"          = "CL:0000788",   # naive B cell
  "B memoria"        = "CL:0000787",   # memory B cell
  "NK CD56bright"    = "CL:0000938",   # CD16-negative, CD56-bright natural killer cell, human
  "Monocito CD16+"   = "CL:0000875",   # non-classical monocyte
  "DC convencional"  = "CL:0002399",   # CD1c-positive myeloid dendritic cell
  "Plasmablasto"     = "CL:0000980",   # plasmablast
  "DC plasmacitoide" = "CL:0000784",   # plasmacytoid dendritic cell
  "Sin asignar"      = NA
)

# La tabla final: una fila por cluster, con la etiqueta, su CL, la confianza y
# la evidencia. Es lo que va al material suplementario de un paper.
anotacion <- resumen |>
  transmute(cluster, n_celulas,
            etiqueta     = unname(etiquetas[cluster]),
            cl_id        = unname(cl_ontologia[etiqueta]),
            confianza,
            cl_monaco    = unname(fina_a_cl[monaco_fina]),
            bp_cluster, bp_celula, coherencia, monaco_fina, marcadores)
print(anotacion[, c("cluster", "etiqueta", "cl_id", "confianza", "cl_monaco")],
      row.names = FALSE)

write.csv(anotacion, file.path(dir_res, "clase3_anotacion_por_cluster.csv"),
          row.names = FALSE)

# Mismo cuidado que en la sección 3.1: `factor()` NO borra los nombres, así que
# hay que renombrar por célula antes de asignar a la metadata.
id_cl <- as.character(pbmc$clusters)
pbmc$tipo_celular <- setNames(factor(etiquetas[id_cl], levels = unique(etiquetas)),
                              colnames(pbmc))
pbmc$cl_ontologia <- setNames(unname(cl_ontologia[etiquetas[id_cl]]), colnames(pbmc))
pbmc$confianza_anotacion <- setNames(
  resumen$confianza[match(id_cl, resumen$cluster)], colnames(pbmc))
table(pbmc$tipo_celular)

# Los mismos colores que en las diapositivas de la clase
COL_CONF <- c(alta = "#0F7173", media = "#E0A526", baja = "#F05D5E")

(DimPlot(pbmc, group.by = "tipo_celular", label = TRUE, repel = TRUE,
         label.size = 3.5) + ggtitle("Anotación final") + NoLegend()) |
(DimPlot(pbmc, group.by = "confianza_anotacion", cols = COL_CONF) +
   ggtitle("Confianza declarada"))

# Verificación: los marcadores canónicos deben ser coherentes con las etiquetas
DotPlot(pbmc, features = marcadores_canonicos, group.by = "tipo_celular") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8), strip.text = element_text(size = 7, angle = 60)) +
  labs(title = "Control final: ¿cada etiqueta enciende los marcadores que le corresponden?")

saveRDS(pbmc, ck3)
tam_archivo(ck3)

# --------------------------------------------------------------------------
## 8. Síntesis ---------------------------------------------------------------
#
# 1. **La anotación es una hipótesis**, no un resultado. Se sostiene con evidencia
#    convergente.
# 2. **Los marcadores canónicos vienen de la proteína**; el ARNm es un proxy imperfecto
#    (dropout, ARNm ≠ proteína, expresión compartida).
# 3. **`celldex` cambió de API**: `fetchReference(nombre, versión)`. Las referencias
#    ahora están versionadas y mapeadas a Cell Ontology.
# 4. **Nunca mires solo la etiqueta**: score, delta y `pruned.labels` son parte del
#    resultado. Un `NA` es información.
# 5. **Anota por célula y por cluster.** Donde no coinciden, hay algo que entender.
# 6. **El RNA ambiental produce anotaciones falsas.** Diagnostica siempre; corrige con
#    cuidado y compara.
# 7. **Declara la confianza con una regla explícita.** "Baja" no es "incorrecta":
#    es donde los métodos no convergen y la etiqueta depende de tu criterio.
# 8. **Usa vocabulario controlado** (Cell Ontology) para que tu anotación sirva más allá
#    de tu propio análisis.
#
# El jueves agregamos una modalidad que ataca directamente el punto 2: **medir la
# proteína de superficie en las mismas células**.
