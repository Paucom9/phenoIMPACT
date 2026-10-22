# Regenerate Figure 1 panels from LOCAL data, and update the supplied PPTX.
# No environmental/demographic model is refitted. Panel b refits only example GAMs.
# source('rebuild_figure1_nature.R')
# fig1 <- rebuild_figure1_nature(pptx_file = file.choose())
rebuild_figure1_nature <- function(
  project_root = "E:/phenoIMPACT project/code/phenoIMPACT",
  pptx_file = file.choose()) {
  pkgs <- c("ggplot2", "xml2", "zip", "dplyr", "tidyr", "here", "sf",
            "rnaturalearth", "data.table", "mgcv", "lubridate", "purrr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
  stopifnot(file.exists(pptx_file), capabilities("cairo"))
  pptx_file <- normalizePath(pptx_file, winslash = "/")
  previous <- getwd(); on.exit(setwd(previous), add = TRUE); setwd(project_root)
  out <- file.path(project_root, "output", "figures", paste0("figure1_nature_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(out, recursive = TRUE)
  # Evaluate source code on its own, removing session clearing and obsolete exports.
  # Original scripts, figures and data remain unchanged.
  run_panel <- function(script, object) {
    txt <- readLines(file.path(project_root, script), warn = FALSE)
    cut <- grep("^#### Save figure", txt)
    if (length(cut)) txt <- txt[seq_len(cut[length(cut)]-1L)]
    txt <- txt[!grepl("^rm\\(list = ls\\(\\)\\)", txt) & !grepl("^library\\(extrafont\\)", txt)]
    txt <- gsub('"Garamond"', '"Arial"', txt, fixed = TRUE)
    # The map's intermediate summary tables also go to the new export folder.
    txt <- gsub('here::here("output", "figures", "fig1_map_tables")',
                'file.path(export_directory, "map_tables")', txt, fixed = TRUE)
    env <- new.env(parent = globalenv()); env$export_directory <- out
    eval(parse(text = txt), envir = env)
    get(object, envir = env)
  }
  message("Regenerating map from local site data")
  a <- run_panel("map_study_sites.R", "fig1a_map_clean")
  message("Regenerating annual curves for the two existing example populations")
  b <- run_panel("pheno_curve_examples.R", "fig1b_pheno_anomaly_curves") +
    ggplot2::labs(colour = "Temperature anomaly (SD)")
  message("Regenerating conceptual hypotheses")
  c <- run_panel("plot_conceptual_plasticity.R", "fig1c_conceptual")
  c$data <- c$data # Layer data carry the facet labels in this source.
  rename_labels <- function(d) {
    if (is.data.frame(d) && "moderator" %in% names(d)) {
      x <- as.character(d$moderator)
      x[x == "Daylength"] <- "Photoperiod"
      x[x == "Mean temperature"] <- "Background temperature"
      x[x == "Temperature trend"] <- "Warming trend"
      d$moderator <- factor(x, levels = c("Photoperiod", "Background temperature",
        "Temperature predictability", "Warming trend"))
    }
    d
  }
  for (j in seq_along(c$layers)) c$layers[[j]]$data <- rename_labels(c$layers[[j]]$data)
  c <- c + ggplot2::labs(subtitle = "Hypothesized responses") +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(family = "Arial", size = 11))
  plots <- list(a = a, b = b, c = c)
  dimensions <- list(a = c(5.8, 6.2), b = c(5.5, 7.2), c = c(9.5, 5.2))
  for (name in names(plots)) {
    size <- dimensions[[name]]
    ggplot2::ggsave(file.path(out, paste0("fig1", name, ".png")), plots[[name]],
      width = size[1], height = size[2], dpi = 600, bg = "white")
    ggplot2::ggsave(file.path(out, paste0("fig1", name, ".pdf")), plots[[name]],
      width = size[1], height = size[2], device = grDevices::cairo_pdf, bg = "white")
  }
  # Reuse the PPTX composition and butterfly photographs unchanged.
  temp <- tempfile("fig1_pptx_"); dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE), add = TRUE)
  utils::unzip(pptx_file, exdir = temp)
  slide_path <- file.path(temp, "ppt", "slides", "slide1.xml")
  rel_path <- file.path(temp, "ppt", "slides", "_rels", "slide1.xml.rels")
  slide <- xml2::read_xml(slide_path); rel <- xml2::read_xml(rel_path)
  ns <- c(p = "http://schemas.openxmlformats.org/presentationml/2006/main",
    a = "http://schemas.openxmlformats.org/drawingml/2006/main",
    r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships")
  # Original supplied figure has duplicate conceptual/curve images; replace all.
  replacements <- c("image1.png" = "c", "image4.png" = "c", "image2.png" = "b",
                     "image7.png" = "b", "image3.png" = "a")
  relation_nodes <- xml2::xml_find_all(rel, "//*[local-name()='Relationship']")
  for (node in relation_nodes) {
    target <- xml2::xml_attr(node, "Target"); base <- basename(target)
    if (base %in% names(replacements)) {
      destination <- file.path(temp, "ppt", "media", base)
      file.copy(file.path(out, paste0("fig1", replacements[[base]], ".png")), destination, overwrite = TRUE)
    }
  }
  for (node in xml2::xml_find_all(slide, ".//a:latin | .//a:ea | .//a:cs", ns))
    xml2::xml_set_attr(node, "typeface", "Arial")
  # Explicit font on every editable run, including panel letters and labels.
  for (node in xml2::xml_find_all(slide, ".//a:rPr | .//a:defRPr | .//a:endParaRPr", ns)) {
    latin <- xml2::xml_find_first(node, "./a:latin", ns)
    if (inherits(latin, "xml_missing"))
      xml2::xml_add_child(node, '<a:latin xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" typeface="Arial"/>')
    else xml2::xml_set_attr(latin, "typeface", "Arial")
  }
  for (node in xml2::xml_find_all(slide, ".//a:t", ns)) {
    t <- xml2::xml_text(node)
    t <- gsub("Unvoltine", "Univoltine", t, fixed = TRUE)
    t <- gsub("Multivolotine", "Multivoltine", t, fixed = TRUE)
    xml2::xml_set_text(node, t)
  }
  # Thin the large outer panel boxes to 1 pt; preserve plot borders and artwork.
  for (shape in xml2::xml_find_all(slide, ".//p:sp", ns)) {
    ext <- xml2::xml_find_first(shape, "./p:spPr/a:xfrm/a:ext", ns)
    line <- xml2::xml_find_first(shape, "./p:spPr/a:ln", ns)
    if (!inherits(ext,"xml_missing") && !inherits(line,"xml_missing")) {
      if (as.numeric(xml2::xml_attr(ext,"cx")) > 3000000 &&
          as.numeric(xml2::xml_attr(ext,"cy")) > 2000000)
        xml2::xml_set_attr(line, "w", "12700")
    }
  }
  xml2::write_xml(slide, slide_path)
  output_pptx <- file.path(out, "Figura1_phenoIMPACT_Arial.pptx")
  zip::zipr(output_pptx, list.files(temp, all.files = TRUE, no.. = TRUE), root = temp)
  writeLines(c("Arial typography throughout regenerated chart panels and editable labels.",
    "Panel b anomaly is standardized: units corrected to SD.",
    "Panel c retains original predictions explicitly labelled as hypotheses, not final results.",
    "Conceptual line type identifies voltinism; explain this in the Figure 1 caption.",
    "Source PPTX composition and butterfly photographs retained. Original files unchanged.",
    "PowerPoint uses 600-dpi raster panels; separate PDFs provide vector panel exports.",
    "Inspect the new PowerPoint for label overlaps and image cropping before submission."),
    file.path(out, "README.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  zip::zipr(file.path(out, "figure1_nature_review.zip"),
    list.files(out, pattern = "\\.(pptx|pdf|txt)$", full.names = TRUE), root = out)
  message("Updated PowerPoint and panels: ", normalizePath(out, winslash = "/"))
  invisible(list(pptx = output_pptx, plots = plots, output_dir = out))
}
