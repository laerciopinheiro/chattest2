# ====== PACOTES ======
# Instale uma vez se necessário:
# install.packages(c(
#   "readxl", "dplyr", "tidyr", "stringr", "ggplot2", "purrr", "boot"
# ))

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(boot)

# ====== FUNÇÕES AUXILIARES ======
limpa_colnames <- function(df) {
  names(df) <- names(df) %>%
    str_trim() %>%
    str_replace("\\.0$", "") %>%
    str_replace("\\.\\.\\.[0-9]+$", "")
  df
}

limpa_tabela_prenatal <- function(path) {
  stopifnot(file.exists(path))
  anos_chr <- as.character(2014:2023)

  read_excel(path) %>%
    limpa_colnames() %>%
    filter(!str_detect(`Região/Unidade da Federação`, "Região")) %>%
    rename(UF = `Região/Unidade da Federação`) %>%
    select(any_of(c("UF", anos_chr, "Total"))) %>%
    mutate(UF = str_trim(UF))
}

limpa_nv_total <- function(path) {
  stopifnot(file.exists(path))
  anos_chr <- as.character(2014:2023)

  read_excel(path) %>%
    limpa_colnames() %>%
    filter(!str_detect(`Região/Unidade da Federação`, "Região")) %>%
    rename(UF = `Região/Unidade da Federação`) %>%
    select(any_of(c("UF", anos_chr, "Total"))) %>%
    mutate(UF = str_trim(UF))
}

padroniza_uf <- function(x) {
  x %>%
    str_replace("^\\.\\.", "") %>%
    str_trim()
}

calc_spearman_boot <- function(df, xvar, yvar, label, n_boot = 5000, seed = 123) {
  df2 <- df %>% drop_na(.data[[xvar]], .data[[yvar]])
  if (nrow(df2) < 3) {
    return(tibble(
      label = label,
      xvar = xvar,
      rho = NA_real_,
      p_value = NA_real_,
      n = nrow(df2),
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }

  spearman <- suppressWarnings(cor.test(df2[[xvar]], df2[[yvar]], method = "spearman"))

  boot_fun <- function(data, indices) {
    d <- data[indices, , drop = FALSE]
    suppressWarnings(cor(d[[xvar]], d[[yvar]], method = "spearman"))
  }

  set.seed(seed)
  boot_res <- boot::boot(df2, statistic = boot_fun, R = n_boot)
  ci <- quantile(boot_res$t, c(0.025, 0.975), na.rm = TRUE)

  tibble(
    label = label,
    xvar = xvar,
    rho = unname(spearman$estimate),
    p_value = spearman$p.value,
    n = nrow(df2),
    ci_low = ci[1],
    ci_high = ci[2]
  )
}

plot_spearman_scatter <- function(df, uf, xvar, xlabel, spearman_row) {
  annotation <- NULL
  if (!is.null(spearman_row) && !is.na(spearman_row$rho) && !is.na(spearman_row$p_value)) {
    rho_fmt <- sprintf("%.3f", spearman_row$rho)
    p_fmt <- format.pval(spearman_row$p_value, digits = 3, eps = .001)
    annotation <- sprintf("ρ = %s\np = %s", rho_fmt, p_fmt)
  }

  df %>%
    filter(UF == uf) %>%
    ggplot(aes(x = .data[[xvar]], y = Taxa_mortalidade_0_27_por1000)) +
    geom_point(size = 2.5, colour = "#1f77b4") +
    geom_smooth(method = "lm", se = TRUE, colour = "#ff7f0e", fill = scales::alpha("#ff7f0e", 0.2)) +
    geom_text(aes(label = Ano), nudge_y = 0.1, size = 3) +
    {
      if (!is.null(annotation)) {
        annotate(
          "label",
          x = -Inf,
          y = Inf,
          label = annotation,
          hjust = -0.05,
          vjust = 1.1,
          size = 3,
          label.size = 0,
          fill = scales::alpha("white", 0.6)
        )
      } else {
        NULL
      }
    } +
    labs(
      x = xlabel,
      y = "Mortalidade 0–27 dias (por 1.000 NV)",
      title = sprintf("%s: Mortalidade vs %s", uf, xlabel)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

plot_dot_ci <- function(spearman_df) {
  if (nrow(spearman_df) == 0) {
    return(NULL)
  }

  spearman_df %>%
    mutate(grupo = factor(paste(UF, label, sep = " · "), levels = rev(paste(UF, label, sep = " · ")))) %>%
    ggplot(aes(x = rho, y = grupo)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, colour = "#636363") +
    geom_point(size = 3, colour = "#2ca02c") +
    labs(
      x = "ρ de Spearman (com IC95% bootstrap)",
      y = NULL,
      title = "Correlação de Spearman por indicador"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# ====== FUNÇÃO PRINCIPAL ======
run_analysis <- function(paths, out_dir = "output", n_boot = 5000, seed = 123) {
  stopifnot(is.list(paths))
  required <- c("taxas", "sem", "pn_1a3", "pn_4a6", "pn_7mais", "nv_total")
  missing_paths <- setdiff(required, names(paths))
  if (length(missing_paths) > 0) {
    stop("Caminhos ausentes: ", paste(missing_paths, collapse = ", "))
  }

  missing_files <- names(paths)[!file.exists(paths)]
  if (length(missing_files) > 0) {
    formatted <- paste0("  ", missing_files, ": ", unlist(paths[missing_files]))
    stop(
      "Os arquivos abaixo não foram encontrados:\n",
      paste(formatted, collapse = "\n"),
      call. = FALSE
    )
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  sem_pn <- limpa_tabela_prenatal(paths$sem)
  pn_1a3 <- limpa_tabela_prenatal(paths$pn_1a3)
  pn_4a6 <- limpa_tabela_prenatal(paths$pn_4a6)
  pn_7mais <- limpa_tabela_prenatal(paths$pn_7mais)
  nv_total <- limpa_nv_total(paths$nv_total)

  sem_pn$UF <- padroniza_uf(sem_pn$UF)
  pn_1a3$UF <- padroniza_uf(pn_1a3$UF)
  pn_4a6$UF <- padroniza_uf(pn_4a6$UF)
  pn_7mais$UF <- padroniza_uf(pn_7mais$UF)
  nv_total$UF <- padroniza_uf(nv_total$UF)

  taxas_raw <- read_excel(paths$taxas)
  taxas_br_pi <- taxas_raw %>%
    transmute(
      Ano = as.integer(Ano),
      Brasil = `Taxa Brasil (por 1000)`,
      `Piauí` = `Taxa Piauí (por 1000)`
    ) %>%
    pivot_longer(-Ano, names_to = "UF", values_to = "Taxa_mortalidade_0_27_por1000")

  anos_chr <- as.character(2014:2023)

  nv_long <- nv_total %>%
    pivot_longer(all_of(anos_chr), names_to = "Ano", values_to = "NV_total") %>%
    mutate(Ano = as.integer(Ano))

  sem_long <- sem_pn %>%
    pivot_longer(all_of(anos_chr), names_to = "Ano", values_to = "n_sem") %>%
    mutate(Ano = as.integer(Ano))

  mais7_long <- pn_7mais %>%
    pivot_longer(all_of(anos_chr), names_to = "Ano", values_to = "n_7mais") %>%
    mutate(Ano = as.integer(Ano))

  prenatal_prop <- nv_long %>%
    left_join(select(sem_long, UF, Ano, n_sem), by = c("UF", "Ano")) %>%
    left_join(select(mais7_long, UF, Ano, n_7mais), by = c("UF", "Ano")) %>%
    mutate(
      pct_sem = 100 * n_sem / NV_total,
      pct_7mais = 100 * n_7mais / NV_total
    )

  dados_br_pi <- prenatal_prop %>%
    filter(UF %in% c("TOTAL", "Piauí")) %>%
    mutate(UF = if_else(UF == "TOTAL", "Brasil", UF)) %>%
    left_join(taxas_br_pi, by = c("UF", "Ano"))

  spearman_targets <- tribble(
    ~xvar, ~label,
    "pct_7mais", "Pré-natal ≥7 consultas (%)",
    "pct_sem", "Sem pré-natal (%)"
  )

  spearman_results <- dados_br_pi %>%
    group_by(UF) %>%
    group_map(~{
      map_dfr(1:nrow(spearman_targets), function(i) {
        target <- spearman_targets[i, ]
        calc_spearman_boot(
          df = .x,
          xvar = target$xvar,
          yvar = "Taxa_mortalidade_0_27_por1000",
          label = target$label,
          n_boot = n_boot,
          seed = seed
        )
      }) %>%
        mutate(UF = unique(.x$UF))
    }) %>%
    list_rbind()

  valid_spearman <- spearman_results %>% filter(!is.na(rho))

  dot_plot <- plot_dot_ci(valid_spearman)
  if (!is.null(dot_plot)) {
    ggsave(
      filename = file.path(out_dir, "spearman_dotplot.png"),
      plot = dot_plot,
      width = 7,
      height = 5,
      dpi = 300
    )
  }

  scatter_plots <- list()
  for (uf in unique(dados_br_pi$UF)) {
    df_uf <- filter(dados_br_pi, UF == uf)
    for (i in seq_len(nrow(spearman_targets))) {
      target <- spearman_targets[i, ]
      spearman_row <- valid_spearman %>%
        filter(UF == uf, xvar == target$xvar) %>%
        slice_head(n = 1)
      if (nrow(spearman_row) == 0) {
        next
      }
      p <- plot_spearman_scatter(df_uf, uf, target$xvar, target$label, spearman_row)
      fname <- sprintf(
        "scatter_%s_%s.png",
        stringr::str_replace_all(tolower(uf), "[^a-z0-9]+", "_"),
        ifelse(target$xvar == "pct_7mais", "pn7mais", "semprenatal")
      )
      ggsave(
        filename = file.path(out_dir, fname),
        plot = p,
        width = 7,
        height = 5,
        dpi = 300
      )
      scatter_plots[[paste(uf, target$xvar, sep = "_")]] <- p
    }
  }

  list(
    dados = dados_br_pi,
    spearman = spearman_results,
    dot_plot = dot_plot,
    scatter_plots = scatter_plots
  )
}

# ====== EXEMPLO DE USO ======
# paths <- list(
#   taxas = "C:/caminho/para/Taxas_de_Mortalidade_Neonatal__0_27_dias__Brasil_x_Piau__2014_2023 (1).xlsx",
#   sem = "C:/caminho/para/28._Nascim ... (1).xlsx",
#   pn_1a3 = "C:/caminho/para/27. _Nascim ... (1).xlsx",
#   pn_4a6 = "C:/caminho/para/26. _Nascim ... (1).xlsx",
#   pn_7mais = "C:/caminho/para/25. _Nascim ... (1).xlsx",
#   nv_total = "C:/caminho/para/15. Nascim ... (2).xlsx"
# )
# resultado <- run_analysis(paths, out_dir = "C:/Users/laerc/Desktop/Trabalho da Tolstenko", n_boot = 5000)
# View(resultado$spearman)
