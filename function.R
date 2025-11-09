library(R6)
library(pwr)
library(ggplot2)
library(dplyr)
library(shiny)

PowerAnalysis <- R6Class(
  "PowerAnalysis",
  public = list(
    # --- Public fields ---
    method = NULL, # "t2" or "f2"
    n_groups = 4,
    sd = 3.8,
    effect_points = 2.0,
    rsq = 0.5,
    alpha = 0.05,
    n_comparisons = 3,
    dropout_rate = 0.18,
    power_target = 0.80,
    n_covariates = 3,
    total_sample = NULL,
    results = list(),

    # --- Constructor ---
    initialize = function(
      method = "t2",
      n_groups = 4,
      sd = 3.8,
      effect_points = 2.0,
      rsq = 0.5,
      alpha = 0.05,
      n_comparisons = 3,
      dropout_rate = 0.18,
      power_target = 0.8,
      n_covariates = 3,
      total_sample = NULL
    ) {
      self$method <- match.arg(method, c("t2", "f2"))
      self$n_groups <- n_groups
      self$sd <- sd
      self$effect_points <- effect_points
      self$rsq <- rsq
      self$alpha <- alpha
      self$n_comparisons <- n_comparisons
      self$dropout_rate <- dropout_rate
      self$power_target <- power_target
      self$n_covariates <- n_covariates
      self$total_sample <- total_sample
      self$compute()
    },

    # --- Update method ---
    update = function(...) {
      args <- list(...)
      for (nm in names(args)) {
        if (!is.null(args[[nm]])) self[[nm]] <- args[[nm]]
      }
      # Validate numeric inputs
      stopifnot(
        self$n_groups > 1,
        self$sd > 0,
        self$effect_points > 0,
        self$alpha > 0 & self$alpha < 1,
        self$dropout_rate >= 0 & self$dropout_rate < 1
      )
      self$method <- match.arg(self$method, c("t2", "f2"))
      self$compute()
      invisible(self)
    },

    # --- Print method ---
    print = function() {
      cat("\n--- Power Analysis Summary ---\n")
      cat("Method:               ", self$method, "\n")
      cat("Groups:               ", self$n_groups, "\n")
      cat("SD:                   ", self$sd, "\n")
      cat("Effect (points):      ", self$effect_points, "\n")
      cat("Cohen's d:            ", round(self$results$cohen_d, 3), "\n")
      cat("Adjusted d:           ", round(self$results$adj_d, 3), "\n")
      cat("R-squared:            ", self$rsq, "\n")
      if (self$method == "t2") {
        cat(
          "Alpha (Bonferroni):   ",
          signif(self$results$alpha_adj, 3),
          " (",
          self$n_comparisons,
          " comparisons)\n"
        )
      } else {
        cat("Alpha:                ", self$alpha, "\n")
      }
      cat("Target power:         ", self$power_target, "\n")
      cat("Dropout rate:         ", self$dropout_rate * 100, "%\n\n")
      cat(
        "Required sample per group (completers): ",
        self$results$n_per_group,
        "\n"
      )
      cat(
        "Total sample (completers):              ",
        self$results$total_n,
        "\n"
      )
      cat(
        "Required per group (recruited):         ",
        self$results$n_recruited,
        "\n"
      )
      cat(
        "Total sample (recruited):               ",
        self$results$total_recruited,
        "\n"
      )
      if (self$method == "f2") {
        cat(
          "Effect size (f²):                       ",
          round(self$results$f2, 4),
          "\n"
        )
      }
      cat("-----------------------------------------\n")
      invisible(self)
    },

    # --- Plot method ---
    plot = function(
      sample_seq = seq(20, 80, 5),
      effect_seq = seq(1, 3, 0.5),
      methods = c("t2", "f2"),
      sd = self$sd,
      rsq = self$rsq,
      alpha = self$alpha,
      n_groups = self$n_groups,
      dropout_rate = self$dropout_rate
    ) {
      df <- private$.prepare_plot_data(
        sample_seq,
        effect_seq,
        methods,
        sd,
        rsq,
        alpha,
        n_groups,
        dropout_rate
      )
      ggplot(
        df,
        aes(
          x = total_sample,
          y = power,
          color = factor(effect_points),
          linetype = method
        )
      ) +
        geom_line(linewidth = 1) +
        geom_hline(
          yintercept = c(0.8, 0.9),
          linetype = "dashed",
          color = "gray"
        ) +
        scale_color_brewer(palette = "Set1") +
        scale_y_continuous(
          breaks = union(seq(0, 1, 0.2), c(0.8, 0.9)),
          labels = scales::label_percent(1, suffix = "")
        ) +
        scale_linetype_manual(
          labels = c('T2', 'F2'),
          breaks = c('t2', 'f2'),
          values = c(1, 2)
        ) +
        labs(
          x = "Total sample size",
          y = "Power (%)",
          color = "Effect size",
          linetype = "Method"
        ) +
        theme_grey(base_size = 16) +
        theme(
          legend.position = "bottom",
          legend.justification = "left",
          legend.key.width = unit(1, "cm")
        ) +
        expand_limits(y = c(0, 1))
    },

    # --- Core compute ---
    compute = function() {
      if (self$method == "t2") {
        private$compute_t2()
      } else {
        private$compute_f2()
      }
    }
  ),

  private = list(
    # --- Shared effect calculation ---
    .common_effect = function() {
      cohen_d <- self$effect_points / self$sd
      adj_d <- cohen_d / sqrt(1 - self$rsq)
      alpha_adj <- self$alpha / self$n_comparisons
      f2 <- (adj_d^2) / self$n_groups
      list(cohen_d = cohen_d, adj_d = adj_d, alpha_adj = alpha_adj, f2 = f2)
    },

    # --- Compute t2 ---
    compute_t2 = function() {
      vals <- private$.common_effect()
      if (!is.null(self$total_sample)) {
        n_per_group <- self$total_sample / self$n_groups
        res <- pwr::pwr.t.test(
          d = vals$adj_d,
          sig.level = vals$alpha_adj,
          n = n_per_group,
          type = "two.sample",
          alternative = "two.sided",
          power = NULL
        )
        total_n <- self$total_sample
      } else {
        res <- pwr::pwr.t.test(
          d = vals$adj_d,
          sig.level = vals$alpha_adj,
          power = self$power_target,
          type = "two.sample",
          alternative = "two.sided",
          n = NULL
        )
        n_per_group <- ceiling(res$n)
        total_n <- n_per_group * self$n_groups
      }
      n_recruited <- ceiling(n_per_group / (1 - self$dropout_rate))
      total_recruited <- n_recruited * self$n_groups

      self$total_sample <- total_n
      self$results <- c(
        vals,
        list(
          power = res$power,
          n_per_group = n_per_group,
          total_n = total_n,
          n_recruited = n_recruited,
          total_recruited = total_recruited
        )
      )
    },

    # --- Compute f2 ---
    compute_f2 = function() {
      vals <- private$.common_effect()
      if (!is.null(self$total_sample)) {
        v <- self$total_sample - self$n_groups - self$n_covariates
        res <- pwr::pwr.f2.test(
          u = 1,
          f2 = vals$f2,
          sig.level = self$alpha,
          v = v,
          power = NULL
        )
        total_n <- self$total_sample
      } else {
        res <- pwr::pwr.f2.test(
          u = 1,
          f2 = vals$f2,
          sig.level = self$alpha,
          power = self$power_target,
          v = NULL
        )
        total_n <- ceiling(res$v + self$n_groups + self$n_covariates)
      }
      n_per_group <- ceiling(total_n / self$n_groups)
      n_recruited <- ceiling(n_per_group / (1 - self$dropout_rate))
      total_recruited <- n_recruited * self$n_groups

      self$total_sample <- total_n
      self$results <- c(
        vals,
        list(
          power = res$power,
          f2 = vals$f2,
          n_per_group = n_per_group,
          total_n = total_n,
          n_recruited = n_recruited,
          total_recruited = total_recruited
        )
      )
    },

    # --- Prepare plot data ---
    .prepare_plot_data = function(
      sample_seq,
      effect_seq,
      methods,
      sd,
      rsq,
      alpha,
      n_groups,
      dropout_rate
    ) {
      grid <- expand.grid(
        total_sample = sample_seq,
        effect_points = effect_seq,
        method = methods,
        stringsAsFactors = FALSE
      )
      df <- lapply(seq_len(nrow(grid)), function(i) {
        row <- grid[i, ]
        obj <- self$clone(deep = TRUE)
        obj$update(
          total_sample = row$total_sample,
          effect_points = row$effect_points,
          method = row$method,
          sd = sd,
          rsq = rsq,
          alpha = alpha,
          n_groups = n_groups,
          dropout_rate = dropout_rate
        )
        data.frame(
          total_sample = row$total_sample,
          effect_points = row$effect_points,
          method = row$method,
          power = obj$results$power,
          stringsAsFactors = FALSE
        )
      })
      df <- do.call(rbind, df)
      df$method <- factor(df$method, levels = c("t2", "f2"))
      df
    }
  )
)

# pa <- PowerAnalysis$new(method = "t2", effect_points = 2)
# pa$print()

# # Update method and recompute
# pa$update(method = "f2", rsq = 0.4, effect_points = 1.5)
# pa$print()

# # Plot power curves
# pa$plot(sample_seq = seq(20, 100, 5), effect_seq = seq(1, 3, 0.5))
