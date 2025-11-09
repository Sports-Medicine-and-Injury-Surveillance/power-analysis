## -- Load packages -------------------------------------------------------
library(shiny)
library(ggplot2)

## -- Source power analysis function -------------------------------------------------------
source("function.R")

## -- Shiny app -------------------------------------------------------
ui <- navbarPage(
  "Power Analysis Tool",

  # ---- TAB 1: Power Calculator ----
  tabPanel(
    "Power Calculator",
    sidebarLayout(
      sidebarPanel(
        # ---- Inputs for base calculation ----
        sliderInput(
          "sd",
          "Standard deviation (sd):",
          min = 1,
          max = 10,
          value = 3.8,
          step = 0.1
        ),
        sliderInput(
          "rsq",
          "R-squared (rsq):",
          min = 0,
          max = 0.9,
          value = 0.5,
          step = 0.05
        ),
        sliderInput(
          "alpha",
          "Significance level (alpha):",
          min = 0.01,
          max = 0.1,
          value = 0.05,
          step = 0.005
        ),
        sliderInput(
          "n_groups",
          "Number of groups:",
          min = 2,
          max = 6,
          value = 4,
          step = 1
        ),
        sliderInput(
          "dropout_rate",
          "Dropout rate:",
          min = 0,
          max = 0.5,
          value = 0.18,
          step = 0.01
        ),

        # ---- Power analysis-specific inputs ----
        numericInput(
          "effect_points_calc",
          "Effect size (points):",
          value = 2.0,
          min = 0.1,
          max = 10,
          step = 0.1
        ),
        sliderInput(
          "power_target",
          "Target power:",
          min = 0.5,
          max = 0.99,
          value = 0.80,
          step = 0.01
        ),

        # ---- Method selection ----
        selectInput(
          "method_calc",
          "Method:",
          choices = c(T2 = "t2", F2 = "f2"),
          selected = "t2"
        )
      ),
      mainPanel(
        h4("Power Analysis Summary"),
        verbatimTextOutput("power_summary")
      )
    )
  ),

  # ---- TAB 2: Power Plot ----
  tabPanel(
    "Power Curve Plot",
    sidebarLayout(
      sidebarPanel(
        # ---- Show inherited parameters from Tab 1 ----
        h5(strong("Parameters from Power Calculator")),
        verbatimTextOutput("param_summary"),
        tags$hr(),

        # ---- Plot configuration ----
        checkboxGroupInput(
          "method_plot",
          "Methods to display:",
          choices = c(T2 = "t2", F2 = "f2"),
          selected = c("t2", "f2")
        ),
        sliderInput(
          "effect_points_range",
          "Effect size range (points):",
          min = 0.5,
          max = 5,
          value = c(1.5, 3),
          step = 0.1
        ),
        numericInput(
          "effect_step",
          "Effect size step:",
          value = 0.5,
          min = 0.1,
          max = 1,
          step = 0.1
        ),
        sliderInput(
          "sample_range",
          "Sample size range:",
          min = 20,
          max = 300,
          value = c(40, 120),
          step = 5
        ),

        helpText(
          "Plot uses parameters (sd, rsq, alpha, n_groups, dropout_rate) from the Power Calculator tab."
        )
      ),
      mainPanel(
        h4("Power Curve Plot"),
        plotOutput("power_plot", height = "750px")
      )
    )
  )
)

server <- function(input, output, session) {
  # ---- Reactive PowerAnalysis object ----
  pa_obj <- reactive({
    PowerAnalysis$new(
      method = input$method_calc,
      sd = input$sd,
      rsq = input$rsq,
      alpha = input$alpha,
      n_groups = input$n_groups,
      dropout_rate = input$dropout_rate,
      effect_points = input$effect_points_calc,
      power_target = input$power_target
    )
  })

  # ---- Summary output ----
  output$power_summary <- renderPrint({
    pa <- pa_obj()
    pa$compute()
    pa$print()
  })

  # ---- Parameter summary for second tab ----
  output$param_summary <- renderPrint({
    cat("SD:", input$sd, "\n")
    cat("R-squared:", input$rsq, "\n")
    cat("Alpha:", input$alpha, "\n")
    cat("Groups:", input$n_groups, "\n")
    cat("Dropout rate:", input$dropout_rate, "\n")
  })

  # ---- Plot output ----
  output$power_plot <- renderPlot({
    pa <- pa_obj()

    sample_seq <- seq(input$sample_range[1], input$sample_range[2], by = 5)
    effect_seq <- seq(
      input$effect_points_range[1],
      input$effect_points_range[2],
      by = input$effect_step
    )

    pa$plot(
      sample_seq = sample_seq,
      effect_seq = effect_seq,
      methods = input$method_plot,
      sd = input$sd,
      rsq = input$rsq,
      alpha = input$alpha,
      n_groups = input$n_groups,
      dropout_rate = input$dropout_rate
    )
  })
}

shinyApp(ui, server)
