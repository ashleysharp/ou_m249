# app.R — Interactive Distributions Explorer (Shiny)
# -------------------------------------------------
# This Shiny app turns the original R Markdown into an interactive playground
# for discrete & continuous distributions, sampling, CIs, hypothesis tests,
# and simple likelihood visualizations.

# Packages ---------------------------------------------------------------
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# Helpers ---------------------------------------------------------------
num_fmt <- function(x, d=3) format(round(x, d), nsmall=d, trim=TRUE)

plot_pmf_pdf <- function(df, x_col="x", y_col="d", discrete=TRUE, title="PMF/PDF"){
  p <- ggplot(df, aes(x=.data[[x_col]], y=.data[[y_col]])) +
    labs(x="x", y=ifelse(discrete, "p(x)", "f(x)"), title=title)
  if (discrete) p + geom_segment(aes(xend=.data[[x_col]], yend=0)) + geom_point()
  else p + geom_line()
}

plot_cdf <- function(df, x_col="x", y_col="p", discrete=TRUE, title="CDF"){
  p <- ggplot(df, aes(x=.data[[x_col]], y=.data[[y_col]])) +
    labs(x="x", y="F(x)", title=title)
  if (discrete) p + geom_step(direction="vh") else p + geom_line()
}

# UI -------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Distributions Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Distribution", choices = c(
        "Discrete Uniform" = "dunif",
        "Bernoulli" = "bern",
        "Binomial" = "binom",
        "Poisson" = "pois",
        "Geometric (shifted, support 1,2,...)" = "geom1",
        "Negative Binomial (failures until r-th success)" = "nbinom",
        "Continuous Uniform" = "cunif",
        "Exponential" = "exp",
        "Normal" = "norm",
        "Student's t" = "t"
      )),

      # Parameter inputs (conditionally shown)
      conditionalPanel("input.dist == 'dunif'",
        numericInput("m_disc", "m (min integer)", 1, step=1),
        numericInput("n_disc", "n (max integer)", 7, step=1)
      ),
      conditionalPanel("input.dist == 'bern'",
        sliderInput("p_bern", "p", min=0, max=1, value=0.2, step=0.001)
      ),
      conditionalPanel("input.dist == 'binom'",
        numericInput("n_binom", "n (trials)", 50, min=1, step=1),
        sliderInput("p_binom", "p", min=0, max=1, value=0.98, step=0.001)
      ),
      conditionalPanel("input.dist == 'pois'",
        numericInput("lambda_pois", "lambda (rate)", 1.05, min=0.0001, step=0.01)
      ),
      conditionalPanel("input.dist == 'geom1'",
        sliderInput("p_geom", "p (success prob)", min=0.0001, max=0.9999, value=0.8, step=0.0001)
      ),
      conditionalPanel("input.dist == 'nbinom'",
        numericInput("r_nbinom", "r (successes)", 5, min=1, step=1),
        sliderInput("p_nbinom", "p (success prob)", min=0.0001, max=0.9999, value=0.5, step=0.0001)
      ),
      conditionalPanel("input.dist == 'cunif'",
        numericInput("a_cunif", "a", 1, step=0.1),
        numericInput("b_cunif", "b", 5, step=0.1)
      ),
      conditionalPanel("input.dist == 'exp'",
        numericInput("lambda_exp", "lambda (rate)", 2, min=0.0001, step=0.1)
      ),
      conditionalPanel("input.dist == 'norm'",
        numericInput("mu_norm", "mu", 7, step=0.1),
        numericInput("sigma_norm", "sigma", 10, min=0.0001, step=0.1)
      ),
      conditionalPanel("input.dist == 't'",
        numericInput("df_t", "df (nu)", 19, min=1, step=1)
      ),

      # Common sampling & inference controls
      hr(),
      numericInput("n_sample", "Sample size n", 100, min=1, step=1),
      sliderInput("alpha", "alpha (for CIs / tests)", min=0.001, max=0.2, value=0.05, step=0.001),

      # Null parameters for tests (contextual)
      conditionalPanel("input.dist == 'bern'",
        numericInput("p0", "H0: p =", 0.3, min=0, max=1, step=0.001)
      ),
      conditionalPanel("input.dist == 'norm'",
        numericInput("mu0", "H0: mu =", 7, step=0.1)
      ),

      actionButton("resample", "Resample")
    ),

    mainPanel(
      tabsetPanel(id="tabs",
        tabPanel("Population",
          fluidRow(
            column(6, plotOutput("plot_pdf")),
            column(6, plotOutput("plot_cdf"))
          ),
          verbatimTextOutput("pop_stats")
        ),
        tabPanel("Sample & Estimates",
          verbatimTextOutput("sample_stats")
        ),
        tabPanel("Likelihood (grid)",
          plotOutput("plot_like"),
          verbatimTextOutput("like_note")
        ),
        tabPanel("Confidence Interval",
          verbatimTextOutput("ci_text")
        ),
        tabPanel("Hypothesis Test",
          verbatimTextOutput("ht_text")
        )
      )
    )
  )
)

# Server ----------------------------------------------------------------
server <- function(input, output, session){

  # Population objects ---------------------------------------------------
  pop <- reactive({
    switch(input$dist,
      dunif = {
        m <- as.integer(input$m_disc); n <- as.integer(input$n_disc)
        if (n < m) { tmp <- m; m <- n; n <- tmp }
        x <- m:n
        d <- rep(1/length(x), length(x))
        p <- (x - m + 1)/length(x)
        list(df=tibble(x, d=d, p=p),
             mean = (n+m)/2,
             var = ((n-m)*(n-m+2))/12,
             discrete=TRUE,
             label = sprintf("Discrete Uniform[%d,%d]", m, n))
      },
      bern = {
        p <- input$p_bern
        x <- 0:1
        d <- dbinom(x, size=1, prob=p)
        P <- c(1-p, 1)
        list(df=tibble(x, d=d, p=P),
             mean = p,
             var = p*(1-p),
             discrete=TRUE,
             label = sprintf("Bernoulli(p=%.3f)", p))
      },
      binom = {
        n <- as.integer(input$n_binom); p <- input$p_binom
        x <- 0:n
        d <- dbinom(x, n, p)
        P <- pbinom(x, n, p)
        list(df=tibble(x, d=d, p=P),
             mean = n*p,
             var = n*p*(1-p),
             discrete=TRUE,
             label = sprintf("Binomial(n=%d, p=%.3f)", n, p))
      },
      pois = {
        lambda <- input$lambda_pois
        x <- 0:max(14, ceiling(lambda*6))
        d <- dpois(x, lambda)
        P <- ppois(x, lambda)
        list(df=tibble(x, d=d, p=P),
             mean = lambda,
             var = lambda,
             discrete=TRUE,
             label = sprintf("Poisson(%.3f)", lambda))
      },
      geom1 = {
        p <- input$p_geom
        x <- 1:ceiling(max(10, 8/p))
        d <- dgeom(x-1, p) # shifted support 1,2,...
        P <- pgeom(x-1, p)
        list(df=tibble(x, d=d, p=P),
             mean = 1/p,
             var = (1-p)/p^2,
             discrete=TRUE,
             label = sprintf("Geom-shifted(p=%.3f)", p))
      },
      nbinom = {
        r <- as.integer(input$r_nbinom); p <- input$p_nbinom
        x <- 0:max(50, ceiling(r*(1-p)/p)*3)
        d <- dnbinom(x, r, p)
        P <- pnbinom(x, r, p)
        list(df=tibble(x, d=d, p=P),
             mean = r*(1-p)/p,
             var = r*(1-p)/p^2,
             discrete=TRUE,
             label = sprintf("NegBin(r=%d, p=%.3f)", r, p))
      },
      cunif = {
        a <- input$a_cunif; b <- input$b_cunif
        if (b < a){ tmp <- a; a <- b; b <- tmp }
        x <- seq(a, b, length.out=201)
        d <- dunif(x, a, b)
        P <- punif(x, a, b)
        list(df=tibble(x, d=d, p=P),
             mean = (a+b)/2,
             var = (b-a)^2/12,
             discrete=FALSE,
             label = sprintf("Uniform(%.3f, %.3f)", a, b))
      },
      exp = {
        lambda <- input$lambda_exp
        x <- seq(0, qexp(0.995, lambda), length.out=500)
        d <- dexp(x, lambda)
        P <- pexp(x, lambda)
        list(df=tibble(x, d=d, p=P),
             mean = 1/lambda,
             var = 1/lambda^2,
             discrete=FALSE,
             label = sprintf("Exponential(%.3f)", lambda))
      },
      norm = {
        mu <- input$mu_norm; s <- input$sigma_norm
        x <- seq(mu - 4*s, mu + 4*s, length.out=401)
        d <- dnorm(x, mu, s)
        P <- pnorm(x, mu, s)
        list(df=tibble(x, d=d, p=P),
             mean = mu,
             var = s^2,
             discrete=FALSE,
             label = sprintf("Normal(%.3f, %.3f)", mu, s))
      },
      t = {
        df <- as.integer(input$df_t)
        x <- seq(-4, 4, length.out=401)
        d <- dt(x, df)
        P <- pt(x, df)
        list(df=tibble(x, d=d, p=P),
             mean = if (df>1) 0 else NA_real_,
             var = if (df>2) df/(df-2) else NA_real_,
             discrete=FALSE,
             label = sprintf("Student's t(df=%d)", df))
      }
    )
  })

  output$plot_pdf <- renderPlot({
    pp <- pop()
    plot_pmf_pdf(pp$df, discrete=pp$discrete, title=paste0(pp$label, ": PMF/PDF"))
  })

  output$plot_cdf <- renderPlot({
    pp <- pop()
    plot_cdf(pp$df, discrete=pp$discrete, title=paste0(pp$label, ": CDF"))
  })

  output$pop_stats <- renderText({
    pp <- pop()
    paste0("Mean = ", num_fmt(pp$mean), ", Variance = ", num_fmt(pp$var))
  })

  # Sampling -------------------------------------------------------------
  sample_data <- reactiveVal(NULL)

  observeEvent(input$resample, {
    set.seed(NULL)
    n <- as.integer(input$n_sample)
    s <- switch(input$dist,
      dunif = {
        m <- as.integer(input$m_disc); nmax <- as.integer(input$n_disc)
        if (nmax < m) { tmp <- m; m <- nmax; nmax <- tmp }
        sample(m:nmax, n, replace=TRUE)
      },
      bern = rbinom(n, 1, input$p_bern),
      binom = rbinom(n, size=1, prob=input$p_binom), # one Bernoulli draw per row; sample mean ~ p
      pois = rpois(n, input$lambda_pois),
      geom1 = rgeom(n, input$p_geom) + 1,
      nbinom = rnbinom(n, size=as.integer(input$r_nbinom), prob=input$p_nbinom),
      cunif = runif(n, input$a_cunif, input$b_cunif),
      exp = rexp(n, input$lambda_exp),
      norm = rnorm(n, input$mu_norm, input$sigma_norm),
      t = rt(n, df=as.integer(input$df_t))
    )
    sample_data(s)
  }, ignoreInit=TRUE)

  output$sample_stats <- renderText({
    x <- sample_data()
    validate(need(!is.null(x), "Click Resample to generate a sample."))

    n <- length(x)
    xbar <- mean(x)
    s2 <- var(x)
    s <- sd(x)

    txt <- paste0(
      "n = ", n, "\n",
      "sample_mean = ", num_fmt(xbar), "\n",
      "sample_var = ", num_fmt(s2), "\n",
      "sample_sd = ", num_fmt(s), "\n"
    )

    # Distribution-specific estimated var of Xbar
    if (input$dist %in% c("bern","binom")){
      est_var_xbar <- (xbar*(1-xbar))/n
      txt <- paste0(txt, "est_var_X_bar = ", num_fmt(est_var_xbar), ", est_sd_X_bar = ", num_fmt(sqrt(est_var_xbar)), "\n")
    } else if (input$dist %in% c("pois")){
      est_var_xbar <- xbar/n
      txt <- paste0(txt, "est_var_X_bar = ", num_fmt(est_var_xbar), ", est_sd_X_bar = ", num_fmt(sqrt(est_var_xbar)), "\n")
    } else if (input$dist %in% c("norm","cunif","exp","geom1","nbinom","dunif","t")){
      est_var_xbar <- s2/n
      txt <- paste0(txt, "est_var_X_bar (using s^2) = ", num_fmt(est_var_xbar), ", est_sd_X_bar = ", num_fmt(sqrt(est_var_xbar)), "\n")
    }

    txt
  })

  # Likelihood (grid over parameter) ------------------------------------
  output$plot_like <- renderPlot({
    x <- sample_data(); validate(need(!is.null(x), "Click Resample first."))
    n <- length(x)

    if (input$dist == "bern" || input$dist == "binom"){
      theta <- seq(0, 1, by=0.001)
      k <- sum(x)  # treat as single binomial observation
      L <- dbinom(k, size=n, prob=theta)
      df <- tibble(theta, L)
      ggplot(df, aes(theta, L)) + geom_line() + labs(title="Likelihood for p", x="p", y="L(p)")

    } else if (input$dist == "pois"){
      theta <- seq(0, max(1e-6, max(x)), length.out=1000)
      L <- exp(-n*theta) * theta^(n*mean(x)) / prod(factorial(x))
      df <- tibble(theta, L)
      ggplot(df, aes(theta, L)) + geom_line() + labs(title="Likelihood for lambda", x="lambda", y="L(lambda)")

    } else if (input$dist == "exp"){
      theta <- seq(0, 3*max(1e-6, 1/mean(x)), length.out=1000)
      L <- theta^n * exp(-theta*n*mean(x))
      df <- tibble(theta, L)
      ggplot(df, aes(theta, L)) + geom_line() + labs(title="Likelihood for lambda", x="lambda", y="L(lambda)")

    } else if (input$dist == "geom1"){
      theta <- seq(0.0001, 0.9999, by=0.001)
      # shifted geometric (support 1,2,...) => P(X=x) = (1-p)^(x-1)p
      L <- (1-theta)^(sum(x)-n) * theta^n
      df <- tibble(theta, L)
      ggplot(df, aes(theta, L)) + geom_line() + labs(title="Likelihood for p (shifted geometric)", x="p", y="L(p)")

    } else {
      ggplot() + annotate("text", x=0.5, y=0.5, label="Likelihood not implemented for this distribution", size=6) +
        theme_void()
    }
  })

  output$like_note <- renderText({
    if (input$dist %in% c("bern","binom","pois","exp","geom1")) return("")
    "Note: likelihood plotting is implemented for Bernoulli/Binomial, Poisson, Exponential, and shifted Geometric."
  })

  # Confidence Intervals -------------------------------------------------
  output$ci_text <- renderText({
    x <- sample_data(); validate(need(!is.null(x), "Click Resample first."))
    n <- length(x); alpha <- input$alpha
    xbar <- mean(x); s <- sd(x)

    if (input$dist %in% c("bern","binom")){
      se <- sqrt(xbar*(1-xbar)/n)
      z <- qnorm(1-alpha/2)
      CI <- c(xbar - z*se, xbar + z*se)
      paste0("Proportion CI (normal approx): [", num_fmt(CI[1]), ", ", num_fmt(CI[2]), "]")

    } else if (input$dist == "pois"){
      se <- sqrt(xbar/n)
      z <- qnorm(1-alpha/2)
      CI <- c(xbar - z*se, xbar + z*se)
      paste0("Mean (lambda) CI (normal approx): [", num_fmt(CI[1]), ", ", num_fmt(CI[2]), "]")

    } else if (input$dist %in% c("norm","cunif","exp","geom1","nbinom","dunif","t")){
      # t-interval for mean using sample sd
      se <- s/sqrt(n)
      tcrit <- qt(1-alpha/2, df=n-1)
      CI <- c(xbar - tcrit*se, xbar + tcrit*se)
      paste0("Mean CI (t): [", num_fmt(CI[1]), ", ", num_fmt(CI[2]), "]")
    }
  })

  # Hypothesis tests -----------------------------------------------------
  output$ht_text <- renderText({
    x <- sample_data(); validate(need(!is.null(x), "Click Resample first."))
    n <- length(x); alpha <- input$alpha
    xbar <- mean(x); s <- sd(x)

    if (input$dist == "bern"){
      p0 <- input$p0
      z <- (xbar - p0)/sqrt(p0*(1-p0)/n)
      pval <- 2*(1-pnorm(abs(z)))
      crit <- qnorm(1-alpha/2)
      paste0(
        "H0: p = ", p0, "\n",
        "z = ", num_fmt(z), ", two-sided p-value = ", num_fmt(pval),
        "\nReject H0? ", ifelse(abs(z) > crit, "Yes", "No")
      )

    } else if (input$dist == "norm"){
      mu0 <- input$mu0
      se <- s/sqrt(n)
      z <- (xbar - mu0)/se
      pval <- 2*(1-pnorm(abs(z)))
      crit <- qnorm(1-alpha/2)
      paste0(
        "H0: mu = ", mu0, "\n",
        "z = ", num_fmt(z), ", two-sided p-value = ", num_fmt(pval),
        "\nReject H0? ", ifelse(abs(z) > crit, "Yes", "No"),
        "\n(Note: uses z with s; for small n consider t-test.)"
      )

    } else {
      # Generic t-test for mean against 0
      mu0 <- 0
      se <- s/sqrt(n)
      tstat <- (xbar - mu0)/se
      pval <- 2*(1-pt(abs(tstat), df=n-1))
      crit <- qt(1-alpha/2, df=n-1)
      paste0(
        "H0: mu = 0\n",
        "t = ", num_fmt(tstat), ", df = ", n-1, ", two-sided p-value = ", num_fmt(pval),
        "\nReject H0? ", ifelse(abs(tstat) > crit, "Yes", "No")
      )
    }
  })
}

shinyApp(ui, server)
