
library(shiny)
library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(PowerTOST)

# Copied from https://bebac.at/articles/Power-Table.phtml
power.table <- function(alpha = 0.05, CV = CV, theta0, theta1, theta2,
                        design, ABE = TRUE, method = "exact", NTID = FALSE,
                        target = 0.9, minpower = 0.7, bal = TRUE, do.rate = 0,
                        regulator = "EMA", nsims = 1e5, TIE = FALSE,
                        details = TRUE) {
  ###########################################################################
  # Arguments:                                                              #
  # alpha    : Nominal level of the test (default 0.05)                     #
  # CV       : within-subject (crossovers, paired design), total (parallel) #
  #            can be a two-element vector for SABE; first element CVw of T #
  #            and second element CVw of R                                  #
  # theta0   : assumed T/R-ratio (default 0.95 for ABE, 0.90 for SABE,      #
  #            0.975 for NTID and ABE if theta1 = 0.9 or theta1 = 1 / 0.9   #
  # theta1   : lower BE-limit (ABE/NTID), lower PE-constraint (SABE)        #
  #            default 0.80 for ABE/NTID, fixed PE-constraint 0.80 for SABE #
  # theta2   : upper BE-limit (ABE/NTID), upper PE-constraint (SABE)        #
  #            default 1.25 for ABE/NTID, fixed PE-constraint 1.25 for SABE #
  # design   : any supported; see known.designs()                           #
  #            defaults: "2x2x2" for ABE, "2x3x3" for SABE (ABEL and RSABE) #
  #                      "2x2x4" for NTIDs                                  #
  # ABE      : TRUE (default) for ABE, FALSE for SABE                       #
  # method   : ABE only; power method (default "exact")                     #
  # NTID     : FALSE (default), TRUE for RSABE (FDA, CDE)                   #
  # target   : target (desired) power (default 0.8)                         #
  # minpower : minimum acceptable power (default 0.7)                       #
  #            analog to pa.ABE(), pa.scABE(), pa.NTID()                    #
  # do.rate  : anticipated dropout rate (default 0); if > 0, the estimated  #
  #            sample size will be adjusted in order to have at least the   #
  #            estimated sample size in the worst case                      #
  # bal      : should only balanced sequences (crossovers, replicates) or   #
  #            equal group sizes (parallel) be shown (default TRUE)?        #
  # nsims    : SABE only; number of simulations for SABE (default 1e5)      #
  # regulator: SABE only; any of "EMA", "HC", "GCC", "FDA" (default "EMA")  #
  # TIE      : SABE only; empiric Type I Error (default FALSE)              #
  # details  : should the aggregated results and data frame be printed      #
  #            (default TRUE)?                                              #
  # ----------------------------------------------------------------------- #
  # Notes                                                                   #
  # - Only for the multiplicative model, i.e., in PowerTOST’s functions for #
  #   ABE logscale = TRUE is applied                                        #
  # - If simulating the empiric Type I Error 1e6 simulations are employed   #
  # - Sample size adjustments are tricky; I hope to have covered them       #
  #   correctly:                                                            #
  #   - In general 12 _eligible_ subjects are required                      #
  #   - The EMA requires for ABEL 12 _eligible_ subjects in a 2-sequence 3- #
  #     period full replicate design in the RR-sequence, i.e., 24 total     #
  #   - In RSABE (FDA and CDE) we need 24 _dosed_ subjects                  #
  #   - there is no minimum sample size for RSABE of NTIDS given in the     #
  #     guidance; I assumed 12 _eligible_ subjects                          #
  ###########################################################################
  require(PowerTOST)
  balance <- function(n, n.seq) {
    # Round up the sample size to obatain equally sized sequences (crossover /
    # replicate designs) or groups (parallel design)
    return(as.integer(n.seq * (n %/% n.seq + as.logical(n %% n.seq))))
  }
  nadj <- function(n, do.rate, n.seq) {
    # Adjust the sample size for the anticipated dropout-rate
    return(as.integer(balance(n / (1 - do.rate), n.seq)))
  }
  # Input checking: Error handling, dealing with missing arguments,
  #                 messages about strange input
  CV <- abs(CV)
  if (NTID) ABE <- FALSE # we apply RSABE
  if (ABE & length(CV) > 1)
    stop("Give \"CV\" as a one element vector for ABE.")
  if (!ABE & length(CV) > 2)
    stop("\"CV\" must have not more than two elements for SABE.")
  if (!ABE & !regulator %in% c("EMA", "FDA", "HC", "GCC"))
    stop("regulator \"", regulator, "\" not supported.")
  if (missing(design)) {
    if (ABE)          design <- "2x2x2"
    if (!ABE & !NTID) design <- "2x3x3"
    if (NTID)         design <- "2x2x4"
  }
  if (!ABE & !design %in% c("2x3x3", "2x2x4", "2x2x3"))
    stop("A replicate design is required for reference-scaling.")
  if (NTID & !design %in% c("2x2x4", "2x2x3"))
    stop("A full replicate design is required for reference-scaling.")
  if (missing(theta1) & missing(theta2)) theta1 <- 0.80
  if (missing(theta1)) theta1 <- 1/theta2
  if (missing(theta2)) theta2 <- 1/theta1
  if (!method %in% c("exact", "owenq", "mvt", "noncentral",
                     "nct", "central", "shifted"))
    stop("method \"", method, "\" not supported.")
  if (target <= 0.5) stop("\"target = ", target, "\" does not make sense.")
  if (minpower <= 0.5) stop("\"minpower = ", minpower,
                            "\" does not make sense.")
  if (minpower >= target) stop("\"minpower\" must be < \"target\".")
  if (nsims < 5e4) message("Note: sample size based on < 50,000",
                           "\n      simulations may be unreliable.")
  do.rate <- abs(do.rate)
  if (do.rate > 0 & do.rate >= 0.5) stop("\"do.rate = ", do.rate,
                                         "\" does not make sense.")
  if (NTID) {
    # Not a good idea to fiddle with the limits!
    if (!theta1 == 0.8) {
      message(sprintf("Note: theta1 = %.4f ", theta1),
              "\n      not compliant with the guidance.")
    }
    if (!theta2 == 1.25) {
      message(sprintf("Note: theta2 = %.4f ", theta2),
              "\n      not compliant with the guidance.")
    }
  }
  # Default theta0 for the various approaches
  if (missing(theta0)) {
    if (NTID | theta1 == 0.9 | theta2 == 1 / 0.9) {
      theta0 <- 0.975
    }else {
      ifelse (ABE, theta0 <- 0.95, theta0 <- 0.9)
    }
  }
  # End of input checking
  if (!ABE) {# for SABE and NTID
    if (length(CV) == 2) {
      CVwT <- CV[1]
      CVwR <- CV[2]
    }else {
      CVwT <- CVwR <- CV
    }
  }
  # Minimum sample sizes acc. to GLs
  if (ABE | NTID) {
    min.n <- 12                          # eligible
  }else {
    if (regulator %in% c("HC", "GCC")) {
      min.n <- 12                        # eligible
    }else {
      if (regulator == "FDA") {
        min.n <- 24                      # dosed
      }else {# EMA
        if (!design == "2x2x3") {
          min.n <- 12                    # eligible
        }else {
          min.n <- 24                    # eligible
        }
      }
    }
  }
  adj     <- FALSE # default if no dropouts or sample size compliant with GLs
  designs <- known.designs()[, c(1:2, 5)]
  n.step  <- n.seq <- designs[designs$design == design, "steps"]
  if (!bal) n.step <- 1L
  if (NTID) {
    tmp <- sampleN.NTID(alpha = alpha, CV = CV, theta0 = theta0,
                        theta1 = theta1, theta2 = theta2, targetpower = target,
                        design = design, nsims = nsims, details = FALSE,
                        print = FALSE)
    n   <- n.orig <- tmp[["Sample size"]]
    pwr.orig      <- tmp[["Achieved power"]]
    if (n < min.n) {# 12; nothing in the guidance - or is it 24 like in RSABE?
      n   <- min.n
      adj <- TRUE
    }
    if (do.rate > 0) {
      n.est <- n
      n     <- nadj(n, do.rate, n.seq)
      adj   <- TRUE
    }
    n.lo <- sampleN.NTID(alpha = alpha, CV = CV, theta0 = theta0,
                         theta1 = theta1, theta2 = theta2,
                         targetpower = minpower, design = design, nsims = nsims,
                         details = FALSE, print = FALSE)[["Sample size"]]
  }else {
    if (ABE) {
      tmp <- sampleN.TOST(alpha = alpha, CV = CV, theta0 = theta0,
                          theta1 = theta1, theta2 = theta2,
                          targetpower = target, design = design,
                          method = method, print = FALSE)
      n   <- n.orig <- tmp[["Sample size"]]
      pwr.orig      <- tmp[["Achieved power"]]
      if (n < min.n) {
        n   <- min.n
        adj <- TRUE
      }
      if (do.rate > 0) {
        n.est <- n
        n     <- nadj(n, do.rate, n.seq)
        adj   <- TRUE
      }
      if (bal) {
        n.lo <- sampleN.TOST(alpha = alpha, CV = CV, theta0 = theta0,
                             theta1 = theta1, theta2 = theta2,
                             targetpower = minpower, design = design,
                             method = method, print = FALSE)[["Sample size"]]
      }else {
        # this section lifted/adapted from PowerTOST / pa.ABE.R
        n.est <- n
        Ns    <- seq(n.est, 12)
        if (n.est == 12) Ns <- seq(n.est, 2 * n.step)
        nNs   <-length(Ns)
        pwrN  <- power.TOST(alpha = alpha, CV = CV, theta0 = theta0,
                            theta1 = theta1, theta2 = theta2,
                            design = design, method = method, n = n)
        n.min <- NULL
        pBEn  <- NULL
        n     <- numeric(n.step)
        ni    <- 1:n.step
        i     <- 0
        while (pwrN >= minpower & i < nNs) {
          i          <- i + 1
          n[-n.step] <- diff(floor(Ns[i] * ni / n.step))
          n[n.step]  <- Ns[i] - sum(n[-n.step])
          pwrN       <- suppressMessages(
            power.TOST(alpha = alpha, CV = CV, theta0 = theta0,
                       theta1 = theta1, theta2 = theta2,
                       design = design, method = method, n = n))
          if (pwrN >= minpower) {
            n.min <- c(n.min, sum(n))
            pBEn <- c(pBEn, pwrN)
          }else {
            break
          }
        }
        n    <- max(n.min)
        n.lo <- min(n.min)
      }
    }else {
      if (regulator == "FDA") {
        tmp <- sampleN.RSABE(alpha = alpha, CV = CV, theta0 = theta0,
                             targetpower = target, design = design,
                             nsims = nsims, details = FALSE, print = FALSE)
        n   <- n.orig <- tmp[["Sample size"]]
        pwr.orig      <- tmp[["Achieved power"]]
        if (n < min.n) {
          n   <- min.n
          adj <- TRUE
        }
        if (do.rate > 0) {
          n.est <- n
          n     <- nadj(n, do.rate, n.seq)
          adj   <- TRUE
        }
        if (bal) {
          n.lo <- sampleN.RSABE(alpha = alpha, CV = CV, theta0 = theta0,
                                targetpower = minpower, design = design,
                                nsims = nsims, details = FALSE,
                                print = FALSE)[["Sample size"]]
        }else {
          # this section lifted/adapted from PowerTOST / pa.scABE.R
          n.est <- n
          Ns    <- seq(n.est, 12)
          if (n.est == 12) Ns <- seq(n.est, 2 * n.step)
          nNs   <-length(Ns)
          pwrN  <- power.RSABE(alpha = alpha, CV = CV, theta0 = theta0,
                               design = design, n = n, nsims = nsims)
          n.min <- NULL
          pBEn  <- NULL
          n     <- numeric(n.step)
          ni    <- 1:n.step
          i     <- 0
          while(pwrN >= minpower & i < nNs){
            i    <- i + 1
            pwrN <- suppressMessages(
              power.RSABE(alpha = alpha, CV = CV, theta0 = theta0,
                          design = design, n = Ns[i], nsims = nsims))
            if (pwrN >= minpower) {
              n.min <- c(n.min, Ns[i])
              pBEn  <- c(pBEn, pwrN)
            }else {
              break
            }
          }
          n    <- max(n.min)
          n.lo <- min(n.min)
        }
      }else {
        if (!regulator == "HC" & design == "2x3x3" & CVwT > CVwR) {
          # Subject simulations if partial replicate and heteroscedasticity,
          # where T is more variable than R
          cat("Patience please. Simulating subject data is time consuming.\n")
          flush.console()
          tmp <- sampleN.scABEL.sdsims(alpha = alpha, CV = CV, theta0 = theta0,
                                       targetpower = target, design = design,
                                       regulator = regulator, nsims = nsims,
                                       details = FALSE,
                                       print = FALSE)
          n   <- n.orig <- tmp[["Sample size"]]
          pwr.orig      <- tmp[["Achieved power"]]
          n.lo <- sampleN.scABEL.sdsims(alpha = alpha, CV = CV, theta0 = theta0,
                                        targetpower = minpower, design = design,
                                        regulator = regulator, nsims = nsims,
                                        details = FALSE,
                                        print = FALSE)[["Sample size"]]
        }else {
          tmp <- suppressWarnings(
            sampleN.scABEL(alpha = alpha, CV = CV, theta0 = theta0,
                           targetpower = target, design = design,
                           regulator = regulator, nsims = nsims,
                           details = FALSE,
                           print = FALSE))
          n   <- n.orig <- tmp[["Sample size"]]
          pwr.orig      <- tmp[["Achieved power"]]
          if (n < min.n) {
            n   <- min.n
            adj <- TRUE
          }
          if (do.rate > 0) {
            n.est <- n
            n     <- nadj(n, do.rate, n.seq)
            adj   <- TRUE
          }
          if (bal) {
            n.lo <- sampleN.scABEL(alpha = alpha, CV = CV, theta0 = theta0,
                                   targetpower = minpower, design = design,
                                   regulator = regulator, nsims = nsims,
                                   details = FALSE,
                                   print = FALSE)[["Sample size"]]
          }else {
            # this section lifted/adapted from PowerTOST / pa.scABE.R
            n.est <- n
            Ns    <- seq(n.est, 12)
            if (n.est == 12) Ns <- seq(n.est, 2 * n.step)
            nNs   <-length(Ns)
            pwrN  <- power.scABEL(alpha = alpha, CV = CV, theta0 = theta0,
                                  design = design, regulator = regulator, n = n,
                                  nsims = nsims)
            n.min <- NULL
            pBEn  <- NULL
            n     <- numeric(n.step)
            ni    <- 1:n.step
            i     <- 0
            while(pwrN >= minpower & i < nNs){
              i    <- i + 1
              pwrN <- suppressMessages(
                power.scABEL(alpha = alpha, CV = CV, theta0 = theta0,
                             design = design, regulator = regulator,
                             n = Ns[i], nsims = nsims))
              if (pwrN >= minpower) {
                n.min <- c(n.min, Ns[i])
                pBEn  <- c(pBEn, pwrN)
              }else {
                break
              }
            }
            n    <- max(n.min)
            n.lo <- min(n.min)
          }
        }
      }
    }
  }
  x <- data.frame(n = seq(n, n.lo, -n.step), do = NA_integer_,
                  do.pct = NA_real_, power = NA_real_)
  if (NTID) x <- cbind(x, pwr.RSABE = NA_real_, pwr.ABE = NA_real_,
                       pwr.sratio = NA_real_)
  if (!ABE & TIE) {
    x        <- cbind(x, TIE = NA_real_, infl = "no ")
    infl.lim <- binom.test(alpha * 1e6, 1e6, alternative = "less")$conf.int[2]
  }
  if (TIE) {
    cat("Patience please. Simulating empiric Type I Error is time consuming.\n")
    flush.console()
  }
  for (i in 1:nrow(x)) {
    if (NTID) {
      x[i, 4:7] <- suppressMessages(
        as.numeric(
          power.NTID(alpha = alpha, CV = CV, theta0 = theta0,
                     theta1 = theta1, theta2 = theta2,
                     design = design, n = x$n[i], nsims = nsims,
                     details = TRUE)[1:4]))
    }else {
      if (ABE) {
        x$power[i] <- suppressMessages(
          power.TOST(alpha = alpha, CV = CV, theta0 = theta0,
                     theta1 = theta1, theta2 = theta2,
                     design = design, method = method,
                     n = x$n[i]))
      }else {
        reg    <- reg_const(regulator = regulator)
        limits <- scABEL(CV = CVwR, regulator = regulator)
        if (regulator == "FDA") {
          x$power[i] <- suppressMessages(
            power.RSABE(alpha = alpha, CV = CV, theta0 = theta0,
                        design = design, n = x$n[i],
                        nsims = nsims))
          if (TIE) {
            x$TIE[i] <- suppressMessages(
              power.RSABE(alpha = alpha, CV = CV,
                          theta0 = limits[["upper"]],
                          design = design, n = x$n[i], nsims = 1e6))
            if (x$TIE[i] > infl.lim) x$infl[i] <- "yes"
          }
        }else {
          if (!regulator == "HC" & design == "2x3x3" & CVwT > CVwR) {
            x$power[i] <- suppressMessages(
              power.scABEL.sdsims(alpha = alpha, CV = CV,
                                  theta0 = theta0,
                                  design = design,
                                  regulator = regulator,
                                  n = x$n[i], nsims = nsims))
            if (TIE) {
              x$TIE[i] <- suppressMessages(
                power.scABEL.sdsims(alpha = alpha, CV = CV,
                                    theta0 = limits[["upper"]],
                                    design = design,
                                    regulator = regulator,
                                    n = x$n[i], nsims = 1e6))
              if (x$TIE[i] > infl.lim) x$infl[i] <- "yes"
            }
          }else {
            x$power[i] <- suppressMessages(
              suppressWarnings(
                power.scABEL(alpha = alpha, CV = CV,
                             theta0 = theta0,
                             design = design,
                             regulator = regulator,
                             n = x$n[i], nsims = nsims)))
            if (TIE) {
              x$TIE[i] <- suppressMessages(
                suppressWarnings(
                  power.scABEL(alpha = alpha, CV = CV,
                               theta0 = limits[["upper"]],
                               design = design,
                               regulator = regulator, n = x$n[i],
                               nsims = 1e6)))
              if (x$TIE[i] > infl.lim) x$infl[i] <- "yes"
            }
          }
        }
      }
    }
  }
  # dropouts (number and in percent)
  x$do     <- max(x$n) - x$n
  x$do.pct <- 100 * x$do / max(x$n)
  # prepare for the output
  desis <- c("two parallel groups", rep("conventional crossover", 2),
             "Higher-Order crossover (Latin Squares)",
             "Higher-order crossover (Williams\u2019)",
             "Higher-Order crossover (Latin Squares or Williams\u2019)",
             "2-sequence 3-period full replicate",
             "2-sequence 4-period full replicate",
             "4-sequence 4-period full replicate",
             "3-sequence 3-period \u2018partial\u2019 replicate",
             "4-sequence 2-period full replicate; Balaam\u2019s",
             "Liu\u2019s 2×2×2 repeated crossover",
             "paired means")
  regs <- data.frame(regulator = c("EMA", "HC", "GCC", "FDA"),
                     reg.name = c("EMA and others", "Health Canada",
                                  "Gulf Cooperation Council",
                                  "FDA or CDE"))
  # cosmetics
  x$do.pct       <- sprintf("%.3f%%", x$do.pct)
  names(x)[2:3]  <- rep("dropouts", 2)
  x[1, 2:3]      <- rep(0, 2) # !! Changed by me
  ifelse (!TIE, x[, 4:ncol(x)] <- signif(x[, 4:ncol(x)], 5),
          x[, 4:(ncol(x)-1)] <- signif(x[, 4:(ncol(x)-1)], 5))
  if (NTID) names(x)[5:7] <- c("RSABE", "ABE", "s-ratio")
  f1   <- "\n%-14s: %s"
  f2   <- "\n%-14s: %.4f"
  f3   <- "\n%-14s: %.7g"
  f4   <- "\n%-14s: %.4f, %.4f"
  f5   <- "\n%-14s: %4.1f%% (anticipated)"
  f6   <- "\n%-14s: %-3s   (achieved power %.5f)"
  txt  <- paste0("design        : ", design,
                 " (", desis[which(design == designs$design)], ")")
  txt <- paste(txt, sprintf(f2, "alpha", alpha))
  if (ABE) {
    txt <- paste(txt, sprintf(f1, "method", "Average Bioequivalence (ABE)"))
    if (!theta1 == 0.9) {
      txt <- paste(txt, sprintf(f1, "regulator", "all jurisdictions"))
    }else {
      txt <- paste0(txt, "; NTID")
      txt <- paste(txt, sprintf(f1, "regulator", "EMA and others"))

    }
    txt <- paste(txt, sprintf(f1, "power method", method))
  }else {
    if (NTID) {
      r_const <- log(1 / 0.9) / 0.1
      # Note: That’s the exact value; in the guidance log(1.1111) / 0.1
      txt <- paste(txt, sprintf(f1, "regulator", "FDA or CDE"))
      txt <- paste(txt, sprintf(f1, "method", "RSABE (NTID)"))
      txt <- paste(txt, sprintf(f1, "simulations", "intra-subject contrasts"))
      txt <- paste(txt, sprintf("\n%-14s:", "  number"),
                   formatC(nsims, format = "d", big.mark = ","))
      txt <- paste(txt, sprintf(f3, "Reg. constant", r_const),
                   sprintf(f2, "\u2018Cap\u2019         ",
                           se2CV(log(theta2) / r_const)))
      txt <- paste(txt, sprintf(f4, "Implied limits",
                                exp(-r_const * CV2se(CVwR)),
                                exp(r_const * CV2se(CVwR))))
      txt <- paste(txt, sprintf(f4, "ABE limits", theta1, theta2))
    }else {# SABE
      txt <- paste(txt, sprintf(f1, "regulator",
                                regs$reg.name[regs$regulator == regulator]))
      if (regulator == "FDA") {
        txt <- paste(txt, sprintf(f1, "method",
                                  "Reference-scaled Average Bioequivalence (RSABE)"))
      }else {
        if (!regulator == "GCC") {
          txt <- paste(txt, sprintf(f1, "method",
                                    "Average Bioequivalence with Expanding Limits (ABEL)"))
        }else {
          txt <- paste(txt, sprintf(f1, "method",
                                    "Average Bioequivalence with widened limits"))
        }
      }
      if (reg$est_method == "ANOVA") {
        txt <- paste(txt, sprintf(f1, "simulations", "ANOVA"))
      }else {
        txt <- paste(txt, sprintf(f1, "simulations", "intra-subject contrasts"))
      }
      txt <- paste(txt, sprintf("\n%-14s:", "  number"),
                   formatC(nsims, format = "d", big.mark = ","))
      txt <- paste(txt, sprintf(f2, "CVswitch", reg$CVswitch),
                   sprintf(f3, "Reg. constant", reg$r_const))
      if (reg$pe_constr) {
        txt <- paste(txt, sprintf(f4, "PE constraints", theta1, theta2))
      }
      if (!regulator %in% c("FDA", "GCC")) {
        txt <- paste(txt, sprintf(f2, "Cap", reg$CVcap))
      }
    }
  }
  if (ABE) {
    if (design == "parallel") {
      txt <- paste(txt, sprintf(f2, "CVtotal", CV))
    }else {
      txt <- paste(txt, sprintf(f2, "CVw", CV))
    }
  }else {# SABE
    if (length(unique(CV)) == 1) {
      txt <- paste(txt, sprintf(f2, "CVwT = CVwR", unique(CV)))
    }else {
      txt <- paste(txt, sprintf(f2, "CVwT", CVwT))
      txt <- paste(txt, sprintf(f2, "CVwR", CVwR))
    }
  }
  txt <- paste(txt, sprintf(f2, "theta0", theta0))
  if (ABE) {
    txt <- paste(txt, sprintf(f2, "theta1", theta1),
                 sprintf(f2, "theta2", theta2))
  }else {# SABE
    if (!NTID) {# already handled above
      if (regulator == "FDA") {
        if (CV2se(CVwR) < 0.294) {
          txt <- paste(txt, sprintf(f2, "theta1", theta1),
                       sprintf(f2, "theta2", theta2))
        }else {
          txt <- paste(txt, sprintf(f4, "Implied limits",
                                    limits[["lower"]], limits[["upper"]]))
        }
      }
      if (regulator %in% c("EMA", "HC")) {
        if (CVwR <= 0.3) {
          txt <- paste(txt, sprintf(f2, "theta1", theta1),
                       sprintf(f2, "theta2", theta2))
        }else {
          txt <- paste(txt, sprintf(f4, "L, U",
                                    limits[["lower"]], limits[["upper"]]))
        }
      }
      if (regulator == "GCC") {
        if (CVwR <= 0.3) {
          txt <- paste(txt, sprintf(f2, "theta1", theta1),
                       sprintf(f2, "theta2", theta2))
        }else {
          txt <- paste(txt, sprintf(f4, "L, U", 0.75, 1 / 0.75))
        }
      }
    }
  }
  txt <- paste(txt, sprintf(f2, "targetpower", target))
  txt <- paste(txt, sprintf(f6, "estimated n", n.orig, pwr.orig))

  if (adj) {
    if (do.rate == 0) {
      txt <- paste(txt, sprintf(f6, "adjusted n", n, x$power[x$n == n]))
      txt <- paste(txt, "\n               ",
                   "increased to comply with regulatory minimum")
    }else {
      txt <- paste(txt, sprintf(f5, "dropout-rate", 100 * do.rate))
      txt <- paste(txt,  sprintf(f6, "adjusted n", n, x$power[x$n == n]))
    }
  }
  txt <- paste(txt, sprintf(f2, "minpower", minpower))
  txt <- paste(txt, sprintf(f6, "minimum n", min(x$n), min(x$power)), "\n\n")
  if (details) {
    cat(txt)
    print(x, row.names = FALSE)
  }
  # Notes about potential problems with droputs
  if ((ABE | NTID) & min(x$n) <= 12) {
    if (x$dropouts[x$n == 12] == "none") {
      msg <- "Note: If there will be at least one dropout, less than"
    }else {
      msg <- paste0("Note: If there will be more than ",
                    x$dropouts[x$n == 12], " dropouts, less than")
    }
    message(paste(msg, "\n      12 eligible subjects; consider to increase the",
                  "\n      sample size in a pivotal study."))
  }
  if (!ABE & regulator == "EMA" & design == "2x2x3" & head(x$n, 1) <= 24) {
    message("Note: If there will be more than ", x$dropouts[x$n == 24],
            " dropouts, possibly",
            "\n      less than 12 eligible subjects in the RR-sequence;",
            "\n      consider to increase the sample size in a pivotal study.")
  }
  return(invisible(x))
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Calculator for BE stuff"),

    # Sidebar
    sidebarLayout(
        sidebarPanel(
            hr(),
            numericInput(
              "CI_lower",
              "Lower limit of reported Conf.Interval:",
              value = 0.81,
              step = .01
            ),
            numericInput(
              "CI_upper",
              "Upper limit of reported Conf.Interval:",
              value = 1.11,
              step = .01
            ),
            hr(),
            sliderInput("N_sub",
                        "Number of subjects:",
                        min = 2,
                        max = 100,
                        value = 24),
            selectInput("design","Reported design:",
                        c("2x2"="2x2",
                          "2x2x4"="2x2x4",
                          "3x3"="3x3")),


            # Download button
            downloadButton("downloadData", "Download results"),

            # Add your discreet message at the bottom
            tags$hr(),  # Horizontal line for separation

            # Footer with discreet message and custom link
            tags$footer(

              # Github link
              tags$p(
                tags$a(href = "https://github.com/MartynK/BE.Shiny.Calculator/", "Github repo"),
                style = "font-size: 80%; color: grey; text-align: center; padding-top: 10px;"
              ),

              # Author's link
              tags$p(
                tags$a(href = "https://www.linkedin.com/in/marton-kiss-md-29308a112/", "by Márton Kiss"),
                style = "font-size: 80%; color: grey; text-align: center; padding-top: 10px;"
              )
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(

          tags$p(
            tags$p(textOutput("orig_val")),
            style = "font-size: 180%; color: grey; text-align: center; padding-top: 10px;"
          ),



          hr(),
          tags$p(
            tags$p("The below table is for a standard 2x2 design"),
            style = "font-size: 180%; color: grey; text-align: center; padding-top: 10px;"
          ),
          tableOutput("end_val")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {


  iscv <- reactive({
    PowerTOST::CI2CV(lower = input$CI_lower,
                     upper = input$CI_upper,
                     n = input$N_sub,
                     design = input$design)
  })

  power_tab <- reactive({
      power.table(CV = iscv(),
                  theta0 = .95,
                  design = input$design) %>%
      as.data.frame() %>%
      `colnames<-`(c("n","dropouts","drp","power")) %>%
       mutate( n_completers = n - dropouts,
               `n with 10% dropout rate` = ceiling((n_completers / 0.9) / 2)*2,
               `Power(%)` = power * 100) %>%
      dplyr::select(`Power(%)`,n_completers, `n with 10% dropout rate`)


  })


  output$orig_val <- renderText({
    paste0( "ISCV based on the provided data: ", round( iscv() * 100, digits = 2), "%")
    })

  output$end_val  <-  renderTable({ power_tab()},
                                  align = "c",
                                  hover = TRUE,
                                  digits = c(2,0,0))

  # Server response for download
  output$downloadData <- downloadHandler(
    filename = function() {
      "calculate_daily_data.xlsx"
    },
    content = function(file) {
      file.copy("calculate_daily_data.xlsx", file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
