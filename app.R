### Generate data with 3 time points of hipp vol measurement at at 45, 52 and 59
#### just using two timepoints for now in case it helps with speed...
#### assuming that motion prob is independent between test/retest, and that prob for a short scan is 1/5 that of long scan....
### launch to shinyapps.io with rsconnect per https://shiny.rstudio.com/articles/shinyapps.html (or just push the publish button in the upper right corner when app is runnign?)

library(shiny)
library(cowplot) # for plot_grid
library(ggplot2)
library(ggpubr)
library("effsize") # inc cohen.d, wasn't working with "psych" package on shinyapps.io
library("MASS")
library(pwr)
library(lme4)
library(tidyr)
library(simstudy)
library(matrixStats) # for rowMedians

set.seed(483726)

# user define function to run simulation ----
simulate <- function(input){
  
  # hard-coded parameters (for now at least)
  N <- 1000 # nSubs
  nTimepoints <- 3
  interval <- 7 # years between time points
  nRapidScan <- 5 # number of rapid scans collected. not fully generalized when calculating avg at end of loop simulating observed scores!!
  nSim <- 100
  
  # ###### for debugging
  # input=data.frame(input_brainVar_mean=4000, input_brainVar_sd=420, input_decline_mean= .001, input_decline_sd= .0005, input_errorSD_long= 50, input_errorSD_short= 100, input_motion_prob=.1, input_motion_mag=100)
  
  # get input
  input_brainVar_mean <- input$input_brainVar_mean # true mean (e.g hippocampal volume)
  input_brainVar_sd   <- input$input_brainVar_sd   # true SD
  input_decline_mean  <- input$input_decline_mean  # true rate of annual brainvar decline
  input_decline_sd    <- input$input_decline_sd    # true SD of brainvar decline
  input_errorSD_long  <- input$input_errorSD_long  # measurement error for traditional long scan
  input_errorSD_short <- input$input_errorSD_short # measurement error for rapid scan
  input_motion_prob <- input$input_motion_prob # probability of large movement
  input_motion_mag <- input$input_motion_mag # magnitude of noise from large movement
  
  # N <- input$N
  # input_noise_sd <- input$input_noise_sd
  
  # run calculations

  ## STEP1: Define the "True" underlying longitudinal dataset for hippocampal volume that you will sample from with measurement error later on
  ## use simstudy package to define and sample simulated data
  ### for distribution of decline, experimeting wiht beta / 4 (to restrict to [0,1] / 4)  
  tdef <- defData(varname = "time0_true", dist = "normal", formula = input_brainVar_mean, variance = input_brainVar_sd^2)
  tdef <- defData(tdef, varname = "tmp_decRate", dist = "beta", formula = input_decline_mean, variance = input_decline_sd^2)
  tdef <- defData(tdef, varname = "decRate", dist = "nonrandom", formula = "tmp_decRate/4")
  tdef <- defData(tdef, varname = "time1_true", dist = "nonrandom", formula = paste0("time0_true*(1-decRate)^", interval))
  # tdef <- defData(tdef, varname = "time2_true", dist = "nonrandom", formula = paste0("time0_true*(1-decRate)^", 2*interval))
  
  ### create table with a column for each time point, i.e. wide
  dtTrial <- genData(N, tdef)
  
  hist(dtTrial$decRate)
  
  ## STEP 2: sample from "true" longitudinal dataset, with measurement error
  ### modeling measurements collected with standard T1 scan (long scan), and rapid T1 (short scan)
  ### ** for now, just picking numbers for error sds, assuming rapid is higher. can explore parameter space and/or find data to refine estimates
  ### ** for now, assuming that errors are independent!! especially relevant for repeated short scans
  
  ## simulate gaussian error and add to true scores
  ### (not sure why i have to use "get" here...)
  ### (might not be necessary to store errors but do that for now)
  for(i in c(0,1)){ # loop over time points
    # long scan
    dtTrial[, paste0("time", i, "_err_long")] <- rnorm(N, 0, input_errorSD_long) + rnorm(N,0,input_motion_mag) * (runif(N,0,1) < input_motion_prob) # draw and store error, including motion via uniform distribution test
    dtTrial[, paste0("time", i, "_obs_long")] <- dtTrial[, get(paste0("time", i, "_err_long"))] + dtTrial[, get(paste0("time", i,"_true"))] # add error to true to get Observed ("obs") 
    ## also simulate a "retest" for the routine long scan so we can calculate reliability
    dtTrial[, paste0("time", i, "_err_long.retest")] <- rnorm(N, 0, input_errorSD_long) + rnorm(N,0,input_motion_mag) * (runif(N,0,1) < input_motion_prob)
    dtTrial[, paste0("time", i, "_obs_long.retest")] <- dtTrial[, get(paste0("time", i, "_err_long.retest"))] + dtTrial[, get(paste0("time", i,"_true"))] 
    # short/rapid scans
    for(j in seq(1,nRapidScan)){ 
      dtTrial[, paste0("time", i, "_err_short", j)] <- rnorm(N, 0, input_errorSD_short) + rnorm(N,0,input_motion_mag) * (runif(N,0,1) < input_motion_prob/5) 
      dtTrial[, paste0("time", i, "_obs_short", j)] <- dtTrial[, get(paste0("time", i, "_err_short", j))] + dtTrial[, get(paste0("time", i,"_true"))] 
      ## short "retest" so we can get reliability of avg short scans
      dtTrial[, paste0("time", i, "_err_short.retest", j)] <- rnorm(N, 0, input_errorSD_short) + rnorm(N,0,input_motion_mag) * (runif(N,0,1) < input_motion_prob/5) 
      dtTrial[, paste0("time", i, "_obs_short.retest", j)] <- dtTrial[, get(paste0("time", i, "_err_short.retest", j))] + dtTrial[, get(paste0("time", i,"_true"))] 
    }
    ## now avg short scans
    ### (there has to be some way to do this dynamically but couldn't get first try to work! dtTrial[, get(paste0("O",i,"_S",seq(1,nRapidScan)))])
    dtTrial[, paste0("time", i, "_obs_short.avg")] <- rowMeans(cbind(dtTrial[, get(paste0("time", i, "_obs_short1"))], 
                                                                     dtTrial[, get(paste0("time", i, "_obs_short2"))], 
                                                                     dtTrial[, get(paste0("time", i, "_obs_short3"))], 
                                                                     dtTrial[, get(paste0("time", i, "_obs_short4"))], 
                                                                     dtTrial[, get(paste0("time", i, "_obs_short5"))]))
    dtTrial[, paste0("time", i, "_obs_short.med")] <- rowMedians(cbind(dtTrial[, get(paste0("time", i, "_obs_short1"))], 
                                                                     dtTrial[, get(paste0("time", i, "_obs_short2"))], 
                                                                     dtTrial[, get(paste0("time", i, "_obs_short3"))], 
                                                                     dtTrial[, get(paste0("time", i, "_obs_short4"))], 
                                                                     dtTrial[, get(paste0("time", i, "_obs_short5"))]))    
    dtTrial[, paste0("time", i, "_obs_short.avg.retest")] <- rowMeans(cbind(dtTrial[, get(paste0("time", i, "_obs_short.retest1"))], 
                                                                            dtTrial[, get(paste0("time", i, "_obs_short.retest2"))], 
                                                                            dtTrial[, get(paste0("time", i, "_obs_short.retest3"))], 
                                                                            dtTrial[, get(paste0("time", i, "_obs_short.retest4"))], 
                                                                            dtTrial[, get(paste0("time", i, "_obs_short.retest5"))]))
    dtTrial[, paste0("time", i, "_obs_short.med.retest")] <- rowMedians(cbind(dtTrial[, get(paste0("time", i, "_obs_short.retest1"))], 
                                                                            dtTrial[, get(paste0("time", i, "_obs_short.retest2"))], 
                                                                            dtTrial[, get(paste0("time", i, "_obs_short.retest3"))], 
                                                                            dtTrial[, get(paste0("time", i, "_obs_short.retest4"))], 
                                                                            dtTrial[, get(paste0("time", i, "_obs_short.retest5"))]))    
  }
  ### difference (/decline) - *** just doing 0 to 1 for now!! ***
  dtTrial$time01_true_dec          <- dtTrial$time1_true          - dtTrial$time0_true
  dtTrial$time01_obs_dec_long      <- dtTrial$time1_obs_long      - dtTrial$time0_obs_long
  dtTrial$time01_obs_dec_short.avg <- dtTrial$time1_obs_short.avg - dtTrial$time0_obs_short.avg
  dtTrial$time01_obs_dec_short.med <- dtTrial$time1_obs_short.med - dtTrial$time0_obs_short.med
  dtTrial$time01_obs_dec_short1    <- dtTrial$time1_obs_short1    - dtTrial$time0_obs_short1
  
  ## look at reliability, lazy right now with cor but ***could update to ICC; could add time1 and time2 but shoudl be same***
  rel_long <- cor(dtTrial[, c("time0_obs_long", "time0_obs_long.retest")])[1,2]
  rel_short.avg <- cor(dtTrial[, c("time0_obs_short.avg", "time0_obs_short.avg.retest")])[1,2]
  rel_short.med <- cor(dtTrial[, c("time0_obs_short.med", "time0_obs_short.med.retest")])[1,2]
  x <- cor(dtTrial[, c("time0_obs_short1", "time0_obs_short2", "time0_obs_short3", "time0_obs_short4", "time0_obs_short5" )])
  rel.avg_short <- mean(x[lower.tri(x)]) # average of single short scan reliabilities
  ## cor with true score
  ### (maybe won't even use the raw measure acc scores...)
  acc_long <- cor(dtTrial[, c("time0_obs_long", "time0_true")])[1,2]
  acc_short1 <- cor(dtTrial[, c("time0_obs_short1", "time0_true")])[1,2]
  acc_short.avg <- cor(dtTrial[, c("time0_obs_short.avg", "time0_true")])[1,2]
  acc_short.med <- cor(dtTrial[, c("time0_obs_short.med", "time0_true")])[1,2]
  acc_dec_long <- cor(dtTrial[, c("time01_true_dec", "time01_obs_dec_long")])[1,2]
  acc_dec_short1 <- cor(dtTrial[, c("time01_true_dec", "time01_obs_dec_short1")])[1,2]
  acc_dec_short.avg <- cor(dtTrial[, c("time01_true_dec", "time01_obs_dec_short.avg")])[1,2]
  acc_dec_short.med <- cor(dtTrial[, c("time01_true_dec", "time01_obs_dec_short.med")])[1,2]
  mean_dec_true      <- mean(dtTrial$time01_true_dec)
  mean_dec_long      <- mean(dtTrial$time01_obs_dec_long)
  mean_dec_short.avg <- mean(dtTrial$time01_obs_dec_short.avg)
  mean_dec_short.med <- mean(dtTrial$time01_obs_dec_short.med)
  mean_dec_short1    <- mean(dtTrial$time01_obs_dec_short1)
  pval_true <- t.test(dtTrial$time01_true_dec)$p.value
  pval_long <- t.test(dtTrial$time01_obs_dec_long)$p.value
  pval_short.avg <- t.test(dtTrial$time01_obs_dec_short.avg)$p.value
  pval_short.med <- t.test(dtTrial$time01_obs_dec_short.med)$p.value
  pval_short1 <- t.test(dtTrial$time01_obs_dec_short1)$p.value
  
  ### (porting these over from Rmd script a little lazily; could clean up / make more concise)
  ### (this data frame is a little awkward with replicating the true score, but doing it for ggplot)
  sim.df <- data.frame( rate = c(dtTrial$decRate, dtTrial$decRate, dtTrial$decRate),
                        true = c(dtTrial$time01_true_dec, dtTrial$time01_true_dec, dtTrial$time01_true_dec),
                        observed = c(dtTrial$time01_obs_dec_long, dtTrial$time01_obs_dec_short.avg, dtTrial$time01_obs_dec_short1),
                        obs_raw0 = c(dtTrial$time0_obs_long, dtTrial$time0_obs_short.avg, dtTrial$time0_obs_short1),
                        obs_raw0.retest = c(dtTrial$time0_obs_long.retest, dtTrial$time0_obs_short.avg.retest, dtTrial$time0_obs_short.retest1),
                        type = c(rep("long",N), rep("short.avg",N), rep("short1",N))
  )
  
  stats.df <- data.frame( rel_long=rel_long, rel_short.avg=rel_short.avg, rel.avg_short=rel.avg_short, rel_short.med=rel_short.med,
                          acc_long=acc_long, acc_short1=acc_short1, acc_short.avg=acc_short.avg,acc_short.med=acc_short.med,
                          acc_dec_long=acc_dec_long, acc_dec_short1=acc_dec_short1, acc_dec_short.avg=acc_dec_short.avg, acc_dec_short.med=acc_dec_short.med,
                          mean_dec_true=mean_dec_true, mean_dec_long=mean_dec_long, mean_dec_short1=mean_dec_short1, mean_dec_short.avg=mean_dec_short.avg,  mean_dec_short.med=mean_dec_short.med,
                          pval_true=pval_true, pval_long=pval_long, pval_short.avg=pval_short.avg, pval_short1=pval_short1, pval_short.med=pval_short.med
                          )
  
  return(list(sim=sim.df, stats=stats.df))
  
}

# Define UI for app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel(h1("Simulation: Detecting individual differences in decline with different measurements", align = "center"), windowTitle="Measuring Decline Simulation"),

  fluidRow(
    column(10, offset=1,
           p(style = "text-align:center;", "With this app, you can ..."),
    )
  ),

  hr(),

  # Sidebar layout with input and output definitions ----
  fluidRow(

    # Sidebar panel for inputs ----
    column(2, style = "background-color:#f4f4f4;",

           # Input sliders
           h4("Parameters"),
           sliderInput(inputId = "input_brainVar_mean", label = "Group Baseline Mean", min = 2000, max = 6000, value = 4322),
           sliderInput(inputId = "input_brainVar_sd",   label = "Group Baseline Variation", min = 100, max = 1000, value = 423),
           sliderInput(inputId = "input_decline_mean",  label = "Group Decline Rate Parameter 1", min = .00001, max = .5, value = .05),
           sliderInput(inputId = "input_decline_sd",    label = "Group Decline Rate Parameter 2", min = .00001, max = 10, value = 9),
           sliderInput(inputId = "input_errorSD_long",  label = "Measurement Error (long scan)", min = 1, max = 1000, value = 50),
           sliderInput(inputId = "input_errorSD_short", label = "Measurement Error (short scan)", min = 1, max = 1000, value = 100),
           sliderInput(inputId = "input_motion_prob", label = "Motion probability (over long length = 5x short length)", min = 0, max = 1, value = .1),
           sliderInput(inputId = "input_motion_mag", label = "Motion noise/error", min = 1, max = 1000, value = 100),
           actionButton("run", label = "Re-run simulation"),
           htmlOutput("text_params"),
           br()
           
    ),

    # Main panel for displaying outputs ----
    column(8,

           h4("Results - distribution of decline rate and 'true' decline values"),
           br(),
           plotOutput(outputId = "hist"),
           br(),
           h4("Results - single simulation of two time points"),
           br(),
           plotOutput(outputId = "distPlot"),
           br(),
           h4("Results - distributions over 100 simulations of two time points"),
           br()
           # plotOutput(outputId = "distPlot100")

    ),

    # Sidebar panel (R) ----
    column(2, style = "background-color:#f4f4f4;",

           h4("Example scenarios"),
           p(strong("1: Modeling hippocampal volume."),  "Parameters based on empirical data...."),
           actionButton("hippoVolume", label = "Run this scenario"),
           br(),
           br()
           # p(strong("2: Reliable, Well-powered."),  "With a reliable measure and lots of subjects, you can be well-powered to detect associations with another variable even if you don't observe group-level activation."),
           # actionButton("RelPow", label = "Run this scenario"),
           # br(),
           # br()

    )

  )
)


# Define server logic required to create outputs ----
server <- function(input, output, session) {
  
  vals <- reactiveValues()
  observe({

    input$run # access re-run button, so that if it's clicked this will run again

    # run 1 iteration
    result <- simulate(input)

    # also run 100 iterations
  # results <- data.frame(iter=integer(), errSD_long=double(), errSD_short=double(),
  # type=character(), reliability=double(), accuracy=double(), dec_accuracy=double(), pval=double())
  
    stats100.df <- data.frame( )
    # for( i in 1:100 ){
    #   stats.df <- simulate(input)$stats
    #   stats100.df <- rbind(stats100.df, stats.df)
    # }

    # save variables
    vals$sim.df <- result$sim
    vals$stats.df <- result$stats
    vals$stats100.df <- stats100.df

  })

  # Update for different example scenarios ----
  ### Mean and variance 4000, 420 per https://doi.org/10.1016/j.nicl.2019.101904
  ### Mean and variance 4322, 423 per Dunedin p45
  observeEvent(input$hippoVolume, {
    updateSliderInput(session, "input_brainVar_mean", value = 4322)
    updateSliderInput(session, "input_brainVar_sd", value = 423)
    updateSliderInput(session, "input_decline_mean", value = .05)
    updateSliderInput(session, "input_decline_sd", value = 9)
    updateSliderInput(session, "input_errorSD_long", value = 50)
    updateSliderInput(session, "input_errorSD_short", value = 100)  
    updateSliderInput(session, "input_motion_prob", value = .1)
    updateSliderInput(session, "input_motion_mag", value = 100)      
  })

  # # observeEvent(input$RelPow, {
  # #   updateSliderInput(session, "input_mean", value = 0)
  # #   updateSliderInput(session, "input_sd", value = 5)
  # #   updateSliderInput(session, "N", value = 2000)
  # #   updateSliderInput(session, "input_noise_sd", value = .8)
  # # })
  # 
  # Main center panel ----
  output$hist <- renderPlot({
    
    N <- input$N
    df <- vals$sim.df
    
    p_rate <- ggplot(df, aes(rate)) + geom_histogram() + ggtitle(paste0("Mean decline rate: ", round(mean(df$rate), 3), " (sd=", round(sd(df$rate), 3), ")"))
    p_true <- ggplot(df, aes(true)) + geom_histogram() + ggtitle(paste0("Mean decline: ", round(mean(df$true), 3), " (sd=", round(sd(df$true), 3), ")"))
    
    plot_grid(p_rate, p_true, ncol = 2, align = "hv")
    
  })
  
  output$distPlot <- renderPlot({

    N <- input$N
    df <- vals$sim.df

    stats.df <- vals$stats.df

    p_dec <- ggplot(df, aes(x=true, y=observed, color=type)) + facet_grid(cols=vars(type)) + geom_point(size=3, alpha=.2) + geom_abline(slope=1, intercept = 0) +
      theme(legend.position = "none") + xlab("True decline") + ylab("Observed decline") +
      ggtitle(paste0("TRUE DECLINE VS OBSERVED DECLINE FROM BASELINE TO TIME2\n",
                    "True mean decline: ", round(stats.df$mean_dec_true, 3), " (p=", round(stats.df$pval_true, 3), ")\n",
                     "Obs mean decline (long): ", round(stats.df$mean_dec_long, 3), " (p=", round(stats.df$pval_long, 3), ") (acc=", round(stats.df$acc_dec_long, 3), ")\n",
                     "Obs mean decline (short.avg): ", round(stats.df$mean_dec_short.avg, 3), " (p=", round(stats.df$pval_short.avg, 3), ") (acc=", round(stats.df$acc_dec_short.avg, 3), ")\n",
                     "Obs mean decline (short.med): ", round(stats.df$mean_dec_short.med, 3), " (p=", round(stats.df$pval_short.med, 3), ") (acc=", round(stats.df$acc_dec_short.med, 3), ")\n",
                     "Obs mean decline (short1): ", round(stats.df$mean_dec_short1, 3), " (p=", round(stats.df$pval_short1, 3), ") (acc=", round(stats.df$acc_dec_short1, 3), ")\n"))
    p_rel <- ggplot(df, aes(x=obs_raw0, y=obs_raw0.retest, color=type)) + facet_grid(cols=vars(type)) + geom_point(size=3, alpha=.5) + geom_abline(slope=1, intercept = 0) +
      xlab("Observed raw measure (baseline)") + ylab("Observed raw measure (baseline retest)") +
      ggtitle(paste0("BASELINE TEST-RETEST RELIABILITY OF MEASURE\n",
                     "Obs reliability (long): ", round(stats.df$rel_long, 3), "\n",
                     "Obs reliability (short.avg): ", round(stats.df$rel_short.avg, 3), "\n",
                     "Obs reliability (short.med): ", round(stats.df$rel_short.med, 3), "\n",
                     "Obs reliability (short1): ", round(stats.df$rel.avg_short, 3)))
    plot_grid(p_dec, p_rel, ncol = 2, align = "hv")

  })
  
  
  # output$distPlot100 <- renderPlot({
  #   
  #   df <- vals$stats100.df
  #   Nsim <- nrow(df)
  #   
  #   df_long <- data.frame(value=c(df$true_pval, df$t1_pval, df$t2_pval, df$d_true, df$d_t1, df$d_t2, df$c1, df$c2),
  #                         type=c(rep("pval", 3*Nsim), rep("d", 3*Nsim), rep("c", 2*Nsim)),
  #                         instance=c(rep(c(rep("Actual/True",Nsim), rep("Obs (Time 1)",Nsim), rep("Obs (Time 2)",Nsim)), 2), rep("Obs (Time 1)",Nsim), rep("Obs (Time 2)",Nsim))
  #                         )
  #    
  #   sigP_true <- sum(df$true_pval<.05)
  #   sigP_t1 <- sum(df$t1_pval<.05)
  #   sigP_t2 <- sum(df$t2_pval<.05)
  #   sigP_t1true <- sum(df$true_pval<.05 & df$t1_pval<.05)
  #   sigP_t2true <- sum(df$true_pval<.05 & df$t2_pval<.05)
  # 
  #   p_pval <- ggplot(df_long[which(df_long$type=="pval"), ], aes(x=instance, y=value)) + geom_violin() + ylab("p-value") + ylim(c(0,1)) +xlab("") +
  #     theme(axis.text.x = element_text(face="bold", color="black", size=12)) +
  #     ggtitle(paste0("GROUP DIFFERENCES DETECTED\n",
  #                    "True # p<.05: ", sigP_true, "\n",
  #                    "Obs # p<.05 (Time1): ", sigP_t1, " (", sigP_t1true, " true)\n",
  #                    "Obs # p<.05 (Time2): ", sigP_t2, " (", sigP_t2true, " true)\n"))
  #   p_d <- ggplot(df_long[which(df_long$type=="d"), ], aes(x=instance, y=value)) + geom_violin() + ylab("Cohen's d") + ylim(c(0,5)) +xlab("") +
  #     theme(axis.text.x = element_text(face="bold", color="black", size=12)) +
  #     ggtitle(paste0("EFFECT SIZES\n",
  #                    "Mean d (true): ", round(mean(df$d_true),3), "\n",
  #                    "Mean d (obs time 1): ", round(mean(df$d_t1),3), "\n",
  #                    "Mean d (obs time 2): ", round(mean(df$d_t2),3), "\n"))
  #   p_c <- ggplot(df_long[which(df_long$type=="c"), ], aes(x=instance, y=value)) + geom_violin() + ylab("Reliability") + ylim(c(0,1)) +xlab("") +
  #     theme(axis.text.x = element_text(face="bold", color="black", size=12)) +
  #     ggtitle(paste0("TEST-RETEST RELIABILITY\n",
  #                    "Mean reliability (group 1): ", round(mean(df$c1),3), "\n",
  #                    "Mean reliability (group 2): ", round(mean(df$c2),3), "\n"))
  #     scale_x_discrete(labels=c("Group 1", "Group 2"))
  # 
  #   plot_grid(p_pval, p_d, p_c, ncol = 3, align = "hv")
  #   
  # })
  # 
  # 
  output$text_params <- renderText({
    paste("<br>",
          "<font size=2><p><b>Expected Group Difference:</b> Need to update!!.</p>",
          "<p><b>True Population Variation:</b> ....</p></font>",
          sep="")
  })
  # 
  # output$text <- renderText({
  #   
  #   df <- vals$stats.df
  # 
  #   paste("<h4>Results</h4>",
  #         "<p>Gp dif (Cohen's d) at time 1: <b>", round(df$d_t1, 3), "</b>.  p=", round(df$t1_pval, 3), "</b><br>",
  #         "Gp dif (Cohen's d) at time 2: <b>", round(df$d_t2, 3), "</b>.  p=", round(df$t2_pval, 3), "</b></p>",
  #         "(The true effect size is: <b>", round(df$d_true, 3), ")</p>",
  #         "<p>Reliability for Gp 1: <b>", round(df$c1, 3), "</b></p>",
  #         "<p>Reliability for Gp 2: <b>", round(df$c2, 3), "</b></p>",
  #         sep="")
  # 
  # })
  
}


shinyApp(ui = ui, server = server)


