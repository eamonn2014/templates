#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
# visualising interactions
# need to fix the fun=exp on contrasts. update march21 interaction plot.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    rm(list=ls())   
    set.seed(3333) # reproducible
    library(shiny) 
    library(shinyWidgets)
    library(shinythemes)  # more funky looking apps
    library(shinyalert)
    library(reshape)
    library(rms)

# design matrix
# interprtation next to forest plots 
    options(max.print=1000000)    
    options(scipen=999)
    fig.width <- 400
    fig.height <- 300
    fig.width1 <- 1380
    fig.height1 <- 700
    fig.width2 <- 1400
    fig.height2 <- 300
    fig.width3 <- 1400  
    fig.height3 <- 600
    fig.width4 <- 1380
    fig.height4 <- 450
    fig.width5 <- 1380
    fig.height5 <- 225
    fig.width6 <- 400
    fig.height6 <- 550
    fig.width7 <- 600
    fig.widthx <- 593
    fig.heightx <- 268
    fig.height7 <- 600
    fig.width9 <- 1380
    fig.height9 <- 679

## convenience functions
    p0 <- function(x) {formatC(x, format="f", digits=0)}
    p1 <- function(x) {formatC(x, format="f", digits=1)}
    p2 <- function(x) {formatC(x, format="f", digits=2)}
    p3 <- function(x) {formatC(x, format="f", digits=3)}
    p4 <- function(x) {formatC(x, format="f", digits=4)}
    p5 <- function(x) {formatC(x, format="f", digits=5)}
    logit <- function(p) log(1/(1/p-1))
    expit <- function(x) 1/(1/exp(x) + 1)
    inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
    is.even <- function(x){ x %% 2 == 0 } # function to id. odd maybe useful
    options(width=200)
    
    getmode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }



############################################################
#### now let us write a plot function to present interactions
###########################################################

i.plot <- function( v , M,  N,   M1,  N1, k1, pv, double) {
    
    # v is the factor of interest
    # M treatment level
    # N treatment level
    # M1 factor level 
    # N1 factor level
    
    # original stats are on log scale
    zz <- k1
    Scorex=as.vector(zz$Contrast)
    lbx =  as.vector(zz$Lower)
    ubx =  as.vector(zz$Upper)
    
    #M=M
    #N=N
    vx=v
    effect = names(k1)[1]
    
    # create a data set where we exponentiate the stats
    df.plot <- data.frame(factor.=c(M1,N1 ),
                          x=paste0("Treatment ",N," - Treatment ",M), 
                          Score=exp(Scorex),
                          lb = exp(lbx),
                          ub =exp(ubx)
    )
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #https://stackoverflow.com/questions/11094822/numbers-in-geometric-progression
    #https://stackoverflow.com/questions/5046026/print-number-as-reduced-fraction-in-r
    v <- 2^seq(-8, 8, by=1)
    v2 <- as.character(MASS::fractions(v)) # labels for axis, mix of fractions and numeric
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # get the y axis the same whether 0,1 or 1,0 ?
    if (M1 < N1) {
        
        df.plot$factor. = factor(df.plot$factor., levels = c(M1,N1 ))
        
    } else {
        
        df.plot$factor. = factor(df.plot$factor., levels = c(N1,M1 ))
        
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # capture interaction effect to present on graph
    # difference on log odds scale, could use abs(diff(Scorex))
    # interaction. <- max(Scorex[2],Scorex[1]) -  min(Scorex[2],Scorex[1])  
    interaction. <-  double$Contrast[[1]]
    
    gp <- ggplot(df.plot, aes(x=factor., y=log(Score), fill="black", group=x))
    gg <- gp + #geom_line(aes(linetype=x), size=.6) + 
        geom_point(aes(shape=x), size=4,  color="blue") + 
        geom_errorbar(aes(ymax=log(ub), ymin=log(lb)), width=0.1, size=1, color="blue") +
        theme(legend.position="none") + ylab("Odds Ratio (OR > 1 better outcomes) ") + xlab( vx) +
        
        theme_bw() +
        
        # tick labels and axis labels
        theme(axis.text.x=element_text(size=14),
              axis.title.x=element_text(size=14,face="bold")) +
        
        theme(axis.text.y=element_text(size=14),
              axis.title.y=element_text(size=14,face="bold")) +
        
        theme(plot.title = element_text(size = 16, face = "bold")) +
        
        theme(legend.position="none") +
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # line of no effect
        geom_hline(yintercept=log(1), linetype="dashed", color = "blue") +
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        scale_y_continuous(
            breaks= log(v)  ,  
            limits = c(log(min(v)),log(max(v))),  
            label=     v2  # created earlier
        ) +
        
        coord_flip() +
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        
        geom_text(aes(   
            y=log(40),
            label = paste0(p3(Score),", 95%CI (" ,p3(lb),", ",p3(ub), ")"), 
            vjust=-1.0), size=5.8, color='black') +
        
        ggtitle( paste0("Adjusted odds ratio for treatment effect (Treatment ",N," - Treatment ",M,") for ",vx,", interaction p-value ",p5(pv) ) )
    
    gg <- gg + labs(caption = c(paste0("The p-value tests for the necessity of the interaction printed in orange, it is the result of a hypothesis test assessing the interaction with treatment alone."))) + 
        
        theme(plot.caption = element_text( size = 14, face = "bold")) 
    
    # Add arrows
    i <- gg + geom_segment(
        x = 1.5, y =  Scorex[1],        # y start of arrow
        xend = 1.5, yend =  Scorex[2],  # end of arrow yend 
        lineend = "round",              # more available arrow types 
        linejoin = "round",
        size = .5, 
        arrow = arrow(length = unit(0.2, "cm")),
        colour = "#EC7014" # Also accepts "red", "blue' etc
    ) 
    
    # now add text , we exponentiate the dif of the log odds ratios and show the interaction form both points of view
    k <- i + geom_text( aes(
        x = 1.4, #y = (Scorex[1]+Scorex[2])/2,
        y=log(75),
        label = paste0("Interaction multiplication factor:\n ",p3(exp(double$Contrast)),", 95%CI (" ,p3(exp(double$Lower)),", ",p3(exp(double$Upper)), ")"), 
        group = NULL,
        vjust = -1, #.3
        hjust = .7 #1
    ), size=5.8 , color="#EC7014") 
    
}

# function ended
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <- fluidPage(theme = shinytheme("journal"), #https://www.rdocumentation.org/packages/shinythemes/versions/1.1.2
                # paper
                useShinyalert(),  # Set up shinyalert
                setBackgroundColor(
                    color = c( "#2171B5", "#F7FBFF"), 
                    gradient = "linear",
                    direction = "bottom"
                ),
                
                h2("Differential Treatment Effects (Double Differences)"), 
                h4("It is often desired to investigate if there is evidence for different treatment effects depending on the levels of baseline factor variables or the level 
                of continuous variables following an RCT. One or more interactions between baseline covariates and treatment are then explored.
                Here we simulate an RCT with a binary response, 3 treatment arms and 11 baseline covariates. Note in reality this objective will be extremely underpowered, typically one wants to detect a
                differential effect that is smaller than the overall detectable treatment effect [1].  'It is important to note that assessing treatment effect 
                in an isolated subgroup defined by a categorical covariate does not establish differential treatment effects and results in unreliable estimates. 
                Differential treatment effect must be demonstrated.' [2]. However, before embarking on this endeavour it is advisable to assess the variation in the response in the treatment arms first [3].
         "), 
                
                h3("  "), 
                
                sidebarLayout(
                    
                    sidebarPanel( width=3 ,
                                  
                                  tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                  
                                  
                                  actionButton(inputId='ab1', label="R Shiny ",   icon = icon("th"),   
                                               onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/differential-treatment-effects/master/dif.trt.effects/app.R', '_blank')"), 
                                  actionButton(inputId='ab1', label="R code",   icon = icon("th"),   
                                               onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/differential-treatment-effects/master/R%20raw%20code.R', '_blank')"),  
                                  actionButton("resample", "Simulate a new sample"),
                                  br(),  
                                  tags$style(".well {background-color:#b6aebd ;}"), 
                                  
                                  h4(" 
                                   "),
                                  div(
                                      
                                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                      
                                      tags$head(
                                          tags$style(HTML('#ab1{background-color:orange}'))
                                      ),
                                      
                                      tags$head(
                                          tags$style(HTML('#resample{background-color:orange}'))
                                      ),
                                      #tags$hr(),
                                      textInput('n',
                                                div(h5(tags$span(style="color:blue", "Sample size"))), value= "1000"),
                                      
                                      
                                      
                                      selectInput("Design",
                                                  div(h5(tags$span(style="color:blue", "Select design preference:"))),
                                                  
                                                  choices=c(  "No-interaction logit-additive model", 
                                                              "Treatment interacts with smoking only" ,
                                                              "Treatment interacts with all variables" 
                                                  ), width='80%'),
                                      
                                      
                                      
                                      selectInput("Model",
                                                  div(h5(tags$span(style="color:blue", "Select modelling preference (impacts Table 1 & tab 10 & 11):"))),
                                                  choices=c(  "No-interaction logit-additive model",
                                                              "Treatment interacts with smoking only" ,
                                                              "Treatment interacts with all variables"
                                                  ), width='80%'),
                                      
                                      div(h5(tags$span(style="color:blue", "Enter treatment coefficients:"))),
                                      
                                      splitLayout(
                                          textInput("v1", div(h5(tags$span(style="color:blue", "Treatment 3 arms"))), value= "1"),
                                          textInput("v2", div(h5(tags$span(style="color:blue", "Age (continuous)"))), value= "1/(65-18)"),
                                          textInput("v3", div(h5(tags$span(style="color:blue", "Smoking (factor)"))), value= "0.4")
                                          
                                      ),
                                      
                                      splitLayout(
                                          textInput("v4", div(h5(tags$span(style="color:blue", "BMI (factor)"))), value= "0"),
                                          textInput("v5", div(h5(tags$span(style="color:blue", "covar3 (biomarker)"))), value= "1/3"),
                                          textInput("v6", div(h5(tags$span(style="color:blue", "covar1 (Blood score)"))), value= "-.5/10")
                                          
                                      ),
                                      
                                      
                                      splitLayout(
                                          textInput("v7", div(h5(tags$span(style="color:blue", "Vas (continuous)"))), value= "0.25/30"),
                                          textInput("v8", div(h5(tags$span(style="color:blue", "Time (continuous)"))), value= "-.1/10"),
                                          textInput("v9", div(h5(tags$span(style="color:blue", "covar2 (Fitness score)"))), value= "-1/50")
                                          
                                      ),
                                      
                                      
                                      splitLayout(
                                          textInput("v10", div(h5(tags$span(style="color:blue", "fact1 (History binary)"))), value= "log(2)"),
                                          textInput("v11", div(h5(tags$span(style="color:blue", "binary2 (Employed)"))), value= "-log(1)"),
                                          textInput("v12", div(h5(tags$span(style="color:blue", "Sex (binary)"))), value= "log(0.5)")
                                          
                                      ),
                                      
                                      div(h4("References:")),
                                      tags$a(href = "https://statmodeling.stat.columbia.edu/2018/03/15/need-16-times-sample-size-estimate-interaction-estimate-main-effect/", tags$span(style="color:blue", "[1] Andrew Gelman"),),
                                      div(p(" ")),
                                      tags$a(href = "https://www.fharrell.com/post/varyor/", tags$span(style="color:blue", "[2] Frank Harrell...much more here"),),
                                      div(p(" ")),
                                      tags$a(href = "https://eamonn.shinyapps.io/responder-non-responder-fallacy-in-RCTs/", tags$span(style="color:blue", "[3] Responder non responder fallacy"),),
                                      div(p(" ")),
                                      tags$a(href = "https://projecteuclid.org/download/pdfview_1/euclid.aoas/1231424214",  tags$span(style="color:blue", "[4] Andrew Gelman"),),
                                      div(p(" ")),
                                      
                                  )
                                  
                    ),
                    
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tab panels
                    mainPanel(width=9,
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              navbarPage(       
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                                  tags$style(HTML("
                            .navbar-default .navbar-brand {color: orange;}
                            .navbar-default .navbar-brand:hover {color: blue;}
                            .navbar { background-color: #b6aebd;}
                            .navbar-default .navbar-nav > li > a {color:black;}
                            .navbar-default .navbar-nav > .active > a,
                            .navbar-default .navbar-nav > .active > a:focus,
                            .navbar-default .navbar-nav > .active > a:hover {color: pink;background-color: purple;}
                            .navbar-default .navbar-nav > li > a:hover {color: black;background-color:yellow;text-decoration:underline;}
                            .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
                            .navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
                            .navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
                   ")),
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  
                                  tabPanel("1 Design and modelling selection",  
                                           fluidRow(
                                               column(width = 5, offset = 0, style='padding:1px;',
                                                      
                                                      h4(paste("An explanation of the inputs and tabs")), 
                                                      h4(paste("The first input left is the sample size for patients randomly 
                                              assigned to treatment arms in a 1:1:1 fashion.
                                                       The second selection allows the choice of 3 designs i) a main effects model, that is a  
                                                       no-interaction logit-additive model that assumes constancy of treatment ORs ii) 
                                                       a model with a treatment X smoking interaction and iii) a model in which all baseline covariates interact 
                                                       with treatment. ")),
                                                      #br(),
                                                      h4(paste("The next selection is a choice of analysis performed and presented in Table 1. 
                                                       There are three choices once again i) a main effects model, that is a 
                                                       no-interaction logit-additive model that assumes constancy of treatment ORs ii) 
                                                       a model with a treatment X smoking interaction and iii) a model in which all baseline covariates interact 
                                                       with treatment. This only impacts what is presented in Table 1 and tab 10/11.")), 
                                                      #br(),
                                                      h4(paste("Twelve input boxes follow and allow the user to specify the coefficients for treatment and 11
                                              baseline covariates on the log odds scale.
                                                       Note a typical change in an input variable would be unlikely to correspond to a change as large 
                                                       as 5 on the logistic scale (which would move the probability from 0.01 to 0.50 or from 0.50 to 0.99) [4].
                                                       Age in years is uniformly distributed between 18 and 65. covar3 is uniformly distributed
                                                       between 0 to 3, covar1 uniformly distributed between 0 to 10, vas between 1 to 30 and time in years uniformly
                                                       distributed between 0 to 10. Smoking and BMI are 3 level factors and fact1, binary2 and sex are binary factors.
                                                       For the factors the coefficient entered describes the true relationship between all adjacent levels. We also 
                                                       add labels to the variables and they appear on some of the outputs."  )), 
                                                      h4(paste("Click the simulate button to generate another data set from the same population.")),
                                                      h4(paste("Tab 1 presents the regression table of the 3 models, the particular model can be selected. ")),
                                                      h4(paste("Tab 2 presents the regression table of the 
                                                        no-interaction logit-additive model that assumes constancy of treatment ORs and provides an explanation.
                                                        ")),
                                                      h4(paste("Tab 3 presents the a forest plot of the no-interaction logit-additive model that assumes constancy of treatment ORs plus 
                                                       a table of the ORs and log odds ratios.")),
                                                      h4(paste("Tab 4 presents the likelihood ratio test assessing each model with each other.")),
                                                      
                                                      h4(paste("Tab 5 presents the forest plots by treatment for the treatment interacting with all baseline covariates model plus 
                                                       tables of the ORs and log odds ratios and the regression tables. We allow the user to adjust reference levels 
                                                       and the range for which continuous predictor effects are estimated and presented. 
                                                       For good measure we also describe a couple of the estimated regression coefficients.")),
                                                      h4(paste("Tab 6 presents the forest plots by treatment for the treatment x smoking interaction model and 
                                                       tables of the ORs and log odds ratios.")),
                                                      h4(paste("Tab 7 presents relative measures of explained variation and the AIC of each model is reported.")),
                                                      h4(paste("Tab 8 tab presents anova table and dot plot.")),
                                                      h4(paste("Tab 9 tab presents double differences on the log odds scale (ie. the interactions) .")),
                                                      
                                                      h4(paste("Tab 10 No new information is presented on this tab, but only the 3 model outputs presented together.")),
                                                      
                                                      h4(paste("The 11th  tab presents another way to estimate treatment effects with interactions using contrast statements.")),
                                                      
                                                      h4(paste("The final tab presents a listing of the simulated data and diagnostics."))
                                               ),
                                               
                                               fluidRow(
                                                   column(width = 1, offset = 0, style='padding:1px;',
                                                          
                                                   ),
                                                   
                                                   fluidRow(
                                                       column(width = 5, offset = 0, style='padding:1px;',
                                                              h4(htmlOutput("textWithNumber2",) ),
                                                              div( verbatimTextOutput("user") )
                                                       ))),
                                           )
                                  ),
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  tabPanel("2 No-interaction model",  
                                           h4(paste("Table 2. No-interaction logit-additive model that assumes constancy of treatment ORs")), 
                                           
                                           fluidRow(
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      
                                                      div( verbatimTextOutput("Cx2") )
                                                      
                                                      
                                               ),
                                               h4(paste("An explanation of the inputs..Do we recover the true population values?")), 
                                               fluidRow(
                                                   
                                                   column(width = 5, offset = 0, style='padding:1px;',
                                                          
                                                          h4(htmlOutput("textWithNumber",) ),
                                                   ))),
                                  ),
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  tabPanel("3 Forest plot, no-interaction model", value=7, 
                                           h4(paste("Figure 1 & Table 3. No-interaction logit-additive model that assumes constancy of treatment ORs")), 
                                           fluidRow(
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      
                                                      div(plotOutput("f.plot3", width=fig.height1, height=fig.height9)),   
                                                      
                                               ) ,
                                               
                                               fluidRow(
                                                   column(width = 6, offset = 0, style='padding:1px;',
                                                          
                                                          div( verbatimTextOutput("int.trt1C" ) ))
                                                   
                                               ))),
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  tabPanel("4 Likelihood ratio test", value=7, 
                                           
                                           fluidRow(
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      h4("Table 4 [No-interaction logit-additive model that assumes constancy of treatment ORs] vrs [Treatment x Smoking interaction model]"), 
                                                      div( verbatimTextOutput("L1c") ),
                                                      h4("Table 5 [No-interaction logit-additive model that assumes constancy of treatment ORs] vrs [Treatment x all predictors interaction model]"), 
                                                      div( verbatimTextOutput("L1b") ),
                                                      h4("Table 6 [Treatment x Smoking interaction model] vrs [Treatment x all predictors interaction model]"), 
                                                      div( verbatimTextOutput("L1a") )
                                                      
                                               ) ,
                                               # h4("Table 2 xxxxxxxxxxxxxxxxx"),
                                               fluidRow(
                                                   column(width = 6, offset = 0, style='padding:1px;',
                                                          br(),br(),
                                                          h4("A small P-Value in the top most table provides evidence against
                                                        the simpler model fitting the data better. The simpler model 
                                                        being the no-interaction logit-additive model that assumes constancy of treatment ORs."),
                                                          br(),br(),br(),br(),br(),br(),br(),br() ,
                                                          h4("A small P-Value in the middle table provides evidence against
                                                        the simpler model fitting the data better. The simpler model 
                                                        being the no-interaction logit-additive model that assumes constancy of treatment ORs."),
                                                          br(),br(),br(),br(),br(),br(),br(),br() ,
                                                          h4("A small P-Value in the bottom table provides evidence against
                                                        the simpler model fitting the data better. The simpler model 
                                                        being the Treatment x Smoking interaction model."),
                                                          
                                                          splitLayout(
                                                              
                                                          ),
                                                   ))),
                                           
                                  ) ,
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  
                                  tabPanel("5 Forest plot, treatment x all variables", value=3, 
                                           h4(paste("Figure 2 Forest plots by treatment for the model in which treatment is interacted with all baseline covariates")),
                                           
                                           h4(paste("The boxes below can be used to adjust the range for which the effect is estimated for continuous predictors. 
                                                Default is from the 25th to 75th percentile of the variables distribution.")),
                                           
                                           splitLayout(
                                               textInput("age.range", div(h5(tags$span(style="color:blue", "Age (continuous)"))), value= "30, 54"),
                                               textInput("biomarker.range", div(h5(tags$span(style="color:blue", "covar3 (biomarker)"))), value= "0.7675, 2.2300"),  #18
                                               textInput("blood.range", div(h5(tags$span(style="color:blue", "covar1 (Blood score)"))), value= "2.5700, 7.7525"),
                                               textInput("vas.range", div(h5(tags$span(style="color:blue", "Vas (continuous)"))), value= "18, 23"),  #1
                                               textInput("time.range", div(h5(tags$span(style="color:blue", "Time (continuous)"))), value= "2.355, 7.420"),
                                               textInput("fitness.range", div(h5(tags$span(style="color:blue", "covar2 (Fitness score)"))), value= "13, 38") 
                                               
                                           ),
                                           
                                           
                                           
                                           
                                           
                                           h4(paste("The boxes below can be used to adjust the factor reference levels (affecting forest plot and presentation of treatment effects at very bottom). The continuous variables are held at sensible values (we did not center the continuous variables in the regression). 
                                       Set the continuous to zero and observe the treatment comparison confidence intervals. Only the treatment bars will change as treatment interacts with all variables. ")),
                                           
                                           splitLayout(
                                               textInput("adj.smoking", div(h5(tags$span(style="color:blue", "Smoking ref (factor)"))), value= "1"),
                                               textInput("adj.age", div(h5(tags$span(style="color:blue", "Age (continuous)"))), value= "40"),  #18
                                               textInput("adj.biomarker", div(h5(tags$span(style="color:blue", "covar3 (biomarker)"))), value=  "1.3"),
                                               textInput("adj.blood", div(h5(tags$span(style="color:blue", "covar1 (Blood score)"))), value= "5"),
                                               textInput("adj.vas", div(h5(tags$span(style="color:blue", "Vas (continuous)"))), value= "17"),  #1
                                               textInput("adj.time", div(h5(tags$span(style="color:blue", "Time (continuous)"))), value= "4")
                                               
                                           ),
                                           
                                           splitLayout(
                                               
                                               
                                               textInput("adj.fitness", div(h5(tags$span(style="color:blue", "covar2 (Fitness score)"))), value= "20"),  #1
                                               textInput("adj.history", div(h5(tags$span(style="color:blue", "fact1 ref (History binary)"))), value= "0"),
                                               textInput("adj.employed", div(h5(tags$span(style="color:blue", "binary2 ref (Employed)"))), value= "0"),
                                               textInput("adj.sex", div(h5(tags$span(style="color:blue", "Sex red (binary)"))), value= "0"),
                                               textInput("adj.BMI", div(h5(tags$span(style="color:blue", "BMI ref (factor)"))), value= "1")
                                             
                                               
                                           ),
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                           fluidRow(
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      div(plotOutput("f.plot1", width=fig.width4, height=fig.height7)),
                                                      
                                               )),
                                           h4(paste("Examples describing the regression coefficients. For continuous variables the default is to estimate the effect based on the 25th and 75th percentile of the variables distribution (we can change this using the top boxes). Better outcomes are indicated on the right of dotted line, worse outcomes to the left.")),
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                           #h4(paste("We can interpret the forest plot effect estimates thus...")),
                                           h4(htmlOutput("textWithNumber3",) ),
                                           
                                           fluidRow(
                                               column(12,
                                                      #  div( verbatimTextOutput("int.trt1" ) ),
                                                      fluidRow(
                                                          column(12,
                                                                 #  div( verbatimTextOutput("int.trt2" ) ),
                                                                 fluidRow(
                                                                     h4(paste("Tables 7, 8 and 9 odds ratios for treatment 1, 2 and 3")), 
                                                                     column(4, 
                                                                            div( verbatimTextOutput("int.trt1" ) )),
                                                                     
                                                                     column(4,
                                                                            div( verbatimTextOutput("int.trt2" ) )),
                                                                     
                                                                     column(4,
                                                                            div( verbatimTextOutput("int.trt3" ) )),
                                                                     
                                                                 )
                                                          )#,
                                                          #column(width = 6,
                                                          # "Fluid 6")
                                                      )
                                               )
                                           ),
                                           
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                           fluidRow(
                                               column(12,
                                                      #  div( verbatimTextOutput("int.trt1" ) ),
                                                      fluidRow(
                                                          column(12,
                                                                 #  div( verbatimTextOutput("int.trt2" ) ),
                                                                 fluidRow(
                                                                     h4(paste("Table 10, 11 and 12 Model treatment interacting with all baseline covariates, using trt reference levels 1,2,3")),  
                                                                     column(4, 
                                                                            div( verbatimTextOutput("Ax1" ) )),
                                                                     
                                                                     column(4,
                                                                            div( verbatimTextOutput("Ax2" ) )
                                                                     ),
                                                                     
                                                                     column(4,
                                                                            div( verbatimTextOutput("Ax3" ) )
                                                                     ),
                                                                     
                                                                 )
                                                          )#,
                                                          #column(width = 6,
                                                          # "Fluid 6")
                                                      )
                                                      
                                                      
                                               )
                                           ),
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                           
                                           
                                           # new march 21
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                           fluidRow(
                                               column(12,
                                             
                                                      
                                                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                      
                                                      sliderTextInput("Varsx",
                                                                      div(h5(tags$span(style="color:blue", "Select the predictor of interest "))), 
                                                                      choices = c("smoking", "age", "bmi", "covar3", "covar1", "vas", "time", 
                                                                                  "covar2", "fact1", "sex","binary2"), #add further months 
                                                                      selected = c("smoking"), #values which will be selected by default
                                                                      animate = FALSE, grid = FALSE, 
                                                                      hide_min_max = FALSE, from_fixed = FALSE,
                                                                      to_fixed = FALSE, from_min = NULL, from_max = NULL, to_min = NULL,
                                                                      to_max = NULL, force_edges = FALSE, width = NULL, pre = NULL,
                                                                      post = NULL, dragRange = TRUE),
                                                      
                                                      
                                                      
                                                      
                                                      splitLayout(
                                                          
                                                          textInput("treatment.level1x", div(h5(tags$span(style="color:blue", "Treatment level (enter 1,2 or 3) "))), value= "1"),
                                                          textInput("treatment.level2x", div(h5(tags$span(style="color:blue", "Treatment level (enter 1,2 or 3)"))), value= "2")
                                                         
                                                      ),
                                                      h4(paste("Select the predictor levels. For factors enter (1,2,3) for binary factors enter (0,1), for continuous variables enter a range for that variable.")),
                                                      splitLayout(
                                                          
                                                          textInput("interest.level1x", div(h5(tags$span(style="color:blue", "lower level  of variable of interest "))), value= "1"),
                                                          textInput("interest.level2x", div(h5(tags$span(style="color:blue", "upper level  of variable of interest"))), value= "2")
                                                          
                                                      ),
                                                      
                                                      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                      
                                                   
                                                      
                                                     # div( verbatimTextOutput("int.trtc" ) ), 
                                                      
                                                      fluidRow(
                                                          column(12, offset = 0, style='padding:1px;',
                                                                 div(plotOutput("plot.trtc", width=fig.width4, height=fig.height7)),
                                                                 fluidRow(
                                                                     h4(paste("Interaction present if the patterns differs between factors. The interaction effect printed in orange text can be seen on the log scale in the tables 10 or 11 or 12. Further the p-value can be seen in the anova table on tab 8 and selecting the 'Treatment interacts with all variables' option.")),   
                                                                     
                                                                     
                                                                     div( verbatimTextOutput("int.trtc" ) ), 
                                                                     div( verbatimTextOutput("int.trtd" ) ), 
                                                               
                                                                     
                                                                 )
                                                          )#,
                                                          #column(width = 6,
                                                          # "Fluid 6")
                                                      )
                                                      
                                                      
                                               )
                                           ), # end new mar 2021
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                     
                                  ),
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  tabPanel("6 Forest plot, treatment x smoking only", value=3, 
                                           h4(paste("Figure 3 Forest plots by treatment for the model in which treatment is interacted with smoking covariate only")), 
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                           fluidRow(
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      div(plotOutput("f.plot2", width=fig.width4, height=fig.height7)),
                                                      
                                               )),
                                           
                                           
                                           
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                           fluidRow(
                                               column(12,
                                                      #  div( verbatimTextOutput("int.trt1" ) ),
                                                      fluidRow(
                                                          column(12,
                                                                 #  div( verbatimTextOutput("int.trt2" ) ),
                                                                 fluidRow(
                                                                     h4(paste("Table 13, 14 and 15 odds ratios treatment 1, 2 and 3")), 
                                                                     column(4, 
                                                                            div( verbatimTextOutput("int.trt1B" ) )), #
                                                                     
                                                                     column(4,
                                                                            div( verbatimTextOutput("int.trt2B" ) )),
                                                                     
                                                                     column(4,
                                                                            div( verbatimTextOutput("int.trt3B" ) )),
                                                                     
                                                                 )
                                                          )#,
                                                          #column(width = 6,
                                                          # "Fluid 6")
                                                      )
                                               )
                                           ),
                                           
                                  ),
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  tabPanel("7 Relative measures of effect", value=3, 
                                           
                                           h4(htmlOutput("textWithNumber1",) ),
                                           fluidRow(
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      
                                               ) ,
                                               
                                               fluidRow(
                                                   column(width = 5, offset = 0, style='padding:1px;',
                                                          
                                                   ))),
                                  ),
                                  tabPanel("8 Anova", 
                                           
                                           fluidRow(
                                               
                                               
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      # h4("Notes"),
                                                      h4("Table 16 Analysis of variance "),
                                                      div( verbatimTextOutput("an") ),
                                                      h4("Use the select design and select modelling preference to alter the output. To approximate the anova 
                                                   table you can perform likelihood ratio tests using models with the component and without the component of interest."),
                                                      
                                                      
                                               ),
                                               column(width = 3, offset = 0, style='padding:1px;',
                                                      h4("Figure 4 Relative importance"),
                                                      div(plotOutput("anovaf", width=fig.height1, height=fig.height9)),
                                               ),
                                           )
                                  ),##end,
                                  
                                  
                                  
                                  
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ dd
                                  tabPanel("9 Double differences", value=3, 
                                           h4(paste("We are testing [(Smoking level 2 - Smoking level 1) in treatment level 2] - [(Smoking level 2 - Smoking level 1)  in treatment level 1] 
                                           this is the double difference, we can check it matches the interaction term in tab 5. This default will match the trt=2 * smoking=2 interaction.
                                        Basically we are pulling out the interaction terms. For the continuous variables, the range will make no difference to the p-value as the effects are modelled as truly linear.
                                                    For the binary variables enter 0 and 1 only.")),
                                           
                                           
                                           sliderTextInput("Vars",
                                                           div(h5(tags$span(style="color:blue", "Select the predictor of interest "))), 
                                                           choices = c("smoking", "age", "bmi", "covar3", "covar1", "vas", "time", 
                                                                       "covar2", "fact1", "sex","binary2"), #add further months 
                                                           selected = c("smoking"), #values which will be selected by default
                                                           animate = FALSE, grid = FALSE, 
                                                           hide_min_max = FALSE, from_fixed = FALSE,
                                                           to_fixed = FALSE, from_min = NULL, from_max = NULL, to_min = NULL,
                                                           to_max = NULL, force_edges = FALSE, width = NULL, pre = NULL,
                                                           post = NULL, dragRange = TRUE),
                                           
                                           
                                           
                                           
                                           splitLayout(
                                               
                                               textInput("treatment.level1", div(h5(tags$span(style="color:blue", "Treatment level (enter 1,2 or 3) "))), value= "1"),
                                               textInput("treatment.level2", div(h5(tags$span(style="color:blue", "Treatment level (enter 1,2 or 3)"))), value= "2")
                                               # textInput("interest.level1", div(h5(tags$span(style="color:blue", "Smoking level A "))), value= "1"),
                                               #textInput("interest.level2", div(h5(tags$span(style="color:blue", "smoking level B"))), value= "2")
                                               
                                           ),
                                           h4(paste("Select the predictor levels. For factors enter (1,2,3) for binary factors enter (0,1), for continuous variables enter a range for that variable.")),
                                           splitLayout(
                                               
                                               #textInput("treatment.level1", div(h5(tags$span(style="color:blue", "Treatment level A "))), value= "1"),
                                               #textInput("treatment.level2", div(h5(tags$span(style="color:blue", "Treatment level B"))), value= "2"),
                                               textInput("interest.level1", div(h5(tags$span(style="color:blue", "lower level  of variable of interest "))), value= "1"),
                                               textInput("interest.level2", div(h5(tags$span(style="color:blue", "upper level  of variable of interest"))), value= "2")
                                               
                                           ),
                                           

                                           
                                           
                                           br(), br(),
                                           fluidRow(
                                               column(12,
                                                      #  div( verbatimTextOutput("int.trt1" ) ),
                                                      fluidRow(
                                                          column(12,
                                                                 #  div( verbatimTextOutput("int.trt2" ) ),
                                                                 fluidRow(
                                                                     h4(paste("Table 17. Double differences. This should match the interaction term (ie double difference in log odds). Note the exponential of which is a multiplicative factor on the odds scale.
                                                                             
                                                                             ")), 
                                                                     column(6, 
                                                                            div( verbatimTextOutput("DD1" ) )),
                                                                     
                                                                     
                                                                     
                                                                 )
                                                          )#,
                                                          #column(width = 6,
                                                          # "Fluid 6")
                                                      )
                                               )
                                           ),
                                           
                               
                                           
                                           #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  ),
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  tabPanel("10 Print all models", value=3,    
                                           
                                           fluidRow(
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      h4(paste("Table 18. No-interaction logit-additive model that assumes constancy of treatment ORs")), 
                                                      div( verbatimTextOutput("Cx") ),
                                                      h4(paste("Table 19. Model treatment x smoking interaction only")), 
                                                      div( verbatimTextOutput("Bx") )
                                                      
                                               ),
                                               
                                               fluidRow(
                                                   column(width = 5, offset = 0, style='padding:1px;',
                                                          h4(paste("Tables 20. Model treatment interacting with all baseline covariates")), 
                                                          div( verbatimTextOutput("Ax") )
                                                   ))),
                                           
                                  ),
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                  tabPanel("11 Contrast example", value=3, 
                                           
                                           
                                           # selectInput("predictors",
                                           #               strong("predictors"),
                                           #               choices=varz),
                                           h4(paste("Table 21. Contrast output (showing design matrix) to see the effect of smoking in each treatment on the log odds scale, compare 
                                              to tab 5")), 
                                           
                                           textInput('levz', 
                                                     strong(div(h5(tags$span(style="color:blue", "Contrast to compare smoking levels")))), "2,1"),
                                           
                                           #div( verbatimTextOutput(print("z99.")   ) ),
                                           div( verbatimTextOutput("z99.") ),
                                           
                                           fluidRow(
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      #  div( verbatimTextOutput(print("z1.")   ) ),
                                               ) ,
                                               
                                               fluidRow(
                                                   column(width = 5, offset = 0, style='padding:1px;',
                                                          #     div( verbatimTextOutput(print("z2.") )),
                                                   ))),
                                  ),
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                             
                                  tabPanel("12 Data & diagnostics", 
                                           
                                           fluidRow(
                                               
                                               
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      # h4("Notes"),
                                                      
                                                      h4("Figure 5 Linear relationship between continuous predictor variables and the logit of outcome. Focus on the continuous predictors."),
                                                      div(plotOutput("residz", width=fig.height1, height=fig.height7)),
                                                      h4("Table 23 Global test of model fit, if 'Expected value|H0' is coincidental with the 'Sum of squared errors'... don't discard model"),
                                                      div( verbatimTextOutput("gofx") ),
                                                      
                                                      
                                               ),
                                               column(width = 6, offset = 0, style='padding:1px;',
                                                      h4("Table 22 Simulated data listing (note the data are simulated and result in very 'clean' data, 
                                                there are no data management errors and the data behave according to the simualted distributions,
                                                of course we can build in missing data etc. if we so wish.)"),
                                                      div( verbatimTextOutput("datx") ),
                                                      #br(),br(),
                                                      
                                                      
                                                      
                                               ),
                                           )
                                  )##end,
                                  
                                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   END NEW   
                              )
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    )
                ) 
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end tab panels 
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

server <- shinyServer(function(input, output   ) {
    
    shinyalert("Welcome! \nModelling Differential Treatment Effects!",
               "Treatment covariate interactions and double differences", 
               type = "info")
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # This is where a new sample is instigated 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    random.sample <- reactive({
        
        foo <- input$resample
        
        n <- as.numeric(input$n )
        
        # writing like this I can write log and fraction into input boxes!
        v1 <- as.numeric(    eval(parse(text= (input$v1)) ) )
        v2 <- as.numeric(    eval(parse(text= (input$v2)) ) )
        v3 <- as.numeric(    eval(parse(text= (input$v3)) ) )    
        v4 <- as.numeric(    eval(parse(text= (input$v4)) ) )   
        v5 <- as.numeric(    eval(parse(text= (input$v5)) ) )  
        v6 <- as.numeric(    eval(parse(text= (input$v6)) ) ) 
        v7 <- as.numeric(    eval(parse(text= (input$v7)) ) )
        v8 <- as.numeric(    eval(parse(text= (input$v8)) ) )
        v9 <- as.numeric(    eval(parse(text= (input$v9)) ) )
        v10 <- as.numeric(    eval(parse(text= (input$v10)) ) )
        v11 <- as.numeric(    eval(parse(text= (input$v11)) ) )
        v12 <- as.numeric(    eval(parse(text= (input$v12)) ) )
        
        check =c(v1 , v2 , v3 , v4 , v5,  v6, v7,  v8 , v9 , v10 , v11 , v12  )
        
        return(list(
            v1=v1, v2=v2, v3=v3, v4=v4, v5=v5, v6=v6, v7=v7, v8=v8, v9=v9, v10=v10, v11=v11, v12=v12 ,
            check=check, n=n
            
        ))
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    randomness <- reactive({
        
        n <- as.numeric(input$n )
        randomi <- runif(n)
        
        return(list(
            randomi=randomi
        ))
        
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    design <- reactive({
        
        randomi <- randomness()$randomi
        
        sample <- random.sample()
        
        n <- as.numeric(input$n )
        v1 <- as.numeric(    eval(parse(text= (input$v1)) ) )
        v2 <- as.numeric(    eval(parse(text= (input$v2)) ) )
        v3 <- as.numeric(    eval(parse(text= (input$v3)) ) )    
        v4 <- as.numeric(    eval(parse(text= (input$v4)) ) )   
        v5 <- as.numeric(    eval(parse(text= (input$v5)) ) )  
        v6 <- as.numeric(    eval(parse(text= (input$v6)) ) ) 
        v7 <- as.numeric(    eval(parse(text= (input$v7)) ) )
        v8 <- as.numeric(    eval(parse(text= (input$v8)) ) )
        v9 <- as.numeric(    eval(parse(text= (input$v9)) ) )
        v10 <- as.numeric(    eval(parse(text= (input$v10)) ) )
        v11 <- as.numeric(    eval(parse(text= (input$v11)) ) )
        v12 <- as.numeric(    eval(parse(text= (input$v12)) ) )
        
        trt.coef       <-  v1     # log odds ratio so 1 -> 2.718, so 1 is LARGE
        age.coef       <-  v2     # log odds of 1 over the age range
        smoke.coef     <-  v3     # this is odds of 1.5
        bmi.coef       <-  v4     # this is an odds of 1..50:50
        covar3.coef    <-  v5     # log odds 1 over range of 3
        covar1.coef    <-  v6     # log odds -.05 per unit change
        vas.coef       <-  v7     # log odds .008 per unit change. log odds .25 over 30 units odds 1.27
        time.coef      <-  v8     # log odds -.01 per year, log odds -.1 over 10 years or odds .90
        covar2.coef    <-  v9     # log odds 0.02 per joint, log odds 1 over 50 units or odds 2.7
        fact1.coef     <-  v10    # log odds 0.693 per change in binary, or odds of 2   
        binary2.coef   <-  v11    # log odds 0 per change in binary, or odds of 1  
        sex.coef       <-  v12    # log odds -0.693 per change in binary, or odds of .5  
        
        intercept <- -5
        
        trt      <- sample(1:3,   n, replace=TRUE)      # trt 3 levels
        age      <- sample(18:65, n, replace=TRUE)      # continuous
        bmi      <- sample(1:3,   n, replace=TRUE)      # assume 3 equal groups?
        smoking  <- sample(1:3,   n, replace=TRUE)      # categorical assume 3 equal groups?
        covar3   <- round(runif(n,0,3),2)
        covar1   <- round(runif(n,0,10),2)
        vas      <- sample(1:30, n, replace=TRUE)
        time     <- round(runif(n,0,10),2)              # years
        covar2   <- sample(1:50, n, replace=TRUE)
        fact1    <- sample(0:1,  n, replace=TRUE)
        binary2  <- sample(0:1,  n, replace=TRUE)
        sex      <- sample(0:1,  n, replace=TRUE)

        return(list(    
            
            trt.coef       =trt.coef ,
            age.coef       =age.coef,
            smoke.coef     =smoke.coef,
            bmi.coef       =bmi.coef,
            covar3.coef    =covar3.coef,
            covar1.coef    =covar1.coef,
            vas.coef       =vas.coef,
            time.coef      =time.coef,
            covar2.coef    =covar2.coef,
            fact1.coef     =fact1.coef,
            binary2.coef   =binary2.coef,
            sex.coef       =sex.coef,
            
            trt= trt, 
            age=age, 
            bmi=bmi, 
            smoking=smoking,
            covar3=covar3,
            covar1=covar1, 
            vas=vas, 
            time=time, 
            covar2=covar2, 
            fact1=fact1 , 
            binary2=binary2, 
            sex=sex,
            
            randomi=randomi))
        
        
    })    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    lp1 <- reactive({
        
        d <- design()
        
        trt      <-d$trt      
        age      <-d$age       
        bmi      <-d$bmi       
        smoking  <-d$smoking  
        covar3      <-d$covar3      
        covar1   <-d$covar1   
        vas      <-d$vas       
        time     <-d$time      
        covar2   <-d$covar2    
        fact1    <-d$fact1     
        binary2 <-d$binary2  
        sex      <-d$sex  
        
        trt.coef      =d$trt.coef 
        age.coef      =d$age.coef
        smoke.coef    =d$smoke.coef
        bmi.coef      =d$bmi.coef
        covar3.coef      =d$covar3.coef
        covar1.coef   =d$covar1.coef
        vas.coef      =d$vas.coef
        time.coef     =d$time.coef
        covar2.coef   =d$covar2.coef
        fact1.coef    =d$fact1.coef
        binary2.coef =d$binary2.coef
        sex.coef      =d$sex.coef
        
        randomi <- d$randomi
        intercept <- -3
        
        
        if ( (input$Design) == "Treatment interacts with all variables" )  {
            
            lp = intercept + trt*trt.coef*(smoking*smoke.coef   +   age*age.coef  + bmi*bmi.coef + covar3*covar3.coef +
                                               covar1*covar1.coef + vas*vas.coef + time*time.coef + covar2*covar2.coef +
                                               fact1*fact1.coef +
                                               binary2*binary2.coef + sex*sex.coef) 
            
        }   else if ( (input$Design ) == "Treatment interacts with smoking only" ) {    
            
            # truth  only smoking interacts  with trt
            lp = intercept + (trt*trt.coef*smoking*smoke.coef)   +   age*age.coef   + bmi*bmi.coef + covar3*covar3.coef +
                covar1*covar1.coef + vas*vas.coef + time*time.coef + covar2*covar2.coef + fact1*fact1.coef +
                binary2*binary2.coef + sex*sex.coef
            
        }   else if ( (input$Design) == "No-interaction logit-additive model" ) {  
            
            # truth no interactions
            lp = intercept + trt*trt.coef + smoking*smoke.coef + age*age.coef  + bmi*bmi.coef + covar3*covar3.coef +
                covar1*covar1.coef + vas*vas.coef + time*time.coef + covar2*covar2.coef + fact1*fact1.coef +
                binary2*binary2.coef + sex*sex.coef
        }
        
        
        y <- ifelse(randomi < plogis(lp), 1, 0)   # one liner RANDOM!!!
        
        dat <- data.frame(cbind(y,  trt ,  smoking, age, covar3, covar1, vas, time, covar2, fact1, binary2, sex, bmi))
        
        return(list(datx=dat))
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    analysis <- reactive({
        
        da <- lp1()$datx 
        
        da$trt <-     factor(da$trt)
        da$smoking <- factor(da$smoking)
        da$fact1 <-   factor(da$fact1)
        da$binary2 <-factor(da$binary2)
        da$sex <-     factor(da$sex)
        da$bmi <-     factor(da$bmi)
        
        label(da$age)                <- 'Age'                       # label is in Hmisc
        label(da$trt)                <- 'Treatment'
        label(da$bmi)                <- 'Body Mass Index'
        label(da$smoking)            <- 'Smoking'
        label(da$covar3)             <- 'Biomarker'
        label(da$covar1)             <- 'Blood score'
        label(da$vas)                <- 'Visual analogue score'
        label(da$time)               <- 'Time since diagnosis'
        label(da$covar2)             <- 'Fitness score'
        label(da$fact1)              <- "History"
        label(da$binary2)            <- "Employed"
        label(da$sex)                <- 'Sex'
        
        dd <<- datadist(da)
        options(datadist="dd")
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # user can change the range at which effects are estimated
        
        Ages <-   (as.numeric(unlist(strsplit(input$age.range,","))))    
        
        dd$limits$age[1] <<- Ages[1]
        dd$limits$age[3] <<- Ages[2]
        
        Ages <-   (as.numeric(unlist(strsplit(input$biomarker.range,","))))    
        
        dd$limits$covar3[1] <<- Ages[1]
        dd$limits$covar3[3] <<- Ages[2]
        
        Ages <-   (as.numeric(unlist(strsplit(input$blood.range,","))))    
        
        dd$limits$covar1[1] <<- Ages[1]
        dd$limits$covar1[3] <<- Ages[2]
        
        
        Ages <-   (as.numeric(unlist(strsplit(input$vas.range,","))))    
        
        dd$limits$vas[1] <<- Ages[1]
        dd$limits$vas[3] <<- Ages[2]
        
        Ages <-   (as.numeric(unlist(strsplit(input$time.range,","))))    
        
        dd$limits$time[1] <<- Ages[1]
        dd$limits$time[3] <<- Ages[2]
        
        
        Ages <-   (as.numeric(unlist(strsplit(input$fitness.range,","))))    
        
        dd$limits$covar2[1] <<- Ages[1]
        dd$limits$covar2[3] <<- Ages[2]
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        A<-lrm(y~   trt * (smoking  + age  + bmi + covar3 + covar1 + vas + time + covar2 + fact1 + binary2 +sex),da, y=TRUE, x=TRUE)   # all interact with trt
        B<-lrm(y~  (trt *  smoking) + age  + bmi + covar3 + covar1 + vas + time + covar2 + fact1 + binary2 +sex, da, y=TRUE, x=TRUE)   # smoking * trt only
        C<-lrm(y~   trt +  smoking  + age +  bmi + covar3 + covar1 + vas + time + covar2 + fact1 + binary2 +sex, da, y=TRUE, x=TRUE)   # main effect
        
        outputx <- input$Model 
        
        if (  (outputx) == "Treatment interacts with all variables" )  {
            f <- A
            
        }   else if (  (outputx) == "Treatment interacts with smoking only" ) {
            
            f <- B
            
        }   else if (  (outputx) == "No-interaction logit-additive model" ) {
            
            f <- C
        }
        
        da$trt <- relevel(da$trt, "2")
        Aref2 <- lrm(y~   trt * (smoking  + age  + bmi + covar3 + covar1 + vas + time + covar2 + fact1 + binary2 +sex),da, y=TRUE, x=TRUE)   # all interact with trt
        da$trt <- relevel(da$trt, "3")
        Aref3 <- lrm(y~   trt * (smoking  + age  + bmi + covar3 + covar1 + vas + time + covar2 + fact1 + binary2 +sex),da, y=TRUE, x=TRUE)   # all interact with trt
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(  A=A, B=B, C=C, f=f, outputx=outputx, Aref2 = Aref2, Aref3 = Aref3)) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    output$textWithNumber2 <- renderText({ 
        
        txt <- analysis()
        HTML(paste0("Table 1. ", tags$span(style="color:black", txt$outputx  )
        ))    
        
    })  
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    output$user <- renderPrint({
        return(print(analysis()$f, digits=3))
    }) 
    
    # relative explained variation on the risk scale, see frank harrell reference
    rexv <- reactive({
        
        X <- analysis() 
        
        L1 <-   var(predict(X$B, type='fitted')) / 
            var(predict(X$A, type='fitted')) 
        
        L2 <-   var(predict(X$C, type='fitted')) / 
            var(predict(X$A, type='fitted')) 
        
        L3 <-   var(predict(X$C, type='fitted')) / 
            var(predict(X$B, type='fitted')) 
        
        AICA <-   AIC(X$A)  
        AICB <-   AIC(X$B)  
        AICC <-   AIC(X$C) 
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(  L1=L1, L2= L2, L3= L3, AICA = AICA, AICB = AICB, AICC = AICC)) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    contrastz<- reactive({
        
        X <- analysis() 
        
        A <- X$A  # trt x all
        B<-  X$B  # trt x smoke
        C <- X$C  # main
        
        k1 <- contrast(A, list(smoking=c(2), trt=c(1:3)), list(smoking=c(1), trt=c(1:3)) , fun=exp) 
        z1 <- print(k1, X=TRUE)   
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(  k1=k1,  z1=z1 )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    cont2 <- reactive({
        
        i <- as.numeric(unlist(strsplit(input$levz,",")))
        #V <-  input$predictors 
        #v <- as.character(v)
        M <- i[1] 
        N <- i[2]
        
        X <- analysis() 
        
        A <- X$A  # trt x all
        
        k99 <- contrast(A, list(smoking=M, trt=c(1:3)), list(smoking=N, trt=c(1:3)) , fun=exp) 
        z99 <- print(k99, X=TRUE)  
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(    z99=z99 )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    output$z99. <- renderPrint({
        return(cont2()$z99)
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    lrtestx<- reactive({
        
        X <- analysis() 
        
        L1 <- (lrtest(X$A, X$B))
        
        L2 <- (lrtest(X$A, X$C))
        
        L3 <- (lrtest(X$B, X$C))
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(  L1=L1, L2= L2, L3= L3)) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    output$L1a <- renderPrint({
        return(print(lrtestx()$L1, digits=3))
    })
    output$L1b <- renderPrint({
        return(print(lrtestx()$L2, digits=3))
    })
    output$L1c <- renderPrint({
        return(print(lrtestx()$L3, digits=3))
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    output$textWithNumber <- renderText({ 
        
        v1 <- as.numeric(    eval(parse(text= (input$v1)) ) )
        v2 <- as.numeric(    eval(parse(text= (input$v2)) ) )
        v3 <- as.numeric(    eval(parse(text= (input$v3)) ) )
        v4 <- as.numeric(    eval(parse(text= (input$v4)) ) )
        v5 <- as.numeric(    eval(parse(text= (input$v5)) ) )
        v6 <- as.numeric(    eval(parse(text= (input$v6)) ) )
        
        v7 <- as.numeric(    eval(parse(text= (input$v7)) ) )
        v8 <- as.numeric(    eval(parse(text= (input$v8)) ) )
        v9 <- as.numeric(    eval(parse(text= (input$v9)) ) )
        
        v10 <- as.numeric(    eval(parse(text= (input$v10)) ) )
        v11 <- as.numeric(    eval(parse(text= (input$v11)) ) )
        v12 <- as.numeric(    eval(parse(text= (input$v12)) ) )
        
        HTML(paste0( "In the case of treatment the mean treatment level 2 coefficient expectation is "
                     , tags$span(style="color:red",  p4( v1*1) ) ,
                     " and treatment level 3 coefficient expectation is "
                     , tags$span(style="color:red",  p4(v1*2) ) ,".",
                     
                     br(), br(),  
                     "For smoking, the true log odds is "
                     , tags$span(style="color:red",  p4(v3) ) , 
                     " so we expect 'smoking=2' to be "
                     , tags$span(style="color:red",  p4(v3*1) ) , 
                     " and 'smoking=3' to be "
                     , tags$span(style="color:red",  p4(v3*2) ) ,".",
                     
                     br(), br(),  
                     " In the case of age, the true effect is a change of "
                     , tags$span(style="color:red",  p4(v2) ) ,
                     " log odds for each unit change in age. Age is uniformly distributed and the minimum and maximum age are quoted in the box.","",
                     
                     br(), br(),  
                     
                     " The true coefficient for BMI a 3 level categorical factor is "
                     , tags$span(style="color:red",  p4(v4) ) ,
                     " so we expect 'BMI=2' to be "
                     , tags$span(style="color:red",  p4(v4*1) ) , 
                     " and 'BMI=3' to be "
                     , tags$span(style="color:red",  p4(v4*2) ) ,".",
                     
                     br(), br(),  
                     
                     " covar3 is a continuous variable and the true coefficient for covar3 is "
                     , tags$span(style="color:red",  p4(v5) ) ,
                     ". So for each unit change in covar3 the log odds of p(y=1|x) increases by  "
                     , tags$span(style="color:red",  p4(v5) ) ,".",
                     
                     br(), br(),  
                     "covar1 is also continuous and the true coefficient is "
                     , tags$span(style="color:red",  p4(v6) ) , 
                     ". So for each unit change in covar1 the log odds of p(y=1|x) changes by "
                     , tags$span(style="color:red",  p4(v6) ) , ".",
                     
                     br(), br(), 
                     " Vas again is continuous and the true coefficient is "
                     , tags$span(style="color:red",  p4(v7) ) , 
                     ". So for each unit change in vas the log odds of p(y=1|x) changes by "
                     , tags$span(style="color:red",  p4(v7) ) , ".",
                     
                     br(), br(), 
                     "Time is continuous and the true coefficient is "
                     , tags$span(style="color:red",  p4(v8) ) , 
                     ". So for each unit change in time the log odds of p(y=1|x) shifts by "
                     , tags$span(style="color:red",  p4(v8) ) , ".",
                     
                     br(), br(),
                     
                     "Covar2 is treated as continuous and the coefficient "
                     , tags$span(style="color:red",  p4(v9) ) , 
                     ". So for each unit change in covar2 the log odds of p(y=1|x) shifts by "
                     , tags$span(style="color:red",  p4(v9) ) , ".",
                     br(), br(),
                     "fact1, binary2 and Sex are binary predictors. For fact1 the default coefficient is "
                     , tags$span(style="color:red",  p4(v10) ) , 
                     ". So  the change to the next level of fact1 results in a "
                     , tags$span(style="color:red",  p4(v10) ) ,  
                     " shift in the log odds of p(y=1|x).",".",
                     br(), br(),
                     "For binary2, in truth there is an effect of "
                     , tags$span(style="color:red",  p4(v11) ) , 
                     ". So the change to the next level results in a "
                     , tags$span(style="color:red",  p4(v11) ) ,   
                     ". In the case of comparing sex =1 to sex=0 in truth the coefficient is "
                     , tags$span(style="color:red",  p4(v12) ) ,"."
                     
        ))    
        
    })
    
    #### pull some estimates to describe what we find
    
    output$textWithNumber3 <- renderText({ 
        
        x <-  zummary()$A1
        x2 <- zummary()$A2
        x3 <- zummary()$A3
        
        ########################
        
        txt1   <- paste0("For ",rownames(x)[1]  , " the average change in odds of the outcome comparing patients aged ",x[2,2]," 
                      with patients aged ",x[2,1] )
        
        ##treatment 1
        if ( p2(x[2,4])  > 1) { 
            
            txt2.1  <- paste0( " on treatment 1  while being identical in all other predictors is ",p2(x[2,4])," with  95%CI (",p2(x[2,6]),", ",p2(x[2,7]),"). 
                    This means that the odds of a good outcome for a patient aged ",x[2,2]," are ",p2(x[2,4])," that of the odds of a good outcome for a patient aged ",x[2,1]," if they take treatment 1. 
                    The odds (and hence probability) of a good outcome are increased for the older patient by taking the treatment 1. ")
            
            txt3.1 <- paste0( "We can also express the increase by saying that the odds are increased by approximately ",p0 ((x[2,4]-1)*100),"%.")
            
        } else { 
            
            txt2.1  <- paste0( " on treatment 1  while being identical in all other predictors is ",p2(x[2,4])," with  95%CI (",p2(x[2,6]),", ",p2(x[2,7]),").  
                    This means that the odds of a good outcome for a patient aged ",x[2,2]," are ",p2(x[2,4])," that of the odds of a good outcome for a patient aged ",x[2,1]," if they take treatment 1. 
                    The odds (and hence probability) of a good outcome are decreased for the older patient by taking the treatment 1. ")
            
            txt3.1 <- paste0( "We can also express the decrease by saying that the odds are decreased by approximately ",p0 (abs(x[2,4]-1)*100),"%,
                      since the odds are reduced by a factor of ",p2(x[2,4]),"")
            
        }
        ############################

        if ( p2(x2[2,4])  > 1) { 
            
            txt2.2  <- paste0( " on treatment 2  while being identical in all other predictors is ",p2(x2[2,4])," with  95%CI (",p2(x2[2,6]),", ",p2(x2[2,7]),"). 
                    This means that the odds of a good outcome for a patient aged ",x2[2,2]," are ",p2(x2[2,4])," that of the odds of a good outcome for a patient aged ",x2[2,1]," if they take treatment 2. 
                    The odds (and hence probability) of a good outcome are increased for the older patient by taking the treatment 2. ") 
            
            txt3.2 <- paste0( "We can also express the increase by saying that the odds are increased by approximately ",p0 ((x2[2,4]-1)*100),"%.")
            
        } else { 
            
            txt2.2  <- paste0( " on treatment 2  while being identical in all other predictors is ",p2(x2[2,4])," with  95%CI (",p2(x2[2,6]),", ",p2(x2[2,7]),"). 
                    This means that the odds of a good outcome for a patient aged ",x2[2,2]," are ",p2(x2[2,4])," that of the odds of a good outcome for a patient aged ",x2[2,1]," if they take treatment 2. 
                    The odds (and hence probability) of a good outcome are decreased for the older patient by taking the treatment 2. ")
            
            txt3.2 <- paste0( "We can also express the decrease by saying that the odds are decreased by approximately ",p0 (abs(x2[2,4]-1)*100),"%,
                      since the odds are reduced by a factor of ",p2(x2[2,4]),"")
            
        }

        ###########################
        
        if ( p2(x3[2,4])  > 1) { 
            
            txt2.3  <- paste0( " on treatment 3  while being identical in all other predictors is ",p2(x3[2,4])," with  95%CI (",p2(x3[2,6]),", ",p2(x3[2,7]),").  
                    This means that the odds of a good outcome for a patient aged ",x3[2,2]," are ",p2(x3[2,4])," that of the odds of a good outcome for a patient aged ",x3[2,1]," if they take treatment 3. 
                    The odds (and hence probability) of a good outcome are increased for the older patient by taking the treatment 3. ") 
            
            txt3.3 <- paste0( "We can also express the increase by saying that the odds are increased by approximately ",p0 ((x3[2,4]-1)*100),"%.")
            
        } else { 
            
            txt2.3  <- paste0( " on treatment 3  while being identical in all other predictors is ",p2(x3[2,4])," with  95%CI (",p2(x3[2,6]),", ",p2(x3[2,7]),").  
                    This means that the odds of a good outcome for a patient aged ",x3[2,2]," are ",p2(x3[2,4])," that of the odds of a good outcome for a patient aged ",x3[2,1]," if they take treatment 3. 
                    The odds (and hence probability) of a good outcome are decreased for the older patient by taking the treatment 3. ")
            
            txt3.3 <- paste0( "We can also express the decrease by saying that the odds are decreased by approximately ",p0 (abs(x3[2,4]-1)*100),"%,
                      since the odds are reduced by a factor of ",p2(x3[2,4]),"")
            
        }
        
        ###second example

        txt1x   <- paste0("For ",rownames(x)[3]  , " the average change in odds of the outcome comparing patients with measurement ",x[4,1]," to patients with measurement ",x[4,2] )
        
        ##treatment 1
        if ( p2(x[4,4])  > 1) { 
            
            txt2.1x  <- paste0( " on treatment 1  while being identical in all other predictors is ",p2(x[4,4])," with  95%CI (",p2(x[4,6]),", ",p2(x[4,7]),"). 
                    This means that the odds of a good outcome for a patient with measurement ",x[4,2]," are ",p2(x[4,4])," that of the odds of a good outcome for a patient with measurement ",x[4,1]," if they take treatment 1. 
                    The odds (and hence probability) of a good outcome are increased for the patient with the higher measurement by taking the treatment 1. ")
            
            txt3.1x <- paste0( "We can also express the increase by saying that the odds are increased by approximately ",p0 ((x[4,4]-1)*100),"%.")
            
        } else { 
            
            txt2.1x  <- paste0( " on treatment 1  while being identical in all other predictors is ",p2(x[4,4])," with  95%CI (",p2(x[4,6]),", ",p2(x[4,7]),").  
                    This means that the odds of a good outcome for a patient with measurement ",x[4,2]," are ",p2(x[4,4])," that of the odds of a good outcome for a patient with measurement ",x[4,1]," if they take treatment 1. 
                    The odds (and hence probability) of a good outcome are decreased for the patient with the higher measurement by taking the treatment 1. ")
            
            txt3.1x <- paste0( "We can also express the decrease by saying that the odds are decreased by approximately ",p0 (abs(x[4,4]-1)*100),"%,
                      since the odds are reduced by a factor of ",p2(x[4,4]),".")
            
        }
        ############################
        
        if ( p2(x2[4,4])  > 1) { 
            
            txt2.2x  <- paste0( " on treatment 2  while being identical in all other predictors is ",p2(x2[4,4])," with  95%CI (",p2(x2[4,6]),", ",p2(x2[4,7]),"). 
                    This means that the odds of a good outcome for a patient with measurement ",x2[4,2]," are ",p2(x2[4,4])," that of the odds of a good outcome for a patient with measurement ",x2[4,1]," if they take treatment 2. 
                    The odds (and hence probability) of a good outcome are increased for the patient with the higher measurement by taking the treatment 2. ")
            
            txt3.2x <- paste0( "We can also express the increase by saying that the odds are increased by approximately ",p0 ((x2[4,4]-1)*100),"%.")
            
        } else { 
            
            txt2.2x  <- paste0( " on treatment 2  while being identical in all other predictors is ",p2(x2[4,4])," with  95%CI (",p2(x2[4,6]),", ",p2(x2[4,7]),").  
                    This means that the odds of a good outcome for a patient with measurement ",x2[4,2]," are ",p2(x2[4,4])," that of the odds of a good outcome for a patient with measurement ",x2[4,1]," if they take treatment 2. 
                    The odds (and hence probability) of a good outcome are decreased for the patient with the higher measurement by taking the treatment 2. ")
            
            txt3.2x <- paste0( "We can also express the decrease by saying that the odds are decreased by approximately ",p0 (abs(x[4,4]-1)*100),"%,
                      since the odds are reduced by a factor of ",p2(x[4,4]),".")
            
        }
        
        ###########################
        
        if ( p2(x3[4,4])  > 1) { 
            
            txt2.3x  <- paste0( " on treatment 3  while being identical in all other predictors is ",p2(x3[4,4])," with  95%CI (",p2(x3[4,6]),", ",p2(x3[4,7]),").  
                    This means that the odds of a good outcome for a patient with measurement ",x3[4,2]," are ",p2(x3[4,4])," that of the odds of a good outcome for a patient with measurement ",x3[4,1]," if they take treatment 3. 
                    The odds (and hence probability) of a good outcome are increased for the patient with the higher measurement by taking the treatment 3. ") 
            
            txt3.3x <- paste0( "We can also express the increase by saying that the odds are increased by approximately ",p0 ((x3[4,4]-1)*100),"%.")
            
        } else { 
            
            txt2.3x  <- paste0( " on treatment 3  while being identical in all other predictors is ",p2(x3[4,4])," with  95%CI (",p2(x3[4,6]),", ",p2(x3[4,7]),").  
                    This means that the odds of a good outcome for a patient with measurement ",x3[4,2]," are ",p2(x3[4,4])," that of the odds of a good outcome for a patient with measurement ",x3[4,1]," if they take treatment 3. 
                    The odds (and hence probability) of a good outcome are decreased for the patient with the higher measurement by taking the treatment 3. ")
            
            txt3.3x <- paste0( "We can also express the decrease by saying that the odds are decreased by approximately ",p0 (abs(x3[4,4]-1)*100),"%,
                      since the odds are reduced by a factor of ",p2(x3[4,4]),".")
            
        }

        #######################
        HTML(paste0(
            br(), txt1, txt2.1 , txt3.1 ,
            br(),
            br(), txt1, txt2.2 , txt3.2 ,
            br(), 
            br(), txt1, txt2.3 , txt3.3 
            
            , br(),br(),
            br(), txt1x, txt2.1x , txt3.1x ,
            br(),
            br(), txt1x, txt2.2x , txt3.2x ,
            br(), 
            br(), txt1x, txt2.3x , txt3.3x 
            
        ) )
        
    })
    
    
    output$textWithNumber1 <- renderText({ 
        
        est <- rexv()
        
        HTML(paste0(  tags$hr(),
                      "Comparing the model with treatment x smoking interaction to 
                   the one with all covarites interacting with treatment the relative explained variation on the risk scale is "  
                      , tags$span(style="color:red",  p4( est$L1)  ),
                      " From this we see that even without penalizing for overfitting, 
                   all the treatment interactions 
                   account for only "  
                      , tags$span(style="color:red",  p4(1- est$L1) ) ,
                      " of the predictive information. The treatment x smoking 
                   model is at least "
                      , tags$span(style="color:red",  p4( est$L1)  ),
                      " adequate on a scale from 0 to 1.", 
                      br(), br(),
                      "Comparing the main effects, no-interaction logit-additive model that assumes constancy of treatment ORs to 
                   the one with all covarites interacting with treatment the relative explained variation on the risk scale is "    
                      , tags$span(style="color:red",  p4( est$L2)  ),
                      " From this we see that even without penalizing for overfitting, 
                   all the treatment interactions 
                   account for only "  
                      , tags$span(style="color:red",  p4(1- est$L2) ) ,
                      " of the predictive information. The no-interaction logit-additive 
                   model that assumes constancy of treatment ORs is at least "
                      , tags$span(style="color:red",  p4( est$L2)  ),
                      " adequate on a scale from 0 to 1. ",
                      br(), br(),
                      "Comparing the main effects, no-interaction logit-additive model that assumes constancy of treatment ORs to 
                   the one with only smoking interacting with treatment the relative explained variation on the risk scale is "    
                      , tags$span(style="color:red",  p4( est$L3)  ),
                      " From this we see that even without penalizing for overfitting, 
                   all the treatment interactions 
                   account for only "  
                      , tags$span(style="color:red",  p4(1- est$L3) ) ,
                      " of the predictive information. The no-interaction logit-additive 
                   model that assumes constancy of treatment ORs is at least "
                      , tags$span(style="color:red",  p4( est$L3)  ),
                      " adequate on a scale from 0 to 1. ",
                      tags$hr(),
                      # br(), br(),
                      "We can also use AIC to assess whether allowing for interactions will likely result in better patient-specific outcome predictions. 
                  The lower the AIC the better:",
                      br(), br(),
                      "The AIC of the  no-interaction logit-additive model that assumes constancy of treatment ORs "  
                      , tags$span(style="color:red",  p4( est$AICC)  ),
                      br(), br(),
                      "The AIC of the  treatment X smoking interaction model is "  
                      , tags$span(style="color:red",  p4( est$AICB)  ),
                      br(), br(),
                      "The AIC of the  treatment X all predictor interaction model is "  
                      , tags$span(style="color:red",  p4( est$AICA)  ),
                      br(), br(),
                      tags$hr()
                      
                      
        ))    
        
    })  
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # FULL INTERACTION MODEL
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    zummary<- reactive({
        
        X <- analysis() 
        
        #hmm why did use isolate function here?
        v0. <- isolate(as.numeric(    eval(parse(text= (input$adj.smoking   )) ) ))
        v1. <- isolate(as.numeric(    eval(parse(text= (input$adj.age       )) ) ))
        v2. <- isolate(as.numeric(    eval(parse(text= (input$adj.biomarker )) ) ))
        v3. <- isolate(as.numeric(    eval(parse(text= (input$adj.blood     )) ) ))
        v4. <- isolate(as.numeric(    eval(parse(text= (input$adj.vas       )) ) ))
        v5. <- isolate(as.numeric(    eval(parse(text= (input$adj.time      )) ) ))
        v6. <- isolate(as.numeric(    eval(parse(text= (input$adj.fitness   )) ) ))
        v7. <- isolate(as.numeric(    eval(parse(text= (input$adj.history   )) ) ))
        v8. <- isolate(as.numeric(    eval(parse(text= (input$adj.employed  )) ) ))
        v9. <- isolate(as.numeric(    eval(parse(text= (input$adj.sex       )) ) ))
        v10. <-isolate(as.numeric(    eval(parse(text= (input$adj.BMI       )) ) ))
        
        # allow the figure to change based on inputs
        v0. <-  (as.numeric(    eval(parse(text= (input$adj.smoking   )) ) ))
        v1. <-  (as.numeric(    eval(parse(text= (input$adj.age       )) ) ))
        v2. <-  (as.numeric(    eval(parse(text= (input$adj.biomarker )) ) ))
        v3. <-  (as.numeric(    eval(parse(text= (input$adj.blood     )) ) ))
        v4. <-  (as.numeric(    eval(parse(text= (input$adj.vas       )) ) ))
        v5. <-  (as.numeric(    eval(parse(text= (input$adj.time      )) ) ))
        v6. <-  (as.numeric(    eval(parse(text= (input$adj.fitness   )) ) ))
        v7. <-  (as.numeric(    eval(parse(text= (input$adj.history   )) ) ))
        v8. <-  (as.numeric(    eval(parse(text= (input$adj.employed  )) ) ))
        v9. <-  (as.numeric(    eval(parse(text= (input$adj.sex       )) ) ))
        v10. <- (as.numeric(    eval(parse(text= (input$adj.BMI       )) ) ))
        
        ##add in means of continuous vars here
        
        d <- design()

        A1 <- summary(X$A, smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                      covar2=v6., fact1=v7., binary2=v8., sex=v9., bmi=v10.,
                      trt=1, est.all=FALSE, vnames=c( "labels"))
        
        A2 <- summary(X$A, smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                      covar2=v6., fact1=v7., binary2=v8., sex=v9., bmi=v10.,
                      trt=2, est.all=FALSE, vnames=c( "labels"))
        
        A3 <- summary(X$A, smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                      covar2=v6., fact1=v7., binary2=v8., sex=v9., bmi=v10.,
                      trt=3, est.all=FALSE, vnames=c( "labels"))
       
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(  A1=A1, A2= A2, A3= A3 )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })

    output$int.trt1 <- renderPrint({
        return(print(zummary()$A1))
    }) 
    output$int.trt2 <- renderPrint({
        return(print(zummary()$A2))
    }) 
    output$int.trt3 <- renderPrint({
        return(print(zummary()$A3))
    }) 
    
    #~~~~~~~~~~~~~~~~~~~new march21~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # is this needed?
    
  zummaryx<- reactive({
        
        X <- analysis() 
        
        d <- doubleDx()
        
        k1     =d$k1
        double =d$double
        pv     =d$pv
        pvalue =d$pvalue
        M      =d$M
        N      =d$N
        M1     =d$M1
        N1     =d$N1
        v      =d$v
        
        v0. <- (as.numeric(    eval(parse(text= (input$adj.smoking   )) ) ))
        v1. <- (as.numeric(    eval(parse(text= (input$adj.age       )) ) ))
        v2. <- (as.numeric(    eval(parse(text= (input$adj.biomarker )) ) ))
        v3. <- (as.numeric(    eval(parse(text= (input$adj.blood     )) ) ))
        v4. <- (as.numeric(    eval(parse(text= (input$adj.vas       )) ) ))
        v5. <- (as.numeric(    eval(parse(text= (input$adj.time      )) ) ))
        v6. <- (as.numeric(    eval(parse(text= (input$adj.fitness   )) ) ))
        v7. <- (as.numeric(    eval(parse(text= (input$adj.history   )) ) ))
        v8. <- (as.numeric(    eval(parse(text= (input$adj.employed  )) ) ))
        v9. <- (as.numeric(    eval(parse(text= (input$adj.sex       )) ) ))
        v10. <-(as.numeric(    eval(parse(text= (input$adj.BMI       )) ) ))
         
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(    k1=k1)) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    # new march21~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # pull in stats, contrast and double differences and fed into a plot function
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$plot.trtc <- renderPlot({   

        X <- analysis() 
        
        d <- doubleDx()
        
        k1     =d$k1
        double =d$double
        pv     =d$pv
        pvalue =d$pvalue
        M      =d$M
        N      =d$N
        M1     =d$M1
        N1     =d$N1
        v      =d$v
        
        p1x <- i.plot(  v=v,  M= M,  N=N,   M1 =M1 ,  N1 = N1, k1=k1, pv=pv, double = double)  
        print(p1x)
      
    })
    
    output$int.trtc <- renderPrint({
        k1 <- zummaryx()$k1
        options(digits=6) 
        return(print(k1, X=TRUE, fun=exp, digits=6)) #*
    }) 
    output$int.trtd <- renderPrint({
        X <- analysis() 
        d <- doubleDx()
        double =d$double
        options(digits=6) # this make sure we see the statement 95 confidence limits not 9 confidence limits
        return(print(double , X=TRUE, fun=exp, digits=6)) #*
    }) 
    ## end new march21~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # SMOKING TRT INTERACTION MODEL
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    zummaryB<- reactive({
        
        X <- analysis() 
        
        A1 <- summary(X$B, smoking=1, age, covar3, covar1, vas, time, covar2, fact1, binary2, sex, bmi=1, trt=1, est.all=FALSE, vnames=c( "labels"))
        
        A2 <- summary(X$B, smoking=1,   trt=2, est.all=FALSE, vnames=c( "labels"))
        
        A3 <- summary(X$B, smoking=1,   trt=3, est.all=FALSE, vnames=c( "labels"))
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(  A1=A1, A2= A2, A3= A3)) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    output$int.trt1B <- renderPrint({
        return(print(zummaryB()$A1))
    }) 
    output$int.trt2B <- renderPrint({
        return(print(zummaryB()$A2))
    }) 
    output$int.trt3B <- renderPrint({
        return(print(zummaryB()$A3))
    }) 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # MAIN EFFECTS MODELL
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    zummaryC<- reactive({
        
        X <- analysis() 
        
        A1 <- summary(X$C, smoking=1, age, covar3, covar1, vas, time, covar2, fact1, binary2, sex, bmi=1, trt=1, est.all=FALSE, vnames=c( "labels"))
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(  A1=A1 )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    output$int.trt1C <- renderPrint({
        return(print(zummaryC()$A1))
    }) 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    output$f.plot1 <- renderPlot({   
        
        X <- analysis() 
        
        A <- X$A
        
        v0. <- as.numeric(    eval(parse(text= (input$adj.smoking)) ) )
        v1. <- as.numeric(    eval(parse(text= (input$adj.age)) ) )
        v2. <- as.numeric(    eval(parse(text= (input$adj.biomarker)) ) )
        v3. <- as.numeric(    eval(parse(text= (input$adj.blood)) ) )
        v4. <- as.numeric(    eval(parse(text= (input$adj.vas)) ) )   
        v5. <- as.numeric(    eval(parse(text= (input$adj.time)) ) ) 
        v6. <- as.numeric(    eval(parse(text= (input$adj.fitness)) ) ) 
        v7. <- as.numeric(    eval(parse(text= (input$adj.history)) ) )
        v8. <- as.numeric(    eval(parse(text= (input$adj.employed)) ) )
        v9. <- as.numeric(    eval(parse(text= (input$adj.sex)) ) )
        v10. <-as.numeric(    eval(parse(text= (input$adj.BMI)) ) )
        
        
        par(mfrow=c(1,3)) 
        
        par(oma=c(3,6,1,1)) 
        
        options(digits=1)
        
        
        plot(summary(A, smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5., 
                     covar2=v6., fact1=v7., binary2=v8., sex=v9., bmi=v10.,  
                     trt=1, est.all=FALSE, vnames=c( "labels")), 
             log=TRUE, xlim=c(log(.01),log(40)),
             q=c(  0.95 ), at=c(.02,0.05,.1,.2,.5,1,2,4,8,20), lwd=3, pch=17,
             col=   rgb(red=.4,green=.1,blue=.5,alpha=c(.5,.3,.2)),
             col.points='black', cex=1, main= "Odds Ratio (Treatment 1)", cex.main=1.8
        )
        
        plot(summary(A, smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5., 
                     covar2=v6., fact1=v7., binary2=v8., sex=v9., bmi=v10.,
                     trt=2, est.all=FALSE, vnames=c( "labels")), 
             log=TRUE, xlim=c(log(.01),log(40)),
             q=c(  0.95 ), at=c(.02,0.05,.1,.2,.5,1,2,4,8,20), lwd=3, pch=17,
             col=   rgb(red=.4,green=.1,blue=.5,alpha=c(.5,.3,.2)),
             col.points='black', cex=1, main= "Odds Ratio (Treatment 2)", cex.main=1.8
        )
        
        plot(summary(A, smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5., 
                     covar2=v6., fact1=v7., binary2=v8., sex=v9., bmi=v10.,
                     trt=3, est.all=FALSE, vnames=c( "labels")), 
             log=TRUE, xlim=c(log(.01),log(40)),
             q=c(  0.95 ), at=c(.02,0.05,.1,.2,.5,1,2,4,8,20), lwd=3, pch=17,
             col=   rgb(red=.4,green=.1,blue=.5,alpha=c(.5,.3,.2)),
             col.points='black', cex=1, main= "Odds Ratio (Treatment 3)", cex.main=1.8
        )
        
        par(mfrow=c(1,1))
        
        
        
    }) 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    output$f.plot2 <- renderPlot({
        
        X <- analysis()
        
        A <- X$B
        
        par(mfrow=c(1,3))
        
        par(oma=c(3,6,1,1))
        
        options(digits=1)
        
        plot(summary(A, smoking=1, age, covar3, covar1, vas, time, covar2, fact1, binary2, sex, bmi=1, trt=1, est.all=FALSE, vnames=c( "labels")),
             log=TRUE, xlim=c(log(.01),log(40)),
             q=c(  0.95 ), at=c(.02,0.05,.1,.2,.5,1,2,4,8,20), lwd=3, pch=17,
             col=   rgb(red=.4,green=.1,blue=.5,alpha=c(.5,.3,.2)),
             col.points='black', cex=1, main= "Odds Ratio (Treatment 1)", cex.main=1.8
        )
        
        plot(summary(A, smoking=1,   trt=2, est.all=FALSE, vnames=c( "labels")),
             log=TRUE, xlim=c(log(.01),log(40)),
             q=c(  0.95 ), at=c(.02,0.05,.1,.2,.5,1,2,4,8,20), lwd=3, pch=17,
             col=   rgb(red=.4,green=.1,blue=.5,alpha=c(.5,.3,.2)),
             col.points='black', cex=1, main= "Odds Ratio (Treatment 2)", cex.main=1.8
        )
        
        plot(summary(A, smoking=1,   trt=3, est.all=FALSE, vnames=c( "labels")),
             log=TRUE, xlim=c(log(.01),log(40)),
             q=c(  0.95 ), at=c(.02,0.05,.1,.2,.5,1,2,4,8,20), lwd=3, pch=17,
             col=   rgb(red=.4,green=.1,blue=.5,alpha=c(.5,.3,.2)),
             col.points='black', cex=1, main= "Odds Ratio (Treatment 3)", cex.main=1.8
        )
        
        par(mfrow=c(1,1))
        
        
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    output$f.plot3 <- renderPlot({
        
        X <- analysis()
        
        A <- X$C
        
        par(oma=c(3,4,1,1))
        
        options(digits=1)
        
        
        plot(summary(A, smoking=1, age, covar3, covar1, vas, time, covar2, fact1, binary2, sex, bmi=1, trt=1, est.all=FALSE, vnames=c( "labels")),
             log=TRUE, xlim=c(log(.2),log(10)),
             q=c( 0.95 ), at=c( .1,.2,.3,.5,.75,1, 1.2,1.5, 2,3,4,6,8,10), lwd=3, pch=17,
             col=   rgb(red=.4,green=.1,blue=.5,alpha=c(.5,.3,.2)),
             col.points='black', cex=1,  main=" <- worse outcomes      Odds Ratio       better outcomes ->                ", cex.main=1.8 
        )
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$datx <- renderPrint({
        return(print(lp1()$datx, digits=3))
    }) 
    output$Ax <- renderPrint({
        return(print(analysis()$A, digits=3))
    }) 
    output$Bx <- renderPrint({
        return(print(analysis()$B, digits=3))
    }) 
    output$Cx <- renderPrint({
        return(print(analysis()$C, digits=3))
    }) 
    output$Cx2 <- renderPrint({
        return(print(analysis()$C, digits=3))
    }) 
    
    ##relevel models
    output$Ax1 <- renderPrint({
        return(print(analysis()$A, digits=3))
    }) 
    output$Ax2 <- renderPrint({
        return(print(analysis()$Aref2, digits=3))
    }) 
    output$Ax3 <- renderPrint({
        return(print(analysis()$Aref3, digits=3))
    }) 
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # anova table output and gof test
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    anov <- reactive({
        
        f <- analysis()$f
        x <- anova(f, india=FALSE, vnames='labels' )
        
        
        GOF <- resid(f, "gof")   # global test of goodness of fit
        #if 'Expected value|H0'  coincidental with the 'Sum of squared errors'...don't discard model
        
        return(list( x=x, GOF= GOF)) 
        
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    output$an <- renderPrint({
        return( print(anov()$x, which=c( 'subscripts' )))
    }) 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # https://stats.stackexchange.com/questions/155246/which-variable-relative-importance-method-to-use
    output$anovaf <- renderPlot({
        
        f <- anov()$x
        
        plot(f, main=
                 "Wald χ2 statistic minus its degrees of freedom for assessing \nthe partial effect of each variable")
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    output$gofx <- renderPrint({
        return( print(anov()$GOF))
    }) 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # model checking
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # http://www.sthda.com/english/articles/36-classification-methods-essentials/148-logistic-regression-assumptions-and-diagnostics-in-r/#logistic-regression-diagnostics
    output$residz <- renderPlot({
        
        f <- analysis()$f
        
        da <- lp1()$datx 
        
        # model checking
        library(tidyverse)
        
        logits <- predict(f)
        
        mydata <- da %>%
            dplyr::select_if(is.numeric) 
        
        predictors <- colnames(mydata)
        
        mydata <- cbind(mydata, logits)
        
        #http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
        L <- gather(mydata, condition, predictor, predictors, factor_key=TRUE)
        
        # LL <- L[L$condition %in% predictors,]
        
        require(ggplot2)
        ggplot(L, aes(logits, predictor))+
            geom_point(size = 0.5, alpha = 0.5) +
            geom_smooth(method = "loess") + 
            theme_bw() + 
            facet_wrap(~condition, scales = "free_y")
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~double differences~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~see harrell ref
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~double differences~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~see harrell ref
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~double differences~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~see harrell ref
    
    
    doubleD <- reactive({
        
        X <- analysis()
        A <- X$A  # trt x all model
        da <- lp1()$datx 
        
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        v  <- input$Vars
        
        L1 <- as.numeric(    eval(parse(text= (input$treatment.level1)) ) )
        L2 <- as.numeric(    eval(parse(text= (input$treatment.level2)) ) )
        
        I1 <- as.numeric(    eval(parse(text= (input$interest.level1)) ) )
        I2 <- as.numeric(    eval(parse(text= (input$interest.level2)) ) )
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        if (v %in% "smoking") {
            
            
            double <- contrast(A, list(trt=L1,  smoking=I1),
                               list(trt=L2,  smoking=I1),
                               list(trt=L1,  smoking=I2),
                               list(trt=L2,  smoking=I2), conf.int=.95)
            
            
            w <- with(da, c(sum(smoking=I1), sum(smoking=I2)))
            
            
            pooled <- contrast(A,    list(trt=L1,  smoking=c(I1,I2)),
                               list(trt=L2,  smoking=c(I1,I2)),
                               type='average', weights=w)
            
            
        } else if (v %in% "age") {
            
            
            double <- contrast(A, list(trt=L1,  age=I1),
                               list(trt=L2,  age=I1),
                               list(trt=L1,  age=I2),
                               list(trt=L2,  age=I2), conf.int=.95)
            
            
            
            pooled <- contrast(A,    list(trt=L1,  age=c(I1,I2)),
                               list(trt=L2,  age=c(I1,I2)), type='average'
            )
            
            
        }   else if (v %in% "bmi") {
            
            
            double <- contrast(A, list(trt=L1,  bmi=I1),
                               list(trt=L2,  bmi=I1),
                               list(trt=L1,  bmi=I2),
                               list(trt=L2,  bmi=I2), conf.int=.95)
            
            
            w <- with(da, c(sum(bmi=I1), sum(bmi=I2)))
            
            
            pooled <- contrast(A,    list(trt=L1,  bmi=c(I1,I2)),
                               list(trt=L2,  bmi=c(I1,I2)),
                               type='average', weights=w)
            
            
            
        } else if (v %in% "covar3") {
            
            
            double <- contrast(A, list(trt=L1,  covar3=I1),
                               list(trt=L2,  covar3=I1),
                               list(trt=L1,  covar3=I2),
                               list(trt=L2,  covar3=I2), conf.int=.95)
            
            
            
            pooled <- contrast(A,    list(trt=L1,  covar3=c(I1,I2)),
                               list(trt=L2,  covar3=c(I1,I2)), type='average'
            )
            
            
        } else if (v %in% "covar1") {
            
            
            double <- contrast(A, list(trt=L1,  covar1=I1),
                               list(trt=L2,  covar1=I1),
                               list(trt=L1,  covar1=I2),
                               list(trt=L2,  covar1=I2), conf.int=.95)
            
            
            
            pooled <- contrast(A,    list(trt=L1,  covar1=c(I1,I2)),
                               list(trt=L2,  covar1=c(I1,I2)), type='average'
            )
            
            
            
        } else if (v %in% "vas") {
            
            
            double <- contrast(A, list(trt=L1,  vas=I1),
                               list(trt=L2,  vas=I1),
                               list(trt=L1,  vas=I2),
                               list(trt=L2,  vas=I2), conf.int=.95)
            
            
            
            pooled <- contrast(A,    list(trt=L1,  vas=c(I1,I2)),
                               list(trt=L2,  vas=c(I1,I2)), type='average'
            )
            
            
        } else if (v %in% "covar2") {
            
            
            double <- contrast(A, list(trt=L1,  covar2=I1),
                               list(trt=L2,  covar2=I1),
                               list(trt=L1,  covar2=I2),
                               list(trt=L2,  covar2=I2), conf.int=.95)
            
            
            
            pooled <- contrast(A,    list(trt=L1,  covar2=c(I1,I2)),
                               list(trt=L2,  covar2=c(I1,I2)), type='average'
            )
            
            
            
            
        } else if (v %in% "time") {
            
            
            double <- contrast(A, list(trt=L1,  time=I1),
                               list(trt=L2,  time=I1),
                               list(trt=L1,  time=I2),
                               list(trt=L2,  time=I2), conf.int=.95)
            
            
            
            pooled <- contrast(A,    list(trt=L1,  time=c(I1,I2)),
                               list(trt=L2,  time=c(I1,I2)), type='average'
            )
            
        }   else if (v %in% "sex") {
            
            
            double <- contrast(A, list(trt=L1,  sex=I1),
                               list(trt=L2,  sex=I1),
                               list(trt=L1,  sex=I2),
                               list(trt=L2,  sex=I2), conf.int=.95)
            
            
            w <- with(da, c(sum(sex=I1), sum(sex=I2)))
            
            
            pooled <- contrast(A,    list(trt=L1,  sex=c(I1,I2)),
                               list(trt=L2,  sex=c(I1,I2)),
                               type='average', weights=w)
            
            
            
        }   else if (v %in% "binary2") {
            
            
            double <- contrast(A, list(trt=L1,  binary2=I1),
                               list(trt=L2,  binary2=I1),
                               list(trt=L1,  binary2=I2),
                               list(trt=L2,  binary2=I2), conf.int=.95)
            
            
            w <- with(da, c(sum(binary2=I1), sum(binary2=I2)))
            
            
            pooled <- contrast(A,    list(trt=L1,  binary2=c(I1,I2)),
                               list(trt=L2,  binary2=c(I1,I2)),
                               type='average', weights=w)
            
            
        }   else if (v %in% "fact1") {
            
            
            double <- contrast(A, list(trt=L1,  fact1=I1),
                               list(trt=L2,  fact1=I1),
                               list(trt=L1,  fact1=I2),
                               list(trt=L2,  fact1=I2), conf.int=.95)
            
            
            w <- with(da, c(sum(fact1=I1), sum(fact1=I2)))
            
            
            pooled <- contrast(A,    list(trt=L1,  fact1=c(I1,I2)),
                               list(trt=L2,  fact1=c(I1,I2)),
                               type='average', weights=w)

        }
      
        
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(    double=double, pooled=pooled )) 
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    output$DD1 <- renderPrint({
        options(digits=6)
        return( print(doubleD()$double, digits=4))
    }) 
    output$DD2 <- renderPrint({
        return( print(doubleD()$pooled, digits=4))
    }) 
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # new function to create contrasts for interactions and double differences as well as 
    # extract the chunk interaction pvalue
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    doubleDx <- reactive({
        
        X <- analysis()
        A <- X$A  # trt x all model
        da <- lp1()$datx 
        
        v0. <- (as.numeric(    eval(parse(text= (input$adj.smoking   )) ) ))
        v1. <- (as.numeric(    eval(parse(text= (input$adj.age       )) ) ))
        v2. <- (as.numeric(    eval(parse(text= (input$adj.biomarker )) ) ))
        v3. <- (as.numeric(    eval(parse(text= (input$adj.blood     )) ) ))
        v4. <- (as.numeric(    eval(parse(text= (input$adj.vas       )) ) ))
        v5. <- (as.numeric(    eval(parse(text= (input$adj.time      )) ) ))
        v6. <- (as.numeric(    eval(parse(text= (input$adj.fitness   )) ) ))
        v7. <- (as.numeric(    eval(parse(text= (input$adj.history   )) ) ))
        v8. <- (as.numeric(    eval(parse(text= (input$adj.employed  )) ) ))
        v9. <- (as.numeric(    eval(parse(text= (input$adj.sex       )) ) ))
        v10. <-(as.numeric(    eval(parse(text= (input$adj.BMI       )) ) ))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        v  <- input$Varsx

        M <- as.numeric(unlist(strsplit(input$treatment.level1x,",")))
        N <- as.numeric(unlist(strsplit(input$treatment.level2x,",")))
        
        M1 <- as.numeric(unlist(strsplit(input$interest.level1x,",")))
        N1 <- as.numeric(unlist(strsplit(input$interest.level2x,",")))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        
        if (v %in% "smoking") {
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=c(M1,N1), age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6., binary2=v8., sex=v9., bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=c(M1,N1), age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6.,  binary2=v8., sex=v9., bmi=v10.,
                                     trt=M) 
            )
            
            
            double <- contrast(A, 
                               list(trt=N,  smoking=N1
                               ),
                               
                               list(trt=M,  smoking=N1 
                               ),
                               
                               list(trt=N,  smoking=M1
                               ),
                               
                               list(trt=M,  smoking=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* smoking", rownames(x)),"P"]
            pvalue <- x[grep("* smoking", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            
        } else if (v %in% "fact1") {
            
            k1 <- rms::contrast(A,   
                                list(fact1=c(M1,N1),  
                                     smoking=v0.,  age=v1., covar3=v2., covar1=v3., vas=v4., time =v5.,
                                     covar2=v6., binary2=v8., sex=v9. , bmi=v10., 
                                     trt=c(N)),
                                
                                list(fact1=c(M1,N1),    
                                     smoking=v0.,  age=v1., covar3=v2., covar1=v3., vas=v4., time =v5.,
                                     covar2=v6., binary2=v8., sex=v9. , bmi=v10., 
                                     trt=c(M))
                                
            )
            
            double <- contrast(A, 
                               list(trt=N,  fact1=N1
                               ),
                               
                               list(trt=M,  fact1=N1 
                               ),
                               
                               list(trt=N,  fact1=M1
                               ),
                               
                               list(trt=M,  fact1=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* fact1", rownames(x)),"P"]
            pvalue <- x[grep("* fact1", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        } else if (v %in% "age") {
            
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=c(M1,N1), covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6., binary2=v8., sex=v9., bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=c(M1,N1), covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6.,  binary2=v8., sex=v9., bmi=v10.,
                                     trt=M) 
            )
            
            
            double <- contrast(A, 
                               list(trt=N,  age=N1
                               ),
                               
                               list(trt=M,  age=N1 
                               ),
                               
                               list(trt=N,  age=M1
                               ),
                               
                               list(trt=M,  age=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* age", rownames(x)),"P"]
            pvalue <- x[grep("* age", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        }   else if (v %in% "bmi") {
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6., binary2=v8., sex=v9., bmi=c(M1,N1),
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6.,  binary2=v8., sex=v9., bmi=c(M1,N1),
                                     trt=M) 
            )
            
            
            double <- contrast(A, 
                               list(trt=N,  bmi=N1
                               ),
                               
                               list(trt=M,  bmi=N1 
                               ),
                               
                               list(trt=N,  bmi=M1
                               ),
                               
                               list(trt=M,  bmi=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* bmi", rownames(x)),"P"]
            pvalue <- x[grep("* bmi", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        } else if (v %in% "covar3") {
            
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=c(M1,N1), covar1=v3., vas=v4., time=v5.,
                                     covar2=v6., binary2=v8., sex=v9., bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=c(M1,N1), covar1=v3., vas=v4., time=v5.,
                                     covar2=v6.,  binary2=v8., sex=v9., bmi=v10.,
                                     trt=M) 
            )
            
            
            double <- contrast(A, 
                               list(trt=N,  covar3=N1
                               ),
                               
                               list(trt=M,  covar3=N1 
                               ),
                               
                               list(trt=N,  covar3=M1
                               ),
                               
                               list(trt=M,  covar3=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* covar3", rownames(x)),"P"]
            pvalue <- x[grep("* covar3", rownames(x)),]
            
            M=M;N=N;M1=M1;N1=N1;v=v;
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        } else if (v %in% "covar1") {
            
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar1=c(M1,N1), covar3=v2., vas=v4., time=v5.,
                                     covar2=v6., binary2=v8., sex=v9., bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar1=c(M1,N1), covar3=v2., vas=v4., time=v5.,
                                     covar2=v6.,  binary2=v8., sex=v9., bmi=v10.,
                                     trt=M) 
            )
            
            
            double <- contrast(A, 
                               list(trt=N,  covar1=N1
                               ),
                               
                               list(trt=M,  covar1=N1 
                               ),
                               
                               list(trt=N,  covar1=M1
                               ),
                               
                               list(trt=M,  covar1=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* covar1", rownames(x)),"P"]
            pvalue <- x[grep("* covar1", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        } else if (v %in% "vas") {
            
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=c(M1,N1), time=v5.,
                                     covar2=v6., binary2=v8., sex=v9., bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=c(M1,N1), time=v5.,
                                     covar2=v6.,  binary2=v8., sex=v9., bmi=v10.,
                                     trt=M) 
            )
            
            double <- contrast(A, 
                               list(trt=N,  vas=N1
                               ),
                               
                               list(trt=M,  vas=N1 
                               ),
                               
                               list(trt=N,  vas=M1
                               ),
                               
                               list(trt=M,  vas=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* vas", rownames(x)),"P"]
            pvalue <- x[grep("* vas", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        } else if (v %in% "covar2") {
            
            
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=c(M1,N1), binary2=v8., sex=v9., bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=c(M1,N1),  binary2=v8., sex=v9., bmi=v10.,
                                     trt=M) 
            )
            
            double <- contrast(A, 
                               list(trt=N,  covar2=N1
                               ),
                               
                               list(trt=M,  covar2=N1 
                               ),
                               
                               list(trt=N,  covar2=M1
                               ),
                               
                               list(trt=M,  covar2=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* covar2", rownames(x)),"P"]
            
            pvalue <- x[grep("* covar2", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        } else if (v %in% "time") {
            
            
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=c(M1,N1),
                                     covar2=v6., binary2=v8., sex=v9., bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=c(M1,N1),
                                     covar2=v6.,  binary2=v8., sex=v9., bmi=v10.,
                                     trt=M) 
            )
            
            
            
            double <- contrast(A, 
                               list(trt=N,  time=N1
                               ),
                               
                               list(trt=M,  time=N1 
                               ),
                               
                               list(trt=N,  time=M1
                               ),
                               
                               list(trt=M,  time=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* time", rownames(x)),"P"]
            pvalue <- x[grep("* time", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        }   else if (v %in% "sex") {
            
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6., binary2=v8., sex=c(M1,N1), bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6.,  binary2=v8., sex=c(M1,N1), bmi=v10.,
                                     trt=M) 
            )
            
            
            
            double <- contrast(A, 
                               list(trt=N,  sex=N1
                               ),
                               
                               list(trt=M,  sex=N1 
                               ),
                               
                               list(trt=N,  sex=M1
                               ),
                               
                               list(trt=M,  sex=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* sex", rownames(x)),"P"]
            pvalue <- x[grep("* sex", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        }   else if (v %in% "binary2") {
            
            
            k1 <- rms::contrast(A,
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6., binary2=c(M1,N1), sex=v9., bmi=v10.,
                                     trt=N),
                                
                                list(fact1=v7.,
                                     smoking=v0., age=v1., covar3=v2., covar1=v3., vas=v4., time=v5.,
                                     covar2=v6.,  binary2=c(M1,N1), sex=v9., bmi=v10.,
                                     trt=M) 
            )
            
            double <- contrast(A, 
                               list(trt=N,  binary2=N1
                               ),
                               
                               list(trt=M,  binary2=N1 
                               ),
                               
                               list(trt=N,  binary2=M1
                               ),
                               
                               list(trt=M,  binary2=M1 
                               ), 
                               conf.int=.95)
            
            
            # lets get the interaction p-value
            x <- anova(A, india=FALSE )
            pv <- x[grep("* binary2", rownames(x)),"P"]
            pvalue <- x[grep("* binary2", rownames(x)),]
            M=M;N=N;M1=M1;N1=N1;v=v;
            
        }
        
      
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        return(list(k1=k1, double=double, pv=pv, pvalue=pvalue, M=M, N=N, M1=M1, N1=N1,v=v))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })

    
})

# Run the application 
shinyApp(ui = ui, server = server)