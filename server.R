
library(shiny)
library(FunPhylo)
library(ggplot2)
library(dplyr)
library(ggthemes)

data('tyson')

communities <- tyson$communities
phylo <- tyson$phylo
spp.list <- tyson$spp.list
demog <- tyson$demo.data
trait.data <- tyson$traits

shinyServer(function(input, output) {
  
  create_trait_list <- reactive({
    
    if(!"Disp_Mech" %in% input$traits &
       "Growth_Form" %in% input$traits){
      traits <- c(input$traits[input$traits!="Growth_Form"],
                  "Stemmed_Herb",
                  "Tree", "Rosette",
                  "Vine", "SubShrub",
                  "Shrub",
                  "Elongated_Leafy_Rhizomatous")
    }
    
    if("Disp_Mech" %in% input$traits &
       "Growth_Form" %in% input$traits){
      traits <- c(input$traits[-c(which(input$traits == 'Disp_Mech' | 
                               input$traits == 'Growth_Form'))],
                  "Stemmed_Herb",
                  "Tree", "Rosette",
                  "Vine", "SubShrub",
                  "Shrub",
                  "Elongated_Leafy_Rhizomatous",
                  "Subterranean",
                  "Unassisted",
                  "Wind",
                  "ExoZoochory",
                  "EndoZoochory",
                  "Ballistic",
                  "Hoarding",
                  "Myrmecochory",               
                  "Water",
                  "Clonal")
    }
    
    if("Disp_Mech" %in% input$traits &
       !"Growth_Form" %in% input$traits){
      traits <- c(input$traits[input$traits!='Disp_Mech'],
                  "Subterranean",
                  "Unassisted",
                  "Wind",
                  "ExoZoochory",
                  "EndoZoochory",
                  "Ballistic",
                  "Hoarding",
                  "Myrmecochory",               
                  "Water",
                  "Clonal")
    }
    
    
    if(!"Disp_Mech" %in% input$traits &
       !"Growth_Form" %in% input$traits){
      traits <- input$traits
    }
    traits
  })
  
  create_local_fig <- reactive({
    out <- numeric()
    
    traits <- create_trait_list()
    
    for(x in unique(demog$Species)){
      phylo.mat <- make_local_phylo_dist(x, communities, phylo)
      fun.mat <- make_local_trait_dist(x, communities, trait.data,
                                       traits = traits,
                                       scale = 'scaledBYrange')
      
      FPD <- rarefy_FPD(x, phylo.mat = phylo.mat,
                        fun.mat = fun.mat,
                        n.rare = 11, a = input$Little.a, p = 2)
      
      if("NND" %in% input$met.phylo){
        out <- c(out, as.numeric(FPD$rare.nnd))
      }
      if('MPD' %in% input$met.phylo){
        out <- c(out, as.numeric(FPD$rare.mpd))
      }
        
    }

    dat <- mutate(demog, out = out)
    
    lmdat <- summary(lm(ESCR2 ~ out + CRBM, data = dat))
    
    slope <- coef(lmdat)[2]
    int <- coef(lmdat)[1]
    textx <- max(dat$out)-.03
    texty <- max(dat$ESCR2) -.5
    
    if(input$Little.a == 0){
      x.lab <- "Functional "
    } else if(input$Little.a == 1){
      x.lab <- "Phylogenetic "
    } else{
      x.lab <- "Functional-Phylogenetic "
    }
    
    Fig <- ggplot(data = dat, aes(x = out, y = ESCR2)) + 
      theme_pander() + 
      geom_point() +
      geom_abline(slope = slope, intercept = int, colour = 'red') +
      annotate("text", 
               label = paste("Adjusted R^2: ", round(lmdat$adj.r.squared, 4),
                             sep = ""),
               x = textx, y = texty) +
      annotate("text",
               label = paste("Pr(>|t|) for FPD: ", 
                             round(lmdat$coefficients[2, 4], 4), sep = ""),
               x = textx, y = texty + .3) +
      scale_x_continuous(paste(x.lab, input$met.phylo, sep = "")) +
      scale_y_continuous("Strength of Competitive Interactions")
    
    
    Fig
    
  })

  create_lin_regional_fig <- reactive({
    out <- numeric()
    traits <- create_trait_list()
    
    
    spp.list$Species <- gsub("-", "\\.", spp.list$Species)
    trait.data <- trait.data[trait.data$Species.Name %in% spp.list$Species, ]
    trait.data$Species.Name <- gsub("-", "\\.", trait.data$Species.Name)
    
    phylo.mat <- make_regional_phylo_dist(trait.data$Species.Name, phylo = phylo)
    fun.mat <- make_regional_trait_dist(trait.data, traits) %>% as.matrix %>%
               as.data.frame()
    
    phylo.mat <- phylo.mat[rownames(phylo.mat) %in% sort(rownames(fun.mat)),
                           names(phylo.mat) %in% sort(names(fun.mat))] %>%
      .[sort(rownames(.)), sort(names(.))]
    fun.mat <- fun.mat[sort(rownames(fun.mat)), sort(names(fun.mat))]
    
    FPD <- func_phy_dist(FDist = as.matrix(fun.mat),
                         PDist = as.matrix(phylo.mat),
                         phyloWeight = input$Little.a,
                         p = 2) %>% data.frame()
    diag(FPD) <- NA
    
    for(x in unique(demog$Species)){
      if("MPD" %in% input$met.phylo){
        out <- c(out, mean(FPD[ ,x], na.rm = TRUE))
      } else{
        out <- c(out, min(FPD[ ,x], na.rm = TRUE))
      }
    }
    
    demog <- mutate(demog, out = out)
    
    lmdat <- summary(lm(ESCR2 ~ out + CRBM, data = demog))
    
    slope <- coef(lmdat)[2]
    int <- coef(lmdat)[1]
    textx <- max(demog$out)-.03
    texty <- max(demog$ESCR2) -.5
    
    if(input$Little.a == 0){
      x.lab <- "Functional "
    } else if(input$Little.a == 1){
      x.lab <- "Phylogenetic "
    } else{
      x.lab <- "Functional-Phylogenetic "
    }
    
    Fig <- ggplot(demog, aes(x = out, y = ESCR2)) + 
           theme_pander() +
           geom_point() + 
           geom_abline(slope = slope, intercept = int, colour = 'red') +
           annotate("text", 
                    label = paste("Adjusted R^2: ",
                                  round(lmdat$adj.r.squared, 4),
                                  sep = ""),
                    x = textx, y = texty) +
           annotate("text",
                    label = paste("Pr(>|t|) for FPD: ", 
                                  round(lmdat$coefficients[2, 4], 4),
                                  sep = ""),
                    x = textx, y = texty + .3) +
           scale_x_continuous(paste(x.lab, input$met.phylo, sep = "")) +
           scale_y_continuous("Strength of Competitive Interactions")
    
    Fig
  })
  
  create_log_regional_fig <- reactive({
    out <- numeric()
    traits <- create_trait_list()
    
    
    Exotics <- dplyr::filter(spp.list, Exotic == 1 & Monocot == 0)
    Exotics$Species <- gsub("-", "\\.", Exotics$Species)
    trait.data <- trait.data[trait.data$Species.Name %in% Exotics$Species, ]
    trait.data$Species.Name <- gsub("-", "\\.", trait.data$Species.Name)
    
    if(input$Little.a != 1){
      Exotics <- Exotics[Exotics$Species %in% trait.data$Species.Name, ]
      # warning('Sample sizes are considerably smaller when incorporating\n',
      #         'functional trait information for expert classification regressions\n',
      #         'as our data set does not contain functional trait information\n',
      #         'for all species at Tyson Research Center. Therefore, these results\n',
      #         'may not be directly comparable to results from regressions based\n',
      #         'purely on phylogenetic information.')
      # 
      fun.mat <- make_regional_trait_dist(trait.data, traits)
      phylo.mat <- make_regional_phylo_dist(Exotics$Species, phylo)
      
      FPD <- func_phy_dist(FDist = as.matrix(fun.mat),
                           PDist = as.matrix(phylo.mat),
                           phyloWeight = input$Little.a,
                           p = 2)
      
      diag(FPD) <- NA
      
      for(x in unique(Exotics$Species)){
        if("MPD" %in% input$met.phylo){
          Exotics[Exotics$Species == x, 'out'] <- mean(FPD[ ,x], na.rm = T)
        } else {
          Exotics[Exotics$Species == x, 'out'] <- min(FPD[ ,x], na.rm = T)
        }
        
      }
      
      x.lab <- "Function-Phylogenetic "
      
    } else {
    
      phylo.mat <- make_regional_phylo_dist(Exotics$Species, phylo)
      phylo.mat <- phylo.mat/max(phylo.mat)
      diag(phylo.mat) <- NA
      for(x in unique(Exotics$Species)){
        if("MPD" %in% input$met.phylo){
          out <- c(out, mean(phylo.mat[ ,x], na.rm = T))
        } else {
          out <- c(out, min(phylo.mat[ ,x], na.rm = T))
        }
        
      }
      x.lab <- "Phylogenetic "
    }
    
    Exotics <- mutate(Exotics, out = out)
    n <- dim(Exotics)[1]
    Exotics$Invasive <- as.numeric(Exotics$Invasive)
    
    if(input$Little.a == 0) x.lab <- "Functional "
    
    lmdat <- summary(glm(as.numeric(Invasive) ~ out, family = binomial(), data = Exotics))
    
    textx <- max(Exotics$out)-.03
    texty <- .1
    
    Fig <- ggplot(data = Exotics, aes(x = out, y = Invasive)) +
           theme_pander() +
           geom_point() + 
           stat_smooth(formula = y ~ x,
                       method = "glm", method.args = list(family = "binomial"),
                       se = FALSE, color = 'red') +
           annotate('text', label = paste("Pr(>|t|):", 
                                          round(coef(lmdat)[2,4], 3)),
                    x = textx, y = texty, size = 4) +
           annotate('text', label = paste("Sample Size: ", n, sep = ""),
                    x = textx, y = texty - .05) +
           scale_x_continuous(paste(x.lab, input$met.phylo, sep = "")) + 
           scale_y_continuous("Pr(Invasive)")
    
  })
  
  regression_switch <- reactive({
    metScaleSwitch <- paste(input$met.inv, input$scale, sep = "_")
    x <- switch(metScaleSwitch,
           "lambda_reg" = create_lin_regional_fig(),
           'lambda_loc' = create_local_fig(),
           'mepp_reg' = create_log_regional_fig(),
           'mepp_loc' = stop('Sample size too small to compute a logistic \n',
                             'regression for focal species only. Please use\n',
                             'regional spatial scale for expert classification\n',
                             'results.'))
    x
    
  })
  
  output$figure1 <- renderPlot({
    
    Fig <- regression_switch()
    
    print(Fig)
      
  })
  
})




