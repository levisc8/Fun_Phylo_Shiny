
library(shiny)
library(FunPhylo)
library(ggplot2)
library(dplyr)
library(ggthemes)

data('tyson')

communities <- tyson$communities
spp.list <- tyson$spp.list
demog <- tyson$demo.data
trait.data <- tyson$traits
demog$Invasive <- NA_character_
for(i in seq_len(dim(demog)[1])) {
  demog$Invasive[i] <- communities$Invasive[communities$community == demog$Species[i]][1]
}

demog$Invasive <- ifelse(demog$Invasive == 1, "Yes", "No")

shinyServer(function(input, output) {
  
  choose_phylo <- reactive(
    switch(input$phylo,
           'tank' = tyson$phylo,
           'all_otb' = tyson$phylo_all,
           'gb_otb' = tyson$phylo_gb)
    
  )
  
  # Creates list of traits based on inputs (very inefficiently!)
  # Growth Form and dispersal mechanism are actually 
  # a combination of many dummy variables that describe the levels
  # they can take
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
  
  # reactive to create local ESCR~Novelty+CRBM data
  create_local_fig <- reactive({
    out <- numeric()
    
    traits <- create_trait_list()
    
    phylo <- choose_phylo()
    
    for(x in unique(demog$Species)){
      
      
      use_com <- filter(communities, exotic_species == x &
                          community %in% phylo$tip.label)
      
      phylo.mat <- make_local_phylo_dist(x, use_com, phylo)
      fun.mat <- make_local_trait_dist(x, use_com, trait.data,
                                       traits = traits,
                                       scale = 'scaledBYrange')
      
      FPD <- rarefy_FPD(x, phylo.mat = phylo.mat,
                        fun.mat = fun.mat,
                        metric = input$met.phylo,
                        n.rare = 11, a = input$Little.a, p = 2,
                        abundance.weighted = input$AW,
                        community.data = communities)
      
      # make sure we save the requested metric
      
      if('NND' %in% input$met.phylo){
        out <- c(out, ifelse(input$log,
                             log(as.numeric(FPD$rare.nnd)),
                             as.numeric(FPD$rare.nnd)))
      }
      if('MPD' %in% input$met.phylo){
        out <- c(out, ifelse(input$log,
                             log(as.numeric(FPD$rare.mpd)),
                             as.numeric(FPD$rare.mpd)))
      }
      
    }
    # create data frame and formula for model
    dat <- mutate(demog, out = out)
    
    if(input$resp.var == 'sig'){
      lm.form <- as.formula(paste0('ESCR2 ~ out + CRBM'))  
    } else {
      lm.form <- as.formula(paste0('ESCR ~ out + CRBM'))  
    }
    
    lm_form <- lm(lm.form, data = dat)
    lmdat <- summary(lm_form)
    
    # extract data needed for plotting
    
    sim_out <- seq(min(dat$out), max(dat$out), length.out = 14)
    sim_crb <- seq(min(dat$CRBM), max(dat$CRBM), length.out = 14)
    
    pred <- predict(lm_form,
                    data.frame(out = sim_out, CRBM = sim_crb),
                    type = 'response',
                    interval = 'confidence',
                    se.fit = TRUE)$fit %>%
      data.frame %>% 
      cbind(., sim_out)
    
    textx <- max(dat$out) - ((max(dat$out) - min(dat$out)) / 8) 
    texty <- max(dat$ESCR2) -.5
    
    # if(input$log){
    #   textx <- textx - (1.1-input$Little.a * .1)
    # }
    
    # X-axis label changes depending on little a
    if(input$Little.a == 0){
      x.lab <- "Functional "
    } else if(input$Little.a == 1){
      x.lab <- "Phylogenetic "
    } else{
      x.lab <- "Functional-Phylogenetic "
    }
    # Plot!
    Fig <- ggplot(data = dat, aes(x = out, y = ESCR2)) + 
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(size = 20,
                                        margin = margin(t = 0,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.title.x = element_text(size = 20, 
                                        margin = margin(t = 20,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.line = element_line(size = 1.5),
            axis.text = element_text(size = 16)) + 
      geom_point(aes(color = Invasive),
                 size = 2) +
      geom_line(data = pred,
                aes(x = sim_out,
                    y = fit),
                colour = 'red') +
      annotate("text", 
               label = paste("Adjusted R^2: ", round(lmdat$adj.r.squared, 4),
                             sep = ""),
               x = textx, y = texty) +
      annotate("text",
               label = paste("Pr(>|t|) for FPD: ", 
                             round(lmdat$coefficients[2, 4], 4), sep = ""),
               x = textx, y = texty + .3) +
      scale_x_continuous(paste(x.lab, input$met.phylo, sep = "")) +
      scale_y_continuous(expression(frac(ln(lambda[CR] + 0.5), ln(lambda[Control] + 0.5))))
    
    
    Fig
    
  })
  
  create_lin_regional_fig <- reactive({
    out <- numeric()
    traits <- create_trait_list()
    phylo <- choose_phylo()
    
    # make sure that trait names match phylogeny names
    spp.list$Species <- gsub("-", "\\.", spp.list$Species)
    trait.data <- trait.data[trait.data$Species.Name %in% spp.list$Species, ]
    trait.data$Species.Name <- gsub("-", "\\.", trait.data$Species.Name)
    
    phylo.mat <- make_regional_phylo_dist(trait.data$Species.Name, phylo = phylo) %>%
      as.matrix() %>%
      as.data.frame()
    
    fun.mat <- make_regional_trait_dist(trait.data, traits) %>% 
      as.matrix %>%
      as.data.frame()
    
    # make sure phylo distances match functional distances
    phylo.mat <- phylo.mat[rownames(phylo.mat) %in% sort(rownames(fun.mat)),
                           names(phylo.mat) %in% sort(names(fun.mat))] %>%
      .[sort(rownames(.)), sort(names(.))]
    fun.mat <- fun.mat[sort(rownames(fun.mat)), sort(names(fun.mat))]
    
    FPD <- func_phy_dist(FDist = as.matrix(fun.mat),
                         PDist = as.matrix(phylo.mat),
                         phyloWeight = input$Little.a,
                         p = 2) %>% data.frame()
    diag(FPD) <- NA
    # Store correct metric
    for(x in unique(demog$Species)){
      if("MPD" %in% input$met.phylo){
        out <- c(out, mean(FPD[ ,x], na.rm = TRUE))
      } else{
        out <- c(out, min(FPD[ ,x], na.rm = TRUE))
      }
    }
    
    if(input$resp.var == 'sig'){
      lm.form <- as.formula(paste0('ESCR2 ~ out + CRBM'))  
    } else {
      lm.form <- as.formula(paste0('ESCR ~ out + CRBM'))  
    }
    
    # create data frame and models
    demog <- dat <- mutate(demog, out = out)
    
    lm_form <- lm(lm.form, data = dat)
    lmdat <- summary(lm_form)
    
    # extract data needed for plotting
    
    sim_out <- seq(min(dat$out), max(dat$out), length.out = 14)
    sim_crb <- seq(min(dat$CRBM), max(dat$CRBM), length.out = 14)
    
    pred <- predict(lm_form,
                    data.frame(out = sim_out, 
                               CRBM = sim_crb),
                    type = 'response',
                    interval = 'confidence',
                    se.fit = TRUE)$fit %>%
      data.frame %>% 
      cbind(., sim_out)
    
    textx <- max(demog$out)-.03
    texty <- max(demog$ESCR2) -.5
    # x-axis label changes depending on little a
    if(input$Little.a == 0){
      x.lab <- "Functional "
    } else if(input$Little.a == 1){
      x.lab <- "Phylogenetic "
    } else{
      x.lab <- "Functional-Phylogenetic "
    }
    # plot!
    Fig <- ggplot(demog, aes(x = out, y = ESCR2)) + 
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(size = 20,
                                        margin = margin(t = 0,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.title.x = element_text(size = 20, 
                                        margin = margin(t = 20,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.line = element_line(size = 1.5),
            axis.text = element_text(size = 16))  +
      geom_point(aes(color = Invasive)) + 
      geom_line(data = pred,
                aes(x = sim_out,
                    y = fit),
                colour = 'red') +
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
      scale_y_continuous(expression(frac(ln(lambda[CR] + 0.5),
                                         ln(lambda[Control] + 0.5))))
    
    Fig
  })
  
  create_log_regional_fig <- reactive({
    out <- numeric()
    traits <- create_trait_list()
    phylo <- choose_phylo()
    
    
    Exotics <- dplyr::filter(spp.list, 
                             Exotic == 1 &
                               Monocot == 0 &
                               Species != 'Cerastium_spp.')
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
      
      phylo.mat <- cophenetic(phylo) %>% data.frame()
      # phylo.mat <- phylo.mat/max(phylo.mat)
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
    n <- length(unique(Exotics$Species))
    Exotics$Invasive <- ifelse(Exotics$Invasive == 1, "Yes", "No")
    
    if(input$Little.a == 0) x.lab <- "Functional "
    
    lmdat <- summary(glm(as.numeric(Invasive) ~ out,
                         family = binomial(),
                         data = Exotics))
    
    textx <- max(Exotics$out)-.03
    texty <- .1
    
    Fig <- ggplot(data = Exotics, aes(x = out, y = Invasive)) +
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(size = 20,
                                        margin = margin(t = 0,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.title.x = element_text(size = 20, 
                                        margin = margin(t = 20,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.line = element_line(size = 1.5),
            axis.text = element_text(size = 16))  +
      geom_point(aes(color = Invasive)) + 
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
  
  create_regional_r2_a_plot <- reactive({
    
    traits <- create_trait_list()
    phylo <- choose_phylo()
    
    a_seq <- seq(0, 1, length.out = input$res)
    mod.data <- demog
    
    mod.data[ , paste0('nna_', a_seq)] <- NA
    mod.data[ , paste0('mpa_', a_seq)] <- NA
    
    R2dat <- data.frame(A = a_seq,
                        NND = rep(NA, input$res),
                        MPD = rep(NA, input$res))
    
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
    
    
    for(a in a_seq){
      i <- which(a_seq == a)
      FPD <- func_phy_dist(FDist = fun.mat, 
                           PDist = phylo.mat,
                           phyloWeight = a,
                           p = 2) %>% data.frame()
      diag(FPD) <- NA
      for(x in unique(demog$Species)){
        mod.data[mod.data$Species == x, paste0('nna_', a)] <- min(FPD[ ,x],
                                                                  na.rm = T)
        mod.data[mod.data$Species == x, paste0('mpa_', a)] <- mean(FPD[ ,x],
                                                                   na.rm = T)
      }
      
      
      if(input$resp.var == 'sig'){
        nnd.form <- as.formula(paste0('ESCR2 ~ nna_', a,'+ CRBM'))
        mpd.form <- as.formula(paste0('ESCR2 ~ mpa_', a,'+ CRBM'))
      } else {
        nnd.form <- as.formula(paste0('ESCR ~ nna_', a,'+ CRBM'))
        mpd.form <- as.formula(paste0('ESCR ~ mpa_', a,'+ CRBM'))
      }
      
      nnd.form <- as.formula(paste0('ESCR2 ~ nna_', a,'+ CRBM'))
      mpd.form <- as.formula(paste0('ESCR2 ~ mpa_', a,'+ CRBM'))
      
      R2dat$NND[i] <- r2_calc(mod.data, nnd.form)
      R2dat$MPD[i] <- r2_calc(mod.data, mpd.form)
    }
    
    maxr2 <- max(R2dat[ ,2:3])
    if(maxr2 %in% R2dat$MPD) {
      maxr2met <- 'MPD'
      maxr2A <- R2dat[which(R2dat$MPD == maxr2), 'A']
    } else if(maxr2 %in% R2dat$NND) {
      maxr2met <- 'NND'
      maxr2A <- R2dat[which(R2dat$NND == maxr2), 'A']
      
    }
    maxr2A <- round(maxr2A, 4)
    
    Fig <- ggplot(data = R2dat, aes(x = A)) +
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(size = 20,
                                        margin = margin(t = 0,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.title.x = element_text(size = 20, 
                                        margin = margin(t = 20,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.line = element_line(size = 1.5),
            axis.text = element_text(size = 16))  +
      geom_point(aes(y = NND, color = 'NND')) +
      geom_line(aes(y = NND, color = 'NND')) + 
      geom_point(aes(y = MPD, color = 'MPD')) +
      geom_line(aes(y = MPD, color = 'MPD')) +
      scale_x_continuous('Phylogenetic Scaling Parameter', limits = c(0,1)) +
      scale_y_continuous(expression(paste('Adjusted ', R^2)),
                         limits = c(-0.2 , 1.2),
                         breaks = seq(0, 1, 0.2)) +
      scale_color_manual('Metric', 
                         breaks = c("NND", "MPD"),
                         values = c("red", "blue")) + 
      annotate('text',
               label = paste0('Maximum R^2: ', round(maxr2, 3)),
               x = .9, y = .95) +
      annotate('text', 
               label = paste0('Best a-value: ', maxr2A),
               x = .9, y = .9) +
      annotate('text', 
               label = paste0('Best Performing Metric: ', maxr2met),
               x = .9, y = .85)
    
  })
  
  create_local_r2_a_plot <- reactive({
    traits <- create_trait_list()
    phylo <- choose_phylo()
    
    a_seq <- seq(0, 1, length.out = input$res)
    mod.data <- demog
    
    mod.data[ , paste0('nna_', a_seq)] <- NA
    mod.data[ , paste0('mpa_', a_seq)] <- NA
    
    R2dat <- data.frame(A = a_seq,
                        NND = rep(NA, input$res),
                        MPD = rep(NA, input$res))
    
    
    for(x in unique(demog$Species)){
      
      
      use_com <- filter(communities, exotic_species == x &
                          community %in% phylo$tip.label)
      
      phylo.mat <- make_local_phylo_dist(x, use_com, phylo)
      fun.mat <- make_local_trait_dist(x, use_com, trait.data,
                                       traits = traits,
                                       scale = 'scaledBYrange')
      
      for(a in a_seq){
        FPD <- rarefy_FPD(x, phylo.mat = phylo.mat,
                          fun.mat = fun.mat,
                          n.rare = 11, a = a, p = 2,
                          abundance.weighted = input$AW,
                          community.data = communities)
        
        mod.data[mod.data$Species == x, paste0('nna_', a)] <- ifelse(input$log, 
                                                                     log(FPD$rare.nnd),
                                                                     FPD$rare.nnd)
        mod.data[mod.data$Species == x, paste0('mpa_', a)] <- ifelse(input$log,
                                                                     log(FPD$rare.mpd),
                                                                     FPD$rare.mpd)
        
      }
    }
    
    for(a in a_seq){
      i <- which(a_seq == a)
      
      if(input$resp.var == 'sig'){
        nnd.form <- as.formula(paste0('ESCR2 ~ nna_', a,'+ CRBM'))
        mpd.form <- as.formula(paste0('ESCR2 ~ mpa_', a,'+ CRBM'))
      } else {
        nnd.form <- as.formula(paste0('ESCR ~ nna_', a,'+ CRBM'))
        mpd.form <- as.formula(paste0('ESCR ~ mpa_', a,'+ CRBM'))
      }
      
      R2dat$NND[i] <- r2_calc(mod.data, nnd.form)
      R2dat$MPD[i] <- r2_calc(mod.data, mpd.form) 
    }
    
    maxr2 <- max(R2dat[ ,2:3])
    if(maxr2 %in% R2dat$MPD) {
      maxr2met <- 'MPD'
      maxr2A <- R2dat[which(R2dat$MPD == maxr2), 'A']
    } else if(maxr2 %in% R2dat$NND) {
      maxr2met <- 'NND'
      maxr2A <- R2dat[which(R2dat$NND == maxr2), 'A']
      
    }
    maxr2A <- round(maxr2A, 4)
    
    Fig <- ggplot(data = R2dat, aes(x = A)) +
      theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(size = 20,
                                        margin = margin(t = 0,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.title.x = element_text(size = 20, 
                                        margin = margin(t = 20,
                                                        l = 0,
                                                        b = 0,
                                                        r = 20)),
            axis.line = element_line(size = 1.5),
            axis.text = element_text(size = 16))  +
      geom_point(aes(y = NND, color = 'NND')) +
      geom_line(aes(y = NND, color = 'NND')) + 
      geom_point(aes(y = MPD, color = 'MPD')) +
      geom_line(aes(y = MPD, color = 'MPD')) +
      scale_x_continuous('Phylogenetic Scaling Parameter', limits = c(0,1)) +
      scale_y_continuous(expression(paste('Adjusted ', R^2)),
                         limits = c(-0.2 , 1.2),
                         breaks = seq(0, 1, 0.2)) + 
      scale_color_manual('Metric', 
                         breaks = c("NND", "MPD"),
                         values = c("red", "blue")) + 
      annotate('text',
               label = paste0('Maximum R^2: ', round(maxr2, 3)),
               x = .9, y = .95) +
      annotate('text', 
               label = paste0('Best a-value: ', maxr2A),
               x = .9, y = .9) +
      annotate('text', 
               label = paste0('Best Performing Metric: ', maxr2met),
               x = .9, y = .85)
    
    
    
  })
  
  UI_Input_switch <- reactive({
    x <- switch(input$plot,
                'lil.a' = r2_a_switch(),
                'FS' = FS_switch())
    
    x
  })
  
  
  r2_a_switch <- reactive({
    scale <- input$scale
    x <- switch(scale,
                'reg' = create_regional_r2_a_plot(),
                'loc' = create_local_r2_a_plot())
    
    x
  })
  
  # Invasive classification or log response ratio
  FS_switch <- reactive({
    metScaleSwitch <- paste(input$met.inv, input$scale, sep = "_")
    x <- switch(metScaleSwitch,
                "lambda_reg" = create_lin_regional_fig(),
                'lambda_loc' = create_local_fig(),
                'mepp_reg' = create_log_regional_fig(),
                'mepp_loc' = stop('Feature not yet added, but is on the way\n',
                                  'Sorry for the inconvenience!'))
    x
    
  })
  
  output$figure1 <- renderPlot({
    
    Fig <- UI_Input_switch()
    
    print(Fig)
    
  })
  
  # Uncomment to render Demography data summary table
  # output$traits <- renderText({
  #   
  #    print(input$traits)
  # })
  
})




