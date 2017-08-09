
# code that has to run regardless of user goes here
# e.g. data assignment (though that may change if
# generalized for other users). 

library(shiny)
library(FunPhylo)
library(ggplot2)
library(dplyr)

plt.blank <-  theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(colour = "black"),
                    legend.title = (element_text(size = 9)))



shinyServer(function(input, output) {
  
  data('tyson')
  
  communities <- tyson$communities
  phylo <- tyson$phylo
  spp.list <- tyson$spp.list
  demog <- tyson$demo.data
  trait.data <- tyson$traits
  
  
  create_data <- reactive({
    out <- numeric()
    
    # if("Growth Form" %in% input$traits){
    #   traits <- c(input$traits[input$traits!="Growth Form"],
    #               "Stemmed_Herb",
    #               "Tree", "Rosette",
    #               "Vine", "SubShrub",
    #               "Shrub",
    #               "Elongated_Leafy_Rhizomatous")
    # }
    # if("Dispersal Mechanism" %in% input$traits){
    #   traits <- c(input$traits[input$traits!='Dis'],
    #               "Subterranean",
    #               "Unassisted",
    #               "Wind",
    #               "ExoZoochory",
    #               "EndoZoochory",
    #               "Ballistic",
    #               "Hoarding",
    #               "Myrmecochory",               
    #               "Water")
    # }
    # 
    # if(!"Dispersal Mechanism" %in% input$traits &
    #    !"Growth Form" %in% input$traits){
    #   traits <- input$traits
    # }
    #   

    

    for(x in unique(demog$Species)){
      phylo.mat <- make_local_phylo_dist(x, communities, phylo)
      fun.mat <- FunPhylo:::make_local_trait_dist(x, communities, trait.data,
                                                  traits = input$traits,
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

    outdat <- mutate(demog, out = out)
    
    outdat
    
  })

   output$figure1 <- renderPlot({
     dat <- create_data()
     lmdat <- summary(lm(ESCR2~ out + CRBM, data = dat))
     
     slope <- coef(lmdat)[2]
     int <- coef(lmdat)[1]
     textx <- max(dat$out)-.03
     texty <- max(dat$ESCR2) -.5
     
     Fig <- ggplot(data = dat, aes(x = out, y = ESCR2)) + 
       plt.blank + 
       geom_point() +
       geom_abline(slope = slope, intercept = int, colour = 'red') +
       annotate("text", 
                label = paste("Adjusted R^2: ", round(lmdat$adj.r.squared, 4),
                      sep = ""),
                x = textx, y = texty) +
       annotate("text",
                label = paste("Pr(>|t|) for FPD: ", 
                      round(lmdat$coefficients[2, 4], 4), sep = ""),
                      x = textx, y = texty + .3)
     
     
     print(Fig)
   })
   
   output$traits <- renderText({
     print(input$traits)
   })
   

})




