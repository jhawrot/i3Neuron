#####---------------------------------------------------------------------######
# 220310 - James Hawrot

# 3 tab shiny app for temporal proteomic data

#####---------------------------------------------------------------------######
### Load libraries and import Expression data 

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggiraph)
library(shiny)
library(shinyWidgets)
library(readr)

## import expression data and set Day as factor
Expression<-read_csv("Expression.csv") %>% 
  mutate(Day=as.factor(Day))

#####---------------------------------------------------------------------######
##### Visualize the proteomic data

### The user interface
ui <- fluidPage(
  sidebarLayout(
    position = "right",
    sidebarPanel(multiInput(
      inputId = "Genes",
      label = "Genes :", 
      choices = Expression$PG.Genes %>% unique(),
      selected = "APOE"),
      p("i3 Neuron: doxycycline-inducible human neurogenin-2"),
      p("Glutamatergic neuron differentiation"),
      p("Day 0 = iPSCs using KOLF2.1J line"),
      p("Shiny app created by James Hawrot")),
    mainPanel(tabsetPanel(type = "tabs",
                          tabPanel("Protein expression", girafeOutput("plot1")),
                          tabPanel("Percentile expression in total proteome", girafeOutput("plot2")),
                          tabPanel("Relative expression", girafeOutput("plot3"))
                          ),
              br(),
              p("Pantazis et al. 2021 : doi.org/10.1101/2021.12.15.472643"),
              p("Reilly et al. 2021 : doi.org/10.1101/2021.11.24.469921"))
  )
)
#####---------------------------------------------------------------------######
### The server
server <- function(input, output, session) { 
  output$plot1 <- renderGirafe({
    myplot1 <- girafe(code = print(Expression%>% 
                                    filter(PG.Genes %in% input$Genes) %>% 
                                    ggplot(aes(Day, measurement)) + 
                                    stat_summary(aes(group=PG.Genes,y=measurement),color="grey", fun = mean, geom = "line", position = position_dodge(width=.1))+
                                    stat_summary(aes(group=PG.Genes,y=measurement),color = "black",shape=21, position = position_dodge(width=.1))+
                                    geom_point_interactive(aes(x=Day,y=measurement_norm,fill= PG.Genes, group=PG.Genes,
                                                                tooltip=paste0(PG.Genes,":",round(measurement_norm))), 
                                                            size=4, color = "black", shape=21, position = position_dodge(width=.1))+
                                    geom_text_repel(data = Expression %>% 
                                                      filter(PG.Genes %in% input$Genes & Day=="28") %>% 
                                                      distinct(PG.Genes, .keep_all=TRUE),
                                                    aes(label = PG.Genes,y=measurement_norm),
                                                    force        = 0.5,
                                                    nudge_x      = 50,
                                                    direction    = "y",
                                                    hjust        = 1,
                                                    segment.size = 0.2,
                                                    max.overlaps = Inf,
                                                    na.rm = TRUE)+
                                    guides(color = "none") +
                                    labs(title = "Protein expression",
                                         x="Days post-differentiation",
                                         y="Protein expression intensity")+
                                    ylim(0,NA)+
                                    scale_x_discrete(expand = expansion(add = c(0.5,1.5)))+
                                    theme_classic()+
                                    theme(
                                      text = element_text(size=16),
                                      plot.title=element_text(size=18,face="bold",hjust = 0.5),
                                      legend.position="none")
    ))
    myplot1
  }
  )
  output$plot2 <- renderGirafe({
    myplot2 <- girafe(code = print(Expression%>% 
                                    filter(PG.Genes %in% input$Genes) %>% 
                                    ggplot(aes(Day, percentile_by_day_mean)) + 
                                    stat_summary(aes(group=PG.Genes,y=percentile_by_day),color="grey", fun = mean, geom = "line", position = position_dodge(width=.1))+
                                    stat_summary(aes(group=PG.Genes,y=percentile_by_day),color = "black",shape=21, position = position_dodge(width=.1))+
                                    geom_point_interactive(aes(x=Day,y=percentile_by_day_mean,fill= PG.Genes, group=PG.Genes,
                                                               tooltip=paste0(PG.Genes,":",round(percentile_by_day_mean),"th")), 
                                                           size=4, color = "black", shape=21, position = position_dodge(width=.1))+
                                    geom_text_repel(data = Expression %>% 
                                                      filter(PG.Genes %in% input$Genes & Day=="28") %>% 
                                                      distinct(PG.Genes, .keep_all=TRUE),
                                                    aes(label = PG.Genes,y=percentile_by_day_mean),
                                                    force        = 0.5,
                                                    nudge_x      = 50,
                                                    direction    = "y",
                                                    hjust        = 1,
                                                    segment.size = 0.2,
                                                    max.overlaps = Inf,
                                                    na.rm = TRUE)+
                                    guides(color = "none") +
                                    labs(title = "Percentile expression in total proteome",
                                         x="Days post-differentiation",
                                         y="Percentile by day")+
                                    ylim(0,100)+
                                    scale_x_discrete(expand = expansion(add = c(0.5,1.5)))+
                                    theme_classic()+
                                    theme(
                                      text = element_text(size=16),
                                      plot.title=element_text(size=18,face="bold",hjust = 0.5),
                                      legend.position="none")
    ))
    myplot2
  }
  )
  output$plot3 <- renderGirafe({
    myplot3 <- girafe(code = print(Expression%>% 
                                    filter(PG.Genes %in% input$Genes) %>% 
                                    ggplot(aes(Day, mean_norm)) + 
                                    stat_summary(aes(group=PG.Genes,y=normalized),color="grey", fun = mean, geom = "line", position = position_dodge(width=.1))+
                                    stat_summary(aes(group=PG.Genes,y=normalized),color = "black",shape=21, position = position_dodge(width=.1))+
                                    geom_point_interactive(aes(x=Day,y=mean_norm,fill= PG.Genes, group=PG.Genes,
                                                               tooltip=paste0(PG.Genes,":",round(mean_norm,3))), 
                                                           size=4, color = "black", shape=21, position = position_dodge(width=.1))+
                                    geom_text_repel(data = Expression %>% 
                                                      filter(PG.Genes %in% input$Genes & Day=="28") %>% 
                                                      distinct(PG.Genes, .keep_all=TRUE),
                                                    aes(label = PG.Genes,y=mean_norm),
                                                    force        = 0.5,
                                                    nudge_x      = 50,
                                                    direction    = "y",
                                                    hjust        = 1,
                                                    segment.size = 0.2,
                                                    max.overlaps = Inf,
                                                    na.rm = TRUE)+
                                    guides(color = "none") +
                                    labs(title="Relative expression",
                                         x="Days post-differentiation",
                                         y="Relative expression to day 0")+
                                    ylim(0,NA)+
                                    scale_x_discrete(expand = expansion(add = c(0.5,1.5)))+
                                    theme_classic()+
                                    theme(
                                      text = element_text(size=16),
                                      plot.title=element_text(size=18,face="bold",hjust = 0.5),
                                      legend.position="none")
    ))
    myplot3
  }
  )
}
#####---------------------------------------------------------------------######
### Run the app
shinyApp(ui, server)

