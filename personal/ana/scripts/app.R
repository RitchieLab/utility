library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)
library(data.table)
library(dplyr)
library(ggiraph)
library(ggforce)
library(ggrepel)
library(dashboardthemes)

age = read.delim("PMBB_AGE_SEX.txt", stringsAsFactors = FALSE)
age$AGE <- as.numeric(age$AGE)
eth <- read.delim("PMBB_RACE_N.txt", stringsAsFactors = FALSE)
eth$PCT <- ifelse(eth$SUBJ_GROUP=="PMBB", eth$N/53611, eth$N/18981)
tti <- read.delim("Top_10_ICD.txt", stringsAsFactors = FALSE)
cat <- read.delim("PMBB_ICD_CAT_N.txt", stringsAsFactors = FALSE)
cat$Category[cat$Category=="NULL"] <- "Other"
lab <- as.data.frame(data.table::fread("PMBB_LABS.txt"))
lab$RESULT_VALUE_NUM <- as.numeric(lab$RESULT_VALUE_NUM)
l <- unique(lab$DATASET)
names(l) <- unique(lab$DATASET)
rec <- read.delim("Department_Recruitment_dummy.txt")

ui <- dashboardPage(
  dashboardHeader(title="PMBB Demographics"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    shinyDashboardThemes(
      theme = "grey_dark"
    ),
    # Boxes need to be put in a row (or column)
    fluidPage(
      fluidRow(
        column(width=4,
               box(
                 plotOutput("plot1", height = 260), 
                 title="Age Distribution",
                 solidHeader = FALSE,
                 status = "info",
                 width=NULL
               ),
               box(
                 plotOutput("plot3", height = 260), 
                 #plotOutput("plot3", height = 260),
                 title="Top 10 ICD-9 Codes",
                 solidHeader = FALSE,
                 status = "info",
                 width=NULL
               )
        ),
        column(width=4,
               box(
                 plotOutput("plot2", height = 260), 
                 title="Self-Reported Ancestry",
                 solidHeader = FALSE,
                 status = "info",
                 width=NULL
               ),
               box(
                 plotOutput("plot4", height=260),
                 title = "Department Recruitment",
                 solidHeader = FALSE,
                 status="info",
                 width=NULL
               )
               
        ),
        column(width = 4,
               box(
                 title = "Subject Group",
                 solidHeader = TRUE,
                 status = "primary",
                 selectInput("select", "Subject Group:", c("All PMBB" = "PMBB", "Genotyped" = "GENOTYPED")),
                 width=NULL
               ),
               box(
                 dataTableOutput("table", height = 250), 
                 title="Number of Patients per ICD Category",
                 solidHeader = FALSE,
                 status = "warning",
                 width=NULL
               )
        )
      ),
      fluidRow(
        column(width=8,
               box(
                 plotOutput("plot5", height = 260),
                 title="Clinical Labs",
                 solidHeader=FALSE,
                 status="info",
                 width=NULL
               )),
        column(width=4,
               box(
                 title = "Lab",
                 solidHeader = TRUE,
                 status = "primary",
                 selectInput("select_lab", "Lab:", l),
                 width=NULL
               )
        )
      )
    )
  )
)

server <- function(input, output) { 
  
  
  output$plot1 <- renderPlot({
    if(input$select == "GENOTYPED"){
      plot_age <- age[age$SUBJ_GROUP==input$select,]
    } else {
      plot_age <- age
    }
    ggplot(data=plot_age, aes(x=AGE, color=GENDER_CODE)) + geom_density() + theme_minimal() + 
      scale_color_brewer(name='Gender', palette = "Dark2") + xlab("Patient Age (Years)") +
      theme(panel.background = element_rect(fill = "#343E48",color="#ffffff"), 
            plot.background = element_rect(fill = "#343E48", color = NA,),
            axis.text = element_text(color="#FFFFFF"),
            axis.title = element_text(color="#ffffff"),
            legend.text = element_text(color="#ffffff"),
            legend.title = element_text(color="#ffffff"),
            panel.grid = element_line(color="#303030", size = 0.3))
    
  })
  
  output$plot2 <- renderPlot({
    plot_eth <- eth[eth$SUBJ_GROUP==input$select,]
    if(input$select=="PMBB"){
      levs=rev(c("WHITE", "BLACK", "UNKNOWN", "OTHER", "ASIAN", "AM IND AK NATIVE", "HI PAC ISLAND"))
    } else {
      levs=rev(c("WHITE", "BLACK", "UNKNOWN", "OTHER", "ASIAN", "HI PAC ISLAND","AM IND AK NATIVE"))
    }
    
    plot_eth$Label <- paste0(plot_eth$RACE_CODE, " - ", signif(plot_eth$PCT*100, digits = 2), "%")
    plot_eth <- plot_eth %>%
      mutate(end = 2 * pi * cumsum(N)/sum(N),
             start = lag(end, default=0),
             middle = 0.5 * (start + end),
             hjust = ifelse(middle > pi, 1, 0),
             vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
    ggplot(plot_eth) + 
      geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                       start = start, end = end, fill = factor(RACE_CODE, levels=rev(levs))), color="white") +
      geom_text_repel(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = Label,
                          hjust = hjust, vjust = vjust), color="#ffffff") +
      coord_fixed() +
      scale_x_continuous(limits = c(-1.3, 1.4),  # Adjust so labels are not cut off
                         name = "", breaks = NULL, labels = NULL) +
      scale_y_continuous(limits = c(-1, 1.2),      # Adjust so labels are not cut off
                         name = "", breaks = NULL, labels = NULL) +
      theme(legend.position = "none", panel.background = element_blank()) +
      scale_fill_brewer(palette="Dark2") +
      theme(panel.background = element_rect(fill = "#343E48", color="#ffffff"), 
            plot.background = element_rect(fill = "#343E48", color = NA))

  })
  
  output$plot3 <- renderPlot({
    plot_tti <- tti[tti$SUBJ_GROUP==input$select,] %>% 
      ungroup() %>% 
      arrange(GENDER,N) %>% 
      mutate(.r=row_number())
    ggplot(data=plot_tti, aes(x=.r, y=N, fill=CATEGORY, tooltip=N)) + geom_bar(stat="identity") + 
      theme_minimal() + scale_fill_brewer(palette="Dark2") + 
      theme(axis.text.x = element_text(angle=45)) + facet_wrap(.~GENDER, scales = "free_y") +
      xlab("Code (Mapped to ICD-9)") + ylab("Number of Patients") +
      scale_x_continuous(breaks = plot_tti$.r, labels=plot_tti$MAPPED_CODE) + coord_flip() +
      theme(panel.background = element_rect(fill = "#343E48", color="#ffffff"), 
            plot.background = element_rect(fill = "#343E48", color = NA),
            axis.text = element_text(color="#FFFFFF"),
            axis.title = element_text(color="#ffffff"),
            legend.text = element_text(color="#ffffff"),
            legend.title = element_text(color="#ffffff"),
            strip.text = element_text(colour = 'white'),
            panel.grid = element_line(color="#303030", size = 0.3))
    #girafe(ggobj=p,   options = list(
    # opts_sizing(rescale = FALSE) ))
  })
  
  output$plot4 <- renderPlot({
    plot_rec <- rec[rec$SUBJ_GROUP==input$select,]
    levs <- unique(plot_rec[order(plot_rec$N), 1])
    ggplot(data=plot_rec, aes(x=DEPARTMENT, y=N, fill=DEPARTMENT)) +
      geom_bar(stat="identity") + theme_minimal() + 
      scale_fill_brewer(palette="Dark2") + coord_flip() +
      theme(panel.background = element_rect(fill = "#343E48",color="#ffffff"), 
            plot.background = element_rect(fill = "#343E48", color = NA,),
            axis.text = element_text(color="#FFFFFF"),
            axis.title = element_text(color="#ffffff"),
            legend.text = element_text(color="#ffffff"),
            legend.title = element_text(color="#ffffff"),
            panel.grid = element_line(color="#303030", size = 0.3))
  })
  
  output$plot5 <- renderPlot({
    if(input$select == "GENOTYPED"){
      plot_lab <- lab[lab$SUBJ_GROUP==input$select & lab$DATASET==input$select_lab,]
      
    } else {
      plot_lab <- lab[lab$DATASET==input$select_lab,]
    }
    
    nobs <- nrow(plot_lab[!is.na(plot_lab$RESULT_VALUE_NUM),])
    nppl <- "N"
    m <- median(plot_lab$RESULT_VALUE_NUM, na.rm = TRUE)
    
    ggplot(data=plot_lab) + 
      geom_boxplot(aes(y=RESULT_VALUE_NUM, x=1.25), width=0.25, color="white", fill=NA) +
      geom_jitter(aes(y=RESULT_VALUE_NUM, x=1), color="#1B9E77", position=position_jitter(width = 0.10, height=0), size=0.15, alpha=0.3) + 
      xlim(0.9,1.4) + coord_flip() + theme_minimal() + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), panel.background = element_blank()) + 
      ylab("Result Value") +
      labs(subtitle=paste0("N Obs. = ", nobs, ", N Patients = ", nppl, ", Median = ", signif(m, digits=3))) + 
      theme(panel.background = element_rect(fill = "#343E48",color="#ffffff"), 
            plot.background = element_rect(fill = "#343E48", color = NA,),
            axis.text = element_text(color="#FFFFFF"),
            axis.title = element_text(color="#ffffff"),
            plot.subtitle = element_text(color="#ffffff"),
            legend.text = element_text(color="#ffffff"),
            legend.title = element_text(color="#ffffff"),
            panel.grid = element_line(color="#303030", size = 0.3))
    
  })
  
  output$table <- renderDataTable({datatable(cat[cat$SUBJ_GROUP==input$select, 1:2])})
}

shinyApp(ui, server)
