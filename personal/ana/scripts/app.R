library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)
library(data.table)
library(dplyr)
library(ggiraph)
library(ggforce)
library(ggrepel)
#library(dashboardthemes)
library(grid)
library(scales)

source("dashboardtheme.R")
demo <- read.delim("PMBB_DEMO.txt", stringsAsFactors = FALSE)
tti <- read.delim("Top_10_ICD.txt", stringsAsFactors = FALSE)
cat <- read.delim("PMBB_ICD_CAT_N.txt", stringsAsFactors = FALSE)
cat$Category[cat$Category=="NULL"] <- "Other"
lab <- as.data.frame(data.table::fread("LABS_STD_SUMMARY_SUBJ.txt"))
lab$RESULT_VALUE <- as.numeric(lab$RESULT_VALUE)
l <- unique(lab$Lab)
m <- c("MEDIAN", "MAX", "MIN")
names(l) <- unique(lab$Lab)
rec <- read.delim("pmbb_recruitment_location.txt")

ui <- dashboardPage(
  dashboardHeader(title="PMBB Demographics"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    shinyDashboardThemes(
      theme = "pmbb"
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
                 title = "Subject Selection",
                 solidHeader = TRUE,
                 status = "primary",
                 selectInput("select", "Group:", c("All PMBB" = "PMBB", "Genotyped" = "GENOTYPED")),
                 width=NULL
               ),
               box(
                 dataTableOutput("table", height = 250), 
                 title="Number of Patients per ICD Category",
                 solidHeader = FALSE,
                 status = "info",
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
        fluidRow(
          column(width=4,
                 box(
                   title = "Clinical Lab Selection",
                   solidHeader = TRUE,
                   status = "primary",
                   selectInput("select_lab", "Lab:", l),
                   width=NULL
                 )
                 ),
        column(width=4,
               box(
                 title = "Clinical Lab Selection",
                 solidHeader = TRUE,
                 status = "primary",
                 selectInput("select_metric", "Metric:", m),
                 width=NULL
               )
               
        )
        )
      )
    )
  )
)

server <- function(input, output) { 
  
  
  output$plot1 <- renderPlot({
    if(input$select=="GENOTYPED"){
      plot_age <- demo[demo$SUBJ_GROUP==input$select,]
    } else {
      plot_age <- demo
    }
    
    ggplot(data=plot_age[plot_age$AGE_20200101>0,], aes(x=AGE_20200101, color=GENDER_CODE)) + geom_density() + theme_minimal() + 
      scale_color_brewer(name='Gender', palette = "Dark2") + xlab("Patient Age (Years)") +
      theme(
            panel.background = element_rect(fill = "#ffffff",color="#f2f2f3"), 
            plot.background = element_rect(fill = "#ffffff", color = NA),
            axis.text = element_text(color="#000000"),
            axis.title = element_text(color="#000000"),
            legend.text = element_text(color="#000000"),
            legend.title = element_text(color="#000000"),
            panel.grid = element_line(color="#cfd0d2", size = 0.3)
      )
    
  })
  
  output$plot2 <- renderPlot({
    if(input$select=="GENOTYPED"){
      plot_eth <- demo[demo$SUBJ_GROUP==input$select,] %>% 
          group_by(RACE) %>% tally() %>% mutate(PCT=n/sum(n))
    } else {
      plot_eth <- demo %>% group_by(RACE) %>% tally() %>% mutate(PCT=n/sum(n))
    }

    levs=rev(c("WHITE", "BLACK", "OTHER", "ASIAN"))
  
      plot_eth$Label <- paste(signif(plot_eth$PCT*100, digits = 2), "%")
      ggplot(plot_eth) + geom_bar(aes(x=factor(RACE, levels=rev(levs)), y=n, fill = factor(RACE, levels=rev(levs))), stat="identity") +
      scale_fill_brewer(palette = "Dark2")  +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "#ffffff", color="#ffffff"), 
            plot.background = element_rect(fill = "#ffffff", color = "#ffffff"),
            axis.text.x = element_text(angle=45)) + 
      guides(fill="none") +
      labs(x="RACE", y="Number of Patients") +
      geom_text(aes(label=Label, x=factor(RACE, levels=rev(levs)), y=n+500), position=position_dodge(0.9), vjust=0) +
      scale_y_continuous(labels = comma)
  })

#  output$table <- renderDataTable({
#    plot_cat <- cat[cat$SUBJ_GROUP==input$select, 1:2]
#    names(plot_cat) <- paste0('<span style="color:',c("white","white"),'">',colnames(plot_cat)[c(1,2)],'</span>')
#    datatable(plot_cat, escape=FALSE) %>%
#      formatStyle(columns=1, color="white") %>%
#      formatStyle(columns=2, color="white")
#  })  
 
  output$table <- renderDataTable(
    datatable(cat[cat$SUBJ_GROUP==input$select, 1:2]),
    selection = list(mode = 'single', server=TRUE)
    )

  output$plot3 <- renderPlot({
    if (length(input$table_rows_selected)){
      print(paste0("this",as.character(input$table_rows_selected)))
      j <- as.integer(as.character(input$table_rows_selected))
      print(j)
      i <- cat[cat$SUBJ_GROUP==input$select, 1][j]
      print(i)
      plot_tti <- tti[tti$SUBJ_GROUP==input$select & tti$CATEGORY==i,] %>% 
        ungroup() %>% 
        arrange(GENDER,N) %>% 
        mutate(.r=row_number())
    } else {
      
    plot_tti <- tti[tti$SUBJ_GROUP==input$select,] %>% 
      ungroup() %>% 
      arrange(GENDER,N) %>% 
      mutate(.r=row_number())
    }
    ggplot(data=plot_tti, aes(x=.r, y=N, fill=CATEGORY, tooltip=N)) + geom_bar(stat="identity") + 
      theme_minimal() + scale_fill_brewer(palette="Dark2") + 
      theme(axis.text.x = element_text(angle=45)) + facet_wrap(.~GENDER, scales = "free_y") +
      xlab("Code (Mapped to ICD-9)") + ylab("Number of Patients") +
      scale_x_continuous(breaks = plot_tti$.r, labels=plot_tti$MAPPED_CODE) + coord_flip() +
      theme(panel.background = element_rect(fill = "#ffffff", color="#f2f2f3"), 
            plot.background = element_rect(fill = "#ffffff", color = NA),
            axis.text = element_text(color="#000000"),
            axis.title = element_text(color="#000000"),
            legend.text = element_text(color="#000000"),
            legend.title = element_text(color="#000000"),
            strip.text = element_text(colour = '#000000'),
            panel.grid = element_line(color="#cfd0d2", size = 0.3))
    #girafe(ggobj=p,   options = list(
    # opts_sizing(rescale = FALSE) ))
  })
  
  output$plot4 <- renderPlot({
    if(input$select=="GENOTYPED"){
      plot_rec <- rec[which(rec$SUBJ_GROUP==input$select),]
    } else {
      plot_rec <- rec
    }
    cols <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E",
              "#E6AB02","#A6761D","#666666","#1B9E77","#D95F02","#7570B3")
    names(cols) <- unique(rec$LOCATION)
    #levs <- unique(plot_rec[order(plot_rec$N), 1])
    ggplot(data=plot_rec, aes(x=LOCATION, fill=LOCATION)) +
      geom_bar() + theme_minimal() + 
      scale_fill_manual(values=cols) + coord_flip() +
      theme(panel.background = element_rect(fill = "#ffffff",color="#f2f2f3"), 
            plot.background = element_rect(fill = "#ffffff", color = NA),
            axis.text = element_text(color="#000000"),
            axis.title = element_text(color="#000000"),
            legend.text = element_text(color="#000000"),
            legend.title = element_text(color="#000000"),
            panel.grid = element_line(color="#cfd0d2", size = 0.3),
            axis.title.y = element_blank()) +
      guides(fill="none") +
      facet_wrap(.~SPECIMEN, scales="free") +
      scale_y_continuous(labels = comma) +
      labs(y="Number of Patients")
  })
  
  output$plot5 <- renderPlot({
    if(input$select == "GENOTYPED"){
      plot_lab <- lab[lab$SUBJ_GROUP==input$select & lab$Lab==input$select_lab & lab$METRIC==input$select_metric,]
      
    } else {
      plot_lab <- lab[lab$Lab==input$select_lab & lab$METRIC==input$select_metric,]
    }
    
    nobs <- nrow(plot_lab[!is.na(plot_lab$RESULT_VALUE),])

    un <- unique(plot_lab$UNIT)
    
    ggplot(data=plot_lab) + 
      geom_boxplot(aes(y=RESULT_VALUE, x=1.25), width=0.25, color="#1B9E77", fill=NA) +
      geom_jitter(aes(y=RESULT_VALUE, x=1), color="#1B9E77", position=position_jitter(width = 0.10, height=0), size=0.5, alpha=0.3) + 
      xlim(0.9,1.4) + coord_flip() + theme_minimal() + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), panel.background = element_blank()) + 
      ylab(paste0("Result Value (", un, ")")) +
      labs(subtitle=paste0("Number of Patients = ", formatC(nobs, big.mark=","))) +
      theme(panel.background = element_rect(fill = "#ffffff",color="#f2f2f3"), 
            plot.background = element_rect(fill = "#ffffff", color = NA),
            axis.text = element_text(color="#000000"),
            axis.title = element_text(color="#000000"),
            plot.subtitle = element_text(color="#000000"),
            legend.text = element_text(color="#000000"),
            legend.title = element_text(color="#000000"),
            panel.grid = element_line(color="#cfd0d2", size = 0.3))
    
  })
  
}

shinyApp(ui, server)
