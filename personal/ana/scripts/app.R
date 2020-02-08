library(shiny)
library(shinydashboard)
library(ggplot2)
library(DT)
library(waffle)
library(data.table)

ui <- dashboardPage(
  dashboardHeader(title="PMBB Demographics"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    # Boxes need to be put in a row (or column)
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
      ),
      fluidRow(
        column(width=12,
          box(
            plotOutput("plot5", height = 1000),
            title="Clinical Labs",
            solidHeader=FALSE,
            status="info",
            width=NULL
          )
        )
      )
    )
  )
)

server <- function(input, output) { 
  age = read.delim("PMBB_AGE_SEX.txt", stringsAsFactors = FALSE)
  age$AGE <- as.numeric(age$AGE)
  eth <- read.delim("PMBB_RACE_N.txt", stringsAsFactors = FALSE)
  eth$PCT <- ifelse(eth$SUBJ_GROUP=="PMBB", eth$N/53611, eth$N/18981)
  tti <- read.delim("Top_10_ICD.txt", stringsAsFactors = FALSE)
  cat <- read.delim("PMBB_ICD_CAT_N.txt", stringsAsFactors = FALSE)
  cat$Category[cat$Category=="NULL"] <- "Other"
  lab <- as.data.frame(data.table::fread("PMBB_LABS.txt"))
  lab$RESULT_VALUE_NUM <- as.numeric(lab$RESULT_VALUE_NUM)
  rec <- read.delim("Department_Recruitment_dummy.txt")

  output$plot1 <- renderPlot({
    if(input$select == "GENOTYPED"){
      plot_age <- age[age$SUBJ_GROUP==input$select,]
    } else {
      plot_age <- age
    }
    ggplot(data=plot_age, aes(x=AGE, color=GENDER_CODE)) + geom_density() + theme_minimal() + 
      scale_color_brewer(name='Gender', palette = "Dark2") + xlab("Patient Age (Years)")
  })
  
  output$plot2 <- renderPlot({
    plot_eth <- eth[eth$SUBJ_GROUP==input$select,]
    if(input$select=="PMBB"){
      levs=rev(c("WHITE", "BLACK", "UNKNOWN", "OTHER", "ASIAN", "AM IND AK NATIVE", "HI PAC ISLAND"))
    } else {
      levs=rev(c("WHITE", "BLACK", "UNKNOWN", "OTHER", "ASIAN", "HI PAC ISLAND","AM IND AK NATIVE"))
    }
    #ggplot(data=plot_eth, aes(x=factor(RACE_CODE, levels=levs), y=N, fill=RACE_CODE)) +
    #  geom_bar(stat="identity") + theme_minimal() + scale_fill_brewer(palette="Dark2") + coord_flip() +
    #  xlab("") + ylab("Number of Patients") + theme(axis.text.x = element_text(angle=45)) + guides(fill=FALSE)
    ggplot(eth[eth$SUBJ_GROUP=="PMBB",], aes(fill=factor(RACE_CODE, levels=levs), values=N)) +  
      geom_waffle(n_rows=10, size=0.33, colour="white", flip=TRUE, make_proportional = TRUE) + 
      theme_enhance_waffle() + scale_fill_brewer(name="RACE", palette="Dark2") + coord_equal() + 
      theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
  })
  
  output$plot3 <- renderPlot({
    plot_tti <- tti[tti$SUBJ_GROUP==input$select,]
    ggplot(data=plot_tti, aes(x=MAPPED_CODE, y=N, fill=CATEGORY)) + geom_bar(stat="identity") + 
      theme_minimal() + scale_fill_brewer(palette="Dark2") + theme(axis.text.x = element_text(angle=45)) + 
      facet_wrap(.~GENDER, scales = "free_x") + xlab("Code (Mapped to ICD-9)") + ylab("Number of Patients")
  })
  
  output$plot4 <- renderPlot({
    plot_rec <- rec[rec$SUBJ_GROUP==input$select,]
    levs <- unique(plot_rec[order(plot_rec$N), 1])
    ggplot(data=plot_rec, aes(x=DEPARTMENT, y=N, fill=DEPARTMENT)) +
      geom_bar(stat="identity") + theme_minimal() + 
      scale_fill_brewer(palette="Dark2") + coord_flip()
  })
  
  output$plot5 <- renderPlot({
    if(input$select == "GENOTYPED"){
      plot_age <- lab[lab$SUBJ_GROUP==input$select,]
    } else {
      plot_age <- lab
    }
    ggplot(data=lab, aes(x=RESULT_VALUE_NUM)) + 
      geom_density( fill="#1B9E77") + facet_wrap(.~DATASET, scales="free") +
      theme_minimal() + xlab("Result Value")
  })
  
  output$table <- renderDataTable({datatable(cat[cat$SUBJ_GROUP==input$select, 1:2])})
  #https://orchid00.github.io/tidy_raincloudplot
}

shinyApp(ui, server)
