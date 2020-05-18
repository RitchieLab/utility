#Author: Anastasia Lucas 
#Date: 2020-05-17
#Code: https://github.com/RitchieLab/utility/blob/master/personal/ana/scripts/app.R

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(ggiraph)
library(ggforce)
library(ggrepel)
library(grid)
library(scales)

###Theme
source("dashboardtheme.R")

###Data
demo <- read.delim("PMBB_DEMO.txt", stringsAsFactors = FALSE)
demo$PT_ID <- as.character(demo$PT_ID)
icd <- as.data.frame(data.table::fread("PMBB_ICD9_map_with_category.txt", stringsAsFactors = FALSE))
icd$PT_ID <- as.character(icd$PT_ID)
nicd <- read.delim("PMBB_ICD_N.txt", stringsAsFactors = FALSE)
lab <- as.data.frame(data.table::fread("LABS_STD_SUMMARY_SUBJ.txt"))
lab$PT_ID <- as.character(lab$PT_ID)
lab$RESULT_VALUE <- as.numeric(lab$RESULT_VALUE)
l <- unique(lab$Lab)
m <- c("MEDIAN", "MAX", "MIN")
names(l) <- unique(lab$Lab)
rec <- read.delim("pmbb_recruitment_location.txt")
rec$PT_ID <- as.character(rec$PT_ID)
i <- c("All", sort(unique(nicd$Code[nicd$n>5 & nicd$Group=="Genotyped"])))
names(i) <- c("All", sort(unique(nicd$Code[nicd$n>5 & nicd$Group=="Genotyped"])))

###UI
ui <- dashboardPage(
  dashboardHeader(title="PMBB Demographics"),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    shinyDashboardThemes(
      theme = "pmbb"
    ),
    fluidPage(
        fluidRow(
          column(width=4,
               box(
                 plotOutput("plot1", height = 260), 
                 title="Age Distribution",
                 solidHeader = FALSE,
                 status = "info",
                 width=NULL
               )),
               
               column(width=4, box(
                 plotOutput("plot2", height = 260), 
                 title="Self-Reported Ancestry",
                 solidHeader = FALSE,
                 status = "info",
                 width=NULL
               )),
            fluidRow(
                column(width=4,
                 box(
                   title = "Subject Selection",
                   solidHeader = TRUE,
                   status = "primary",
                   selectInput("select", "Group:", c("All PMBB" = "PMBB", "Genotyped" = "Genotyped")),
                   width=NULL
                 )),
                
               column(width=4, box(
                 title = "ICD-9 Code Filter",
                 solidHeader = TRUE,
                 status = "primary",
                 selectInput("selecticd", "Code:", i),
                 width=NULL
               )))),
        
      fluidRow(
        
        column(width=4, box(
          plotOutput("plot3", height = 260), 
          #plotOutput("plot3", height = 260),
          title="Top 10 Multimorbidities",
          solidHeader = FALSE,
          status = "info",
          width=NULL
        )
      ),
        column(width=4,
               
               box(
                 plotOutput("plot4", height=260),
                 title = "Recruitment Location",
                 solidHeader = FALSE,
                 status="info",
                 width=NULL
               )
        ),
        column(width = 4,
               box(
                 tableOutput("table"), 
                 title="Number of Patients per Code",
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

###Server
server <- function(input, output) { 
  
###AGE  
  output$plot1 <- renderPlot({
    if(input$select=="Genotyped"){
      if(input$selecticd=="All"){
        plot_age <- demo[which(demo$SUBJ_GROUP==input$select),]
      } else { 
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        plot_age <- demo[which(demo$SUBJ_GROUP==input$select & demo$PT_ID %in% ids ),] 
      }
    } else {
      if(input$selecticd=="All"){
        plot_age <- demo
      } else {
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        plot_age <- demo[which(demo$PT_ID %in% ids ),] 
      }
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

  ###RACE  
  output$plot2 <- renderPlot({
    if(input$select=="Genotyped"){
      if(input$selecticd=="All"){
        eth <- demo[which(demo$SUBJ_GROUP==input$select),] 
      } else { 
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        eth <- demo[which(demo$SUBJ_GROUP==input$select & demo$PT_ID %in% ids ),] 
      }
    } else {
      if(input$selecticd=="All"){
        eth <- demo
      } else {
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        eth <- demo[which(demo$PT_ID %in% ids ),] 
      }
    }
    plot_eth <- eth %>% 
         group_by(RACE) %>% tally() %>% mutate(PCT=n/sum(n))

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
      geom_text(aes(label=Label, x=factor(RACE, levels=rev(levs)), y=n+(max(plot_eth$n)/20)), position=position_dodge(0.9), vjust=0) +
      scale_y_continuous(labels = comma, breaks = pretty_breaks())
  })
  
###TABLE
  output$table <- renderTable(
    #print(input$selecticd)
     if(input$selecticd=="All"){
        data.frame(Group=c("All PMBB", "Genotyped"),
                   Codes=c(formatC(length(unique(icd$GEM_ICD9)), big.mark = ","), formatC(length(unique(icd$GEM_ICD9[icd$SUBJ_GROUP=="Genotyped"])), big.mark=",")),
                   Patients=c(formatC(length(unique(icd$PT_ID)), big.mark = ","), formatC(length(unique(icd$PT_ID[icd$SUBJ_GROUP=="Genotyped"])), big.mark=",")))
     } else { 
       if(nchar(as.character(input$selecticd))==6){
         tab_nicd <- nicd[which(nicd$Code==input$selecticd | nicd$Code==substr(input$selecticd, 1, nchar(input$selecticd)-1) | nicd$Code==substr(input$selecticd, 1, nchar(input$selecticd)-3)), ]
       } else {
         if(nchar(as.character(input$selecticd))==5){
           tab_nicd <- nicd[which(nicd$Code==input$selecticd | nicd$Code==substr(input$selecticd, 1, nchar(input$selecticd)-2)) , ]
         } else {
             tab_nicd <- nicd[which(nicd$Code==input$selecticd), ]
           }
       }
       tab_nicd$Patients <- formatC(tab_nicd$n, big.mark=",")
       tab_nicd$Rollup <- paste0(nchar(gsub("\\.", "", tab_nicd$Code)), "-digit")
       tab_nicd[, c(1,5,4,3)]
     } 
  )

  ###ICD COMORBIDITY  
  output$plot3 <- renderPlot({
    
    #if (length(input$table_rows_selected)){
    #  print(paste0("this",as.character(input$table_rows_selected)))
    #  j <- as.integer(as.character(input$table_rows_selected))
    #  print(j)
    #  i <- cat[cat$SUBJ_GROUP==input$select, 1][j]
    #  print(i)
    #  plot_icd <- icd[icd$SUBJ_GROUP==input$select & icd$CATEGORY==i,] %>% 
    #    ungroup() %>% 
    #    arrange(GENDER,N) %>% 
    #    mutate(.r=row_number())
    #} else {
    #  
    #  plot_icd <- icd[icd$SUBJ_GROUP==input$select,] %>% 
    #    ungroup() %>% 
    #    arrange(GENDER,N) %>% 
    #    mutate(.r=row_number())
    #}
    if(input$select=="Genotyped"){
      if(input$selecticd=="All"){
        pre_icd <- icd[which(demo$SUBJ_GROUP==input$select),]
      } else { 
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        pre_icd <- icd[which(demo$SUBJ_GROUP==input$select & icd$PT_ID %in% ids ),] 
      }
    } else {
      if(input$selecticd=="All"){
        pre_icd <- icd
      } else {
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        pre_icd <- icd[which(icd$PT_ID %in% ids ),] 
      }
    }
    plot_icd <- head(unique(pre_icd[,1:3]) %>% group_by(GEM_ICD9, Category) %>% tally() %>% arrange(desc(n)), n=10)
    ilevs <- plot_icd$GEM_ICD9
    
    ggplot(data=plot_icd, aes(x=factor(plot_icd$GEM_ICD9, levels = rev(ilevs)), y=n, fill=Category)) + geom_bar(stat="identity") + 
      theme_minimal() + scale_fill_brewer(palette="Dark2") + 
      theme(axis.text.x = element_text(angle=45)) + #facet_wrap(.~GENDER, scales = "free_y") +
      xlab("Code (Mapped to ICD-9)") + ylab("Number of Patients") +
      #scale_x_continuous(breaks = plot_icd$.r, labels=plot_icd$GEM_ICD9) + 
      coord_flip() +
      theme(panel.background = element_rect(fill = "#ffffff", color="#f2f2f3"), 
            plot.background = element_rect(fill = "#ffffff", color = NA),
            axis.text = element_text(color="#000000"),
            axis.title = element_text(color="#000000"),
            legend.text = element_text(color="#000000"),
            legend.title = element_text(color="#000000"),
            strip.text = element_text(colour = '#000000'),
            panel.grid = element_line(color="#cfd0d2", size = 0.3)) +
      scale_y_continuous(labels = comma)
  })
  
  ###RECRUITMENT
  output$plot4 <- renderPlot({
    if(input$select=="Genotyped"){
      if(input$selecticd=="All"){
        plot_rec <- rec[which(rec$SUBJ_GROUP==input$select),]
      } else { 
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        plot_rec <- rec[which(rec$SUBJ_GROUP==input$select & rec$PT_ID %in% ids ),] 
      }
    } else {
      if(input$selecticd=="All"){
        plot_rec <- rec
      } else {
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        plot_rec <- rec[which(rec$PT_ID %in% ids ),] 
      }
    }
    cols <- c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E",
              "#E6AB02","#A6761D","#666666","#1B9E77","#D95F02","#7570B3")
    names(cols) <- unique(rec$LOCATION)
    #levs <- unique(plot_rec[order(plot_rec$N), 1])
    p <- ggplot(data=plot_rec, aes(x=LOCATION, fill=LOCATION)) +
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
      scale_y_continuous(labels = comma) +
      labs(y="Number of Patients")
    if(length(unique(plot_rec$SPECIMEN))>1){ p <- p + facet_wrap(.~SPECIMEN, scales="free") }
    p
  })
  
###LABS  
  output$plot5 <- renderPlot({
    if(input$select=="Genotyped"){
      if(input$selecticd=="All"){
        plot_lab <- lab[which(lab$SUBJ_GROUP==input$select & lab$Lab==input$select_lab & lab$METRIC==input$select_metric),]
      } else { 
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        plot_lab <- lab[which(lab$SUBJ_GROUP==input$select & lab$Lab==input$select_lab & lab$PT_ID %in% ids  & lab$METRIC==input$select_metric),] 
      }
    } else {
      if(input$selecticd=="All"){
        plot_lab <- lab[lab$Lab==input$select_lab & lab$Lab==input$select_lab & lab$METRIC==input$select_metric,]
      } else {
        ids <- unique(icd$PT_ID[icd$GEM_ICD9==input$selecticd])
        plot_lab <- lab[which(lab$PT_ID %in% ids & lab$Lab==input$select_lab & lab$METRIC==input$select_metric),] 
      }
    }
    
    nobs <- length(unique(plot_lab$PT_ID[!is.na(plot_lab$RESULT_VALUE)]))
    
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
