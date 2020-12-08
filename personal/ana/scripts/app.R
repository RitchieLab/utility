#Author: Anastasia Lucas 
#Date: 2020-12-08
#Code: https://github.com/RitchieLab/utility/blob/master/personal/ana/scripts/app.R

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(ggforce)
library(ggrepel)
library(grid)
library(scales)

###Theme
source("dashboardtheme.R")

###Data
# demo <- read.delim("PMBB_DEMO_mapped_Dec-2020.txt", stringsAsFactors = FALSE)
demo <- as.data.frame(vroom::vroom("PMBB_DEMO_mapped_Dec-2020.txt"))
demo$PMBB_ID <- as.character(demo$PMBB_ID)
#icd <- as.data.frame(data.table::fread("PMBB_ICD9_map_with_category_with_desc_Dec-2020.txt", stringsAsFactors = FALSE, quote=""))
icd <- as.data.frame(vroom::vroom("PMBB_ICD9_map_with_category_with_desc_Dec-2020.txt"))
#icd$PMBB_ID <- as.character(icd$PMBB_ID)

#nicd <- read.delim("PMBB_ICD9_N_with_rollup_with_desc_Dec-2020.txt", stringsAsFactors = FALSE)
nicd <- as.data.frame(vroom::vroom("PMBB_ICD9_N_with_rollup_with_desc_Dec-2020.txt"))
#lab <- as.data.frame(data.table::fread("LABS_STD_SUMMARY_SUBJ_Dec-2020.txt"))
lab <- as.data.frame(vroom::vroom("LABS_STD_SUMMARY_SUBJ_Dec-2020.txt"))
#lab$PMBB_ID <- as.character(lab$PMBB_ID)
lab$RESULT_VALUE <- as.numeric(lab$RESULT_VALUE)
l <- sort(unique(lab$Lab))
m <- c("MEDIAN", "MAX", "MIN")
names(l) <- sort(unique(lab$Lab))
#rec <- read.delim("pmbb_recruitment_location_mapped_Dec-2020.txt")
rec <- as.data.frame(vroom::vroom("pmbb_recruitment_location_mapped_Dec-2020.txt"))
#rec$PMBB_ID <- as.character(rec$PMBB_ID)
cl <- c("All", sort(unique(nicd$Desc[nicd$n>5 & nicd$Group=="GenotypeOrExome" & nicd$Rollup=="Exact Match"])))
#names(i) <- c("All", sort(unique(icd$Desc[icd$n>5 & icd$Group=="Genotyped" & icd$Rollup=="Exact Match"])))
i <- unlist(lapply(strsplit(cl, " "), FUN=`[[`, 1))
names(i) <- cl

###UI
ui <- dashboardPage(
  dashboardHeader(title="PennMed BioBank"),
  dashboardSidebar(disable = TRUE),
  #dashboardBody(tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "style.css")),
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
                   title = "Cohort Selection",
                   solidHeader = TRUE,
                   status = "primary",
                   selectInput("select", "Select entire Penn Medicine BioBank or patients with genetic data", 
                               c("PMBB" = "PMBB", 
                                 "Genotype" = "Genotype",
                                 "Exome" = "Exome",
                                 "All genotype & exome" = "Genotype|Exome")),
                   width=NULL
                 )),
          
          column(width=4, box(
            title = "Filter by ICD-9 Code",
            solidHeader = TRUE,
            status = "primary",
            selectInput("selecticd", "ICD codes were mapped to ICD-9 using CMS General Equivalence Mappings", i),
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
                   title = "Clinical Labs",
                   solidHeader = TRUE,
                   status = "primary",
                   selectInput("select_lab", "Currently data on 34 labs are available", l),
                   width=NULL
                 )
          ),
          column(width=4,
                 box(
                   title = "Clinical Labs",
                   solidHeader = TRUE,
                   status = "primary",
                   selectInput("select_metric", "Metric (per patient aggregates)", m),
                   width=NULL
                 )
                 
          )
        )
      ),
      fluidRow(
        column(width = 12,
               box(
                 title = "FAQ", width = NULL, background = "navy",status = "primary", solidHeader = TRUE,
                 "In the case that a patient in Penn Medicine BioBank has multiple tissue specimens, they will be included in counts for both 'Blood' and 'Tissue.'\n
                 All ICD codes were mapped to ICD-9 using the General Equivalency Mappings (GEMs) developed by the Centers for Medicare and Medicaid Services and CDC's National Center for Health Statistics.\n
                 ICD multimorbidities were calculated by the 10 codes with highest number of individuals given the specified cohort and ICD filters.\n
                 Currently it is possible to filter by all ICD-9 codes with >= 5 patients in the smallest cohort (i.e. genotyped individuals) in order to generate a distribution.\n 
                 Clinical lab values were converted into standard units when possible (ex. mg/dL can be converted to g/dL mathematically) and excluded when no such conversion could be made.\n
                 Individual clinical lab values were not otherwise removed.\n
                 ***This dashboard is includes genotyped samples included in the 2019 genotype data release and the September 2020 release of exome data.***"
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
    if(input$select!="PMBB"){
      if(input$selecticd=="All"){
        plot_age <- demo[which(grepl(input$select,demo$SUBJ_GROUP)),]
      } else { 
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        plot_age <- demo[which(grepl(input$select, demo$SUBJ_GROUP) & demo$PMBB_ID %in% ids ),] 
      }
    #} else if(input$select=="GenotypeOrExome") {
    #  if(input$selecticd=="All"){
    #    plot_age <- demo[which(grepl("Genotype|Exome",demo$SUBJ_GROUP)),]
    #  } else { 
    #    ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
    #    plot_age <- demo[which(grepl(input$select, demo$SUBJ_GROUP) & demo$PMBB_ID %in% ids ),] 
    #  }
    } else {
      if(input$selecticd=="All"){
        plot_age <- demo
      } else {
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        plot_age <- demo[which(demo$PMBB_ID %in% ids ),] 
      }
    }
    
    acols <- c("#458CFD", "#7A0403")
    names(acols) <- c("F", "M")
    ggplot(data=plot_age[plot_age$AGE>0,], aes(x=AGE, color=GENDER_CODE)) + geom_density() + theme_minimal() + 
      #scale_color_brewer(name='Gender', palette = "Dark2") + 
      scale_color_manual(values=acols) +
      xlab("Patient Age (Years)") +
      theme(
        panel.background = element_rect(fill = "#ffffff",color="#f2f2f3"), 
        plot.background = element_rect(fill = "#ffffff", color = NA),
        axis.text = element_text(color="#000000"),
        axis.title = element_text(color="#000000"),
        legend.text = element_text(color="#000000"),
        legend.title = element_text(color="#000000"),
        text=element_text(size=15),
        panel.grid = element_line(color="#cfd0d2", size = 0.3)
      )
    
  })
  
  ###RACE  
  output$plot2 <- renderPlot({
    if(input$select!="PMBB"){
      if(input$selecticd=="All"){
        eth <- demo[which(grepl(input$select, demo$SUBJ_GROUP)),] 
      } else { 
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        eth <- demo[which(grepl(input$select, demo$SUBJ_GROUP) & demo$PMBB_ID %in% ids ),] 
      }
    } else {
      if(input$selecticd=="All"){
        eth <- demo
      } else {
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        eth <- demo[which(demo$PMBB_ID %in% ids ),] 
      }
    }
    plot_eth <- eth %>% 
      group_by(RACE) %>% tally() %>% mutate(PCT=n/sum(n))
    
    rlevs=rev(c("White", "Black", "Other", "Asian", "N/A"))
    #rcols <- c("#E7298A","#D95F02","#1B9E77","#7570B3")
    rcols <- c("#30123B", "#28BBEB", "#A2FC3C", "#FB8021", "#7A0403")
    names(rcols) <- rlevs
    
    plot_eth$Label <- paste(signif(plot_eth$PCT*100, digits = 2), "%")
    ggplot(plot_eth) + geom_bar(aes(x=factor(RACE, levels=rev(rlevs)), y=n, fill = factor(RACE, levels=rev(rlevs))), stat="identity") +
      scale_fill_manual(values =  rcols)  +
      theme_minimal() +
      theme(panel.background = element_rect(fill = "#ffffff", color="#ffffff"), 
            plot.background = element_rect(fill = "#ffffff", color = "#ffffff"),
            axis.text.x = element_text(angle=30),
            text=element_text(size=15),
            axis.text = element_text(color="#000000"),
            axis.title = element_text(color="#000000")) + 
      guides(fill="none") +
      labs(x="Race", y="Number of Patients") +
      geom_text(aes(label=Label, x=factor(RACE, levels=rev(rlevs)), y=n+(max(n)/20)), position=position_dodge(0.9), vjust=0) +
      scale_y_continuous(labels = comma, breaks = pretty_breaks())
  })
  
  ###TABLE
  output$table <- renderTable(
    if(input$selecticd=="All"){
      if(input$select!="PMBB"){
        data.frame(Codes=c(formatC(length(unique(icd$GEM_ICD9[grepl(input$select, icd$SUBJ_GROUP)])), big.mark = ",")),
                   Patients=c(formatC(length(unique(icd$PMBB_ID[grepl(input$select, icd$SUBJ_GROUP)])), big.mark=",")),
                   Group=c(as.character(input$select)))
      } else {
        data.frame(Codes=c(formatC(length(unique(icd$GEM_ICD9)), big.mark = ",")),
                   Patients=c(formatC(length(unique(icd$PMBB_ID)), big.mark=",")),
                   Group=c("PMBB")
        )
      }
    } else { 
      if(nchar(as.character(input$selecticd))==6){
        tab_nicd <- rbind( nicd[which(nicd$Code==input$selecticd & nicd$Rollup=="Exact Match" & grepl(input$select, nicd$Group)) , ],
                           nicd[which(nicd$Code==substr(input$selecticd, 1, nchar(input$selecticd)-1) & nicd$Rollup=="4-digit" & grepl(input$select, nicd$Group)),],
                           nicd[which(nicd$Code==substr(input$selecticd, 1, nchar(input$selecticd)-3) & nicd$Rollup=="3-digit" & grepl(input$select, nicd$Group)),]
        )
      } else {
        if(nchar(as.character(input$selecticd))==5){
          tab_nicd <- rbind(nicd[which(nicd$Code==input$selecticd & nicd$Rollup=="Exact Match" & grepl(input$select, nicd$Group)), ],
                            nicd[which(nicd$Code==input$selecticd & nicd$Rollup=="4-digit" & grepl(input$select, nicd$Group)), ],
                            nicd[which(nicd$Code==substr(input$selecticd, 1, nchar(input$selecticd)-2) & nicd$Rollup=="3-digit" & grepl(input$select, nicd$Group)), ]
                            
          )
        } else {
          tab_nicd <- rbind(nicd[which(nicd$Code==input$selecticd & nicd$Rollup=="Exact Match" & grepl(input$select, nicd$Group)), ],
                            nicd[which(nicd$Code==input$selecticd & nicd$Rollup=="3-digit" & grepl(input$select, nicd$Group)), ]
          )
          
        }
      }
      tab_nicd$Patients <- formatC(tab_nicd$n, big.mark=",")
      tab_nicd$Description <- sub(".*? ", "", tab_nicd$Desc)
      tab_nicd[, c(1,2,6,7)]
    } 
  )
  
  ###ICD COMORBIDITY  
  output$plot3 <- renderPlot({
    ##Order factor labels of facet
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
    if(input$select!="PMBB"){
      if(input$selecticd=="All"){
        if(input$select=="Genotype") {
          # fill table using following code
          # unique(icd[grepl("Genotype", icd$SUBJ_GROUP), -2]) %>% group_by(GEM_ICD9, Desc, Category) %>% tally() %>% arrange(desc(n))
          plot_icd <-structure(list(GEM_ICD9 = c(), 
                                    Category = c(), 
                                   n = c()), 
                              class = "data.frame", row.names = c(NA, -10L))
        } else if(input$select=="Exome"){ 
          # fill table using following code
          # unique(icd[grepl("Exome", icd$SUBJ_GROUP), -2]) %>% group_by(GEM_ICD9, Desc, Category) %>% tally() %>% arrange(desc(n))
          plot_icd <-structure(list(GEM_ICD9 = c(), 
                                    Category = c(), 
                                    n = c()), 
                               class = "data.frame", row.names = c(NA, -10L))          
        } else { # Genotype | Exome
          # fill table using following code
          # unique(icd[grepl("Genotype|Exome", icd$SUBJ_GROUP), -2]) %>% group_by(GEM_ICD9, Desc, Category) %>% tally() %>% arrange(desc(n))
          plot_icd <-structure(list(GEM_ICD9 = c(), 
                                    Category = c(), 
                                    n = c()), 
                               class = "data.frame", row.names = c(NA, -10L))
        }
      } else { 
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        pre_icd <- icd[which(grepl(input$select, icd$SUBJ_GROUP) & icd$PMBB_ID %in% ids ),] 
        plot_icd <- head(unique(pre_icd[,c("GEM_ICD9", "Category", "PMBB_ID")]) %>% group_by(GEM_ICD9, Category) %>% tally() %>% arrange(desc(n)), n=10)
      }
    } else {
      if(input$selecticd=="All"){
        # fill table using following code
        # unique(icd[, -2]) %>% group_by(GEM_ICD9, Desc, Category) %>% tally() %>% arrange(desc(n))
        plot_icd <- structure(list(GEM_ICD9 = c(), 
                                   Category = c(), 
                                   n = c()), 
                              class = "data.frame", row.names = c(NA, -10L))
        
      } else {
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        pre_icd <- icd[which(icd$PMBB_ID %in% ids ),] 
        plot_icd <- head(unique(pre_icd[,c("GEM_ICD9", "Category", "PMBB_ID")]) %>% group_by(GEM_ICD9, Category) %>% tally() %>% arrange(desc(n)), n=10)
      }
    }
    #dput(plot_icd)
    #plot_icd <- head(unique(pre_icd[,1:3]) %>% group_by(GEM_ICD9, Category) %>% tally() %>% arrange(desc(n)), n=10)
    ilevs <- plot_icd$GEM_ICD9
    
    #icols <- c("#30123B", "#3F3E9C", "#4666DD", "#458CFD", "#2FB2F4", "#1AD4D0",
    #          "#22EBAA", "#52FA7A", "#8BFF4B", "#B7F735", "#DBE236", "#F5C53A",
    #          "#FEA130", "#F9751D", "#E84B0C", "#CE2D04", "#A91601", "#7A0403")
    icols <- c("#CE2D04", "#458CFD", "#52FA7A", "#B7F735", "#F5C53A", "#A91601",
               "#8BFF4B", "#1AD4D0", "#DBE236", "#E84B0C", "#FEA130", "#F9751D",
               "#4666DD", "#7A0403", "#2FB2F4", "#30123B", "#3F3E9C", "#22EBAA")
    names(icols) <- c("Circulatory","Endocrine/metabolic","Morbidity/mortality","Digestive","Symptoms/findings","Musculoskeletal",
                      "Blood/immune","Nervous","Genitourinary","Respiratory","Neoplasms","Skin",
                      "Injury/external","Infectious/parasitic", "Other","Congenital/chromosomal","Mental/behavioural","Perinatal")
    
    ggplot(data=plot_icd, aes(x=factor(GEM_ICD9, levels = rev(ilevs)), y=n, fill=Category)) + geom_bar(stat="identity") + 
      theme_minimal() + scale_fill_manual(values = icols) + 
      theme(axis.text.x = element_text(angle=30)) + #facet_wrap(.~GENDER, scales = "free_y") +
      xlab("") + ylab("Number of Patients") +
      #scale_x_continuous(breaks = plot_icd$.r, labels=plot_icd$GEM_ICD9) + 
      coord_flip() +
      theme(panel.background = element_rect(fill = "#ffffff", color="#f2f2f3"), 
            plot.background = element_rect(fill = "#ffffff", color = NA),
            axis.text = element_text(color="#000000"),
            axis.title = element_text(color="#000000"),
            legend.text = element_text(color="#000000"),
            legend.title = element_text(color="#000000"),
            strip.text = element_text(colour = '#000000'),
            text=element_text(size=15),
            panel.grid = element_line(color="#cfd0d2", size = 0.3)) +
      scale_y_continuous(labels = comma)
  })
  
  ###RECRUITMENT
  output$plot4 <- renderPlot({
    if(input$select!="PMBB"){
      if(input$selecticd=="All"){
        plot_rec <- rec[which(grepl(input$select, rec$SUBJ_GROUP)),]
      } else { 
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        plot_rec <- rec[which(grepl(input$select, rec$SUBJ_GROUP) & rec$PMBB_ID %in% ids ),] 
      }
    } else {
      if(input$selecticd=="All"){
        plot_rec <- rec
      } else {
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        plot_rec <- rec[which(rec$PMBB_ID %in% ids ),] 
      }
    }
    cols <- c("#458CFD", "#A2FC3C", "#7A0403")
    names(cols) <- unique(rec$SPECIMEN)
    #levs <- unique(plot_rec[order(plot_rec$N), 1])
    p <- ggplot(data=plot_rec, aes(x=forcats::fct_infreq(CATEGORY), fill=SPECIMEN)) +
      geom_bar() + theme_minimal() + 
      scale_fill_manual(values=cols) + coord_flip() +
      theme(panel.background = element_rect(fill = "#ffffff",color="#f2f2f3"), 
            plot.background = element_rect(fill = "#ffffff", color = NA),
            axis.text = element_text(color="#000000"),
            axis.text.x = element_text(angle=30),
            axis.title = element_text(color="#000000"),
            legend.text = element_text(color="#000000"),
            legend.title = element_text(color="#000000"),
            text=element_text(size=15),
            panel.grid = element_line(color="#cfd0d2", size = 0.3),
            axis.title.y = element_blank()) +
      guides(fill="none") +
      scale_y_continuous(labels = comma) +
      labs(y="Number of Patients")
    #if(length(unique(plot_rec$SPECIMEN))>1){ p <- p + facet_wrap(.~SPECIMEN, scales="free_y") }
    p
  })
  
  ###LABS  
  output$plot5 <- renderPlot({
    if(input$select!="PMBB"){
      if(input$selecticd=="All"){
        plot_lab <- lab[which(grepl(input$select, lab$SUBJ_GROUP) & lab$Lab==input$select_lab & lab$METRIC==input$select_metric),]
      } else { 
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        plot_lab <- lab[which(grepl(input$select, lab$SUBJ_GROUP) & lab$Lab==input$select_lab & lab$PMBB_ID %in% ids & lab$METRIC==input$select_metric),] 
      }
    } else {
      if(input$selecticd=="All"){
        plot_lab <- lab[lab$Lab==input$select_lab & lab$METRIC==input$select_metric,]
      } else {
        ids <- unique(icd$PMBB_ID[icd$GEM_ICD9==input$selecticd])
        plot_lab <- lab[which(lab$PMBB_ID %in% ids & lab$Lab==input$select_lab & lab$METRIC==input$select_metric),] 
      }
    }
    
    nobs <- length(unique(plot_lab$PMBB_ID[!is.na(plot_lab$RESULT_VALUE)]))
    
    un <- unique(plot_lab$UNIT)
    
    ggplot(data=plot_lab) + 
      geom_boxplot(aes(y=RESULT_VALUE, x=1.25), width=0.25, color="#00144d", fill=NA) +
      geom_jitter(aes(y=RESULT_VALUE, x=1), color="#00144d", position=position_jitter(width = 0.10, height=0), alpha=0.3) + 
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
            text=element_text(size=15),
            panel.grid = element_line(color="#cfd0d2", size = 0.3))
    
  })
  
}

shinyApp(ui, server)

