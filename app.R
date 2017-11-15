# Load packages ----

library(shiny)
library(limma)
library(scatterplot3d)
library(xlsx)
library(gplots)
library(seqinr)

options(shiny.maxRequestSize=30*1024^2)

# Demo data loading

targets_demo <- readTargets('data/Targets.txt')

rg_demo <- read.maimages(targets_demo, source = 'bluefuse',path = 'data/')

# User interface ----
ui <- fluidPage(
  titlePanel("Microarray Analisys"),
  
  sidebarLayout(
    sidebarPanel(
      
      helpText('Is your experiment single channel?'),
      
      selectInput('single_ch', 'Single Channel experiment?',
                  choices = list('Yes' = T,
                                 'No' = F)),
      
      helpText("Select a target file. Remember that the target file
               MUST have these columns: FileName Cy3 Cy5"),
      # Input: Select a file ----
      
      fileInput('target', 'Choose the target file'),
      
      # Raw data selection ----
      
      
      helpText('Now upload here all the raw data files'),
      
      fileInput('raw_files', 'Upload Raw files', multiple = T),
      
      # Scanning sofware ---- 
      
      helpText('Select the scanning software used to generate raw files'),
      
      selectInput('source', 'Scanning program',
                  choices = list('Agilent: Feature Extraction' = 'agilent',
                                 'Bluefuse' = 'bluefuse',
                                 'Others' = 'other')),
      # Reference channel ----
      
      
      
      helpText('Which channel is the reference? 
               This information is at the targets file.'),
      
      selectInput('ref', 'Reference Channel', 
                  choices = list('Cy3' = 'Cy3',
                                 'Cy5' = 'Cy5')),
      
      # Statistical significance ----
      
      helpText('Select the Statistical significance threshold'),
      
      selectInput('fdr_pvalue', 'Select FDR
                  or p-value ',
                  choices = list('FDR (Recomended)' = 'adj.P.Val',
                            'P Value (Strongly not recomended)' = 'P.Value')),
      
      sliderInput('alpha', 'Select the significance threshold',
                  min = 0.001, max = 0.5, value = 0.05, step = 0.01),
      
      sliderInput('fc', ' Select the |log2(FC)| threshold',
                  
                  min = 0, max = 5, value = 1, step = 0.25 )
    

    
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        ## Intro panel ----
      
      tabPanel("Introduction",
               
               
                  
              h3("Introduction"),
                  
              p("This Shiny App is an interactive version of the limma
                package. This tool is designed for those users that 
                are not famirialized with R or programming. The main
                objective of this tool is to make microarray analysis
                using the power of limma, without knowing R programming.
                "),
                  
                h4("Target file"),
                  
                p("The first input that has to be defined in this tool
                  is a Target.txt file. This file MUST contain the
                  following header: FileName Cy3 Cy5. The Cy5 column is
                  only required if the experiment is two color channel.
                  In the FileName colum the user should include the 
                  files names that will be included in the analysis.
                  If the analysis is based in a two color chip, one of the 
                  Cy3 or Cy5 column should contain the word Ref in all
                  the rows, because that channel contains the reference 
                  DNA. The other column should have the treatment
                  conditions. You can see an example of the structure 
                  of a target file for a two channel and single channel
                  experiment below."),
              
              h5('Two channel example (Reference channel is Cy3):'),
                  
                tableOutput("target_example"),
              
              h5('Single channel example:'),
              
              tableOutput('target_example_single'),
              
                h4("Raw files and source program"),
              
              p('The next input in this tool are the raw files that are
                specified in the target file. The format of these files
                can be different from one experiment to another, since
                the chips are different. Because of this, the next input
                to the tool is a selection box where you can choose the 
                program that has generated those raw files. By default,
                the tool asumes that Agilent Feature Extraction Software
                has been used.'),
              
              h4('Reference channel and single channel experiments'),
              
              p('If your experiment is based on a two color chip, 
                one of the channels should contain the reference DNA. In
                our example, the Cy3 channel contains the Ref signature in 
                the target file, so Cy3 is the reference channel.'),
              
              p('If the experiment is based on a single channel chip, 
                leave the selector in Cy3 and remember that the 
                target file MUST contain that column.'),
              
              h4('Statistical significance thresholds'),
              
              p('Finally you will find that you can use FDR or raw pvalue
                in order to decide which features are statistically
                different from one condition to the other. By default, FDR
                is selected for testing significance, and we encourage you
                not to use raw pvalues, since it increase type I errors.
                Below of this selection box, you will fin two sliders that
                set the values of the alpha significance threshold and 
                log2(Fold Change) threshold. These values will be used 
                when generating the xlsx file at the end.'),
              
              br()
              
               
               ),
      
      ## Raw data tab ----
      
      tabPanel("Raw data",
      
      h3('Raw Data exploration'),
              
              plotOutput('raw_box'),
              plotOutput('MA_raw'),
              plotOutput('density_raw')
              
          ),
      tabPanel('Normalized data',
               h3('Normalized data'),
               plotOutput('norm_box'),
               plotOutput('MA_norm'),
               plotOutput('density_norm'),
               
               selectInput('pca_shape', 'PCA dimensions', 
                           choices = list('3D' = '3D',
                                          '2D' = '2D')),
               
               plotOutput('pca')
               ),
      
      
      ## Differential expression tab ----
      
      tabPanel('Differential Expression Analysis',
               
               h3('Differential Expression Analysis'),
               p('First of all, select the comparisons that you want to do.
                 The comparissons can be only be done if the group is defined
                 in the Target file. You can also select if you want a single
                 comparison or a multiple comparison. If this is the case,
                 you must write the comparison this way: A-B,C-B. Remember,
                 DONT use spaces. You can see which groups you can
                 compare below:'),
               tableOutput('groups'),
               
               selectInput('mul_comp', 'Is a multiple comparison?', 
                           choices = list('Yes' = T,
                                          'No' = F)),
               
               textInput('contrasts', 'Comparison(s),
                         remember to not use spaces:',
                         value = 'AC-Cntrl'),
               
               p('Now you can see some results according to the FDR and 
                 |log2(FC)| selected. Type below the comparison you want to 
                 see (the 20 most significants genes). For example you can type
                 AC-Cntrl if you are running the demo data'),
               
               textInput('top_contrast', 'Comparison of interest:',
                         value = 'AC-Cntrl'),
               
               p('In this table below you will find the too 20 genes of the
                 comparison. The column which contains the genes ID will be
                 used in further steps.'),
               
               tableOutput('top_table'),
               
               p('Here you can dowload the tables in a single xlsx file. If
                 there are multiple comparisons, they will be located in
                 different sheets of the file.'),
               
               downloadButton("downloadData", "Download the table in xlsx"),
               
               br()
               
               
               ),
      
      ### Fasta annot tab ----
      
      tabPanel('FASTA files',
               p('First of all, you will need to take a look to the 
                 annotation file. You will have to answer some questions
                 about the format and content.'),
               
               p('Which is the character used for delimiting the fields?
                 Usually is a tab but can be also a semicolon (;)'),
               
               selectInput('del', 'Delimitation character',
                           choices = list('Tab' = '\t',
                                          'Semicolon' = ';',
                                          'Comma' = ',')),
               
               
               
               
               
               numericInput('skip', 'How many lines do the app have to
                            skip reading in the annotation?',
                            min = 0, max = Inf, value = 14),
               
               p('Take a look into the table of the top 20 genes
                 in the previous table. Which column contains the ID of the
                 genes?'),
               
               textInput('top_id', 'Insert the name of the column in 
                         which you can found the names of the genes of
                         interest.'),
               
               fileInput('annotation', 'Insert your annotation file'),
               
               tableOutput('columns'),
               
               textInput('col_id','Which column contains the ID name?',
                         value = 'Reporter.Database.Entry..genbank.'),
               
               textInput('col_seq', 'Which column contains the sequences?',
                         value = "Reporter.Sequence"),
               
               p('Now select the comparisson you are interested in. 
                 From this comparisson we will extract the DEGs with
                 the selected statistical significance and log2(FC). 
                 You can see the groups you can compare in the previous tab
                 (Differential Expression).'),
               
               textInput('group_fasta', 'Select the comparisson of interest.
                         Only one.', value = 'AC-Cntrl'),
               
               
               downloadButton('fasta_down', 'Download your DEGs and their 
                              sequence')
               
               ),
      
      ### Advanced results tab -----------
      
      tabPanel('Advanced Results',
               p('The goal of this section is to generate advanced results
                 like volcano plots or heatmaps'))
  )
)
)
)

# Server logic -----

server <- function(input, output) {
  
  
  # Making posible to load other files than the demo.
  # The reactive function allows to use the demo files or the others
  
  rg <- reactive({
    
    print(input$source)
    
    # If there is no data input we will show only a few samples

     if(is.null(input$target)){
  
     return(rg_demo)
  
   } else {
  
     # but if the input is defined, we will use the input files
  
     targets <- readTargets(input$target$datapath)
  
     if(input$single_ch){
       print('one color used')
  
       rg <- read.maimages(input$raw_files$datapath,
                           source = input$source, green.only = T)
       
       rg$targets <- targets
       
       colnames(rg$E) <- rownames(targets)
  
       return(rg)
  
     }else{
  
       rg <- read.maimages(input$raw_files$datapath,
                           source = input$source)
     }
  
  
     rg$targets <- targets
  
     colnames(rg$G) <- rownames(targets)
     colnames(rg$R) <- rownames(targets)
  
     return(rg)
   }
  })
  
  # definition of normalized data ------
  
  norm <- reactive({
    
    if(input$single_ch){ # If the experiment is single channel do the next:
      
      # Background correction
      
      norm <- backgroundCorrect(rg(), method = 'normexp')
      
      # Between array normalization
      
      norm <- normalizeBetweenArrays(norm, method="quantile")
      
      return(norm)
      
    }else{
      
      print('two color normalization')
      
      
      # Background correction
      
      norm <- backgroundCorrect(rg(), method="normexp", offset=50)  
      
      # within array normalization
      
      norm <- normalizeWithinArrays(norm, method="loess")
      
      
      
      
      # between arrays normalization
      
      norm <- normalizeBetweenArrays(norm, method="Aquantile")
      
      print(class(norm))
      
      return(norm)
      
    }
    
  })
  
  
  # Definition of the outputs ----
  
  # Target example ----
  
  output$target_example <- renderTable({
    
    data.frame( FileName = c('Sample_1_AC.xlsx', 'Sample_2_AC.xlsx',
                             'Sample_3_cntr.xlsx', 'Sample_4_cntr.xlsx'),
                Cy3 = rep('Ref', 4),
                Cy5 = c('AC', 'AC', 'Control', 'Control'))
    
  })
  
  output$target_example_single <- renderTable({
    
    data.frame( FileName = c('Sample_1_AC.xlsx', 'Sample_2_AC.xlsx',
                             'Sample_3_cntr.xlsx', 'Sample_4_cntr.xlsx'),
                Cy3 = c('AC', 'AC', 'Control', 'Control'))
  })
  
  # Raw boxplot ----
  
  output$raw_box <- renderPlot({
    
    
    if(input$single_ch){
      
      boxplot(log2(rg()$E), main = 'Raw data boxplot',
              ylab = 'log2(Intensity)', xaxt='n',
              col = rainbow(ncol(rg()$E)))
      
      
    }else{

    boxplot(cbind(log2(rg()$G), log2(rg()$R)), main = 'Raw data boxplot',
            ylab = 'log2(Intensity)', xaxt='n',
            col = c(rep('green', dim(rg()$G)[2]), rep('red', dim(rg()$G)[2])))
      }
    })
  
  output$norm_box <- renderPlot({
    
    if(input$single_ch){
      
      boxplot(norm()$E, main = 'Normalized data boxplot',
              ylab = 'log2(Intensity)', xaxt='n',
              col = rainbow(ncol(norm()$E)))
      
    }else{
      
      boxplot(norm()$M, main = 'Normalized data boxplot 
              \n (M values)',
              ylab = 'log2(Intensity)', xaxt='n',
              col = rainbow(ncol(norm()$M)) )
      
      
    }
    
    
    
  })
  
  
  # Raw density plot -----
  
  output$density_raw <- renderPlot({
    
    plotDensities(rg(), main = 'Raw density plot')
  })
  

  # raw MA plot ----
  
  output$MA_raw <- renderPlot({
    limma::plotMA(rg(), main = 'MAplot of raw data')
  })
  
  # Normalized MA plot ----
  
  output$MA_norm <- renderPlot({
    
    limma::plotMA(norm(), main = 'MAplot of normalized data')
    
  })
  
  # normalized density plot ----
  
  output$density_norm <- renderPlot({
    
    print(class(norm()))
  
      
      limma::plotDensities(object = norm(),
                                  main = 'Normalized density plot')
 
  })
  
  # PCA analysis ----
  
  # First we store the PCa information in a variable
  
  pca.filt <- reactive({
    
    if(input$single_ch){
      
      prcomp(t(norm()$E), scale = TRUE )
      
    }else{
      
      prcomp(t(norm()$A), scale = TRUE )
      
    }
 
    
    })
  
  output$pca <- renderPlot({
    
    # First we choose the colors of the plot based on the targets argument
    
    if(input$single_ch){
      
      groups <- as.factor(rg()$targets$Cy3)
      
    }else{
      
      col_targets <- colnames(rg()$targets)
      
      
      
      no_select <- c(input$ref, 'FileName')
      
      
      
      sel_colunm <- col_targets[!(col_targets %in% no_select)]
      
      
      
      groups <- as.factor(rg()$targets[,sel_colunm])
      
    }
    
    
    # we choose as many colors as categories are
    
    colors <- rainbow(length(levels(groups)))
    
    
    
    group_color <- rep(NA, length(groups))
    
    
    
    for(i in 1:length(levels(groups))){
      
      
      group_color[groups == levels(groups)[i]] <- colors[i] 

      
    }
    
    # We plot the PCA depending on the type of PCA plot that we want
    
    
    if(input$pca_shape == '3D'){
      
    pca3d<-scatterplot3d(x=pca.filt()$x[,1],
                         y=pca.filt()$x[,2],
                         z=pca.filt()$x[,3],
                         xlab='PC1', ylab='PC2', zlab='PC3',
                         main='PCA Analysis',
                         pch=16,col.grid="lightblue",
                         color = group_color)
    
    }else{
    
    plot(x=pca.filt()$x[,1],y=pca.filt()$x[,2],
         xlab='PC1', ylab='PC2',
         main='PCA Analysis',
         pch=16,
         col = group_color)
      }
    
    # And finally we add a legend in order to know which color correspond with
    # each  group.
    
    legend('topright', legend = levels(groups),
           pch = 16, col = colors)
    
    
  })
  
  # plot a table with the groups that the user can use 
  # to make the comparison.
  
  output$groups <- renderTable({
    
    if(input$single_ch){
      
      groups <- as.factor(rg()$targets$Cy3)
      
    }else{
      
      col_targets <- colnames(rg()$targets)
      
      
      
      no_select <- c(input$ref, 'FileName')
      
      
      
      sel_colunm <- col_targets[!(col_targets %in% no_select)]
      
      
      
      groups <- as.factor(rg()$targets[,sel_colunm])
      
    }
    
    df <- data.frame(levels(groups))
    
    colnames(df) <- 'Groups'
    
    df
  })
  
  # Definition of the contrast matrix regarding if we have
  # a multiple comparison or a single comparison.
  
  DE_fit <- reactive({
    
    if(input$single_ch){
      
      design <- model.matrix(~ 0+factor(as.integer(
        as.factor(norm()$targets$Cy3)
                                                   )))
      
      colnames(design) <- levels(as.factor(norm()$targets$Cy3))
      
    }else{
      
      design <- modelMatrix(rg()$targets, ref = 'Ref')
      
    }
    
    
    fit <- lmFit(norm(), design)
    
    if(input$mul_comp){
      
      comparisons <- as.character(unlist(strsplit(input$contrasts,
                                                   split = ',')))
    
    }else{
      
      comparisons <- input$contrasts
      
      
    }
    
    contrast.matrix <- makeContrasts(contrasts = comparisons,
                                     levels = design)
    
    fit2 <- contrasts.fit(fit, contrast.matrix)
    
    fit2 <- eBayes(fit2)
    
    return(fit2)
    
  })
  
  
  ## Head toptable ----
  
  output$top_table <- renderTable({
    
    # if((input$single_ch & !(input$mul_comp))){
      
      # top_table <- toptable(DE_fit(), number = Inf,
      #                       genelist = DE_fit()$genes)
      # 
      # interest <- top_table[ abs(top_table$logFC) > input$fc &
      #                          top_table[,input$fdr_pvalue] < input$alpha,]
      # 
      # 
      # head(interest, n = 20)
      
    # } else {
      
      top_table <- toptable(DE_fit(), number = Inf,
                            coef = input$top_contrast,
                            genelist = DE_fit()$genes)
      
      interest <- top_table[ abs(top_table$logFC) > input$fc &
                               top_table[,input$fdr_pvalue] < input$alpha,]
      
      
      head(interest
           # [,c(7:12)]
           , n = 20)
    # }
    

    
    
  })
  
  # download table ---- 
  
  output$downloadData <- downloadHandler(
    
    filename = 'data.xlsx',

    
    content = function(file){
      
      if(input$mul_comp){
        
        comparisons <- as.character(unlist(strsplit(input$contrasts,
                                                    split = ',')))
        print(paste('comparisons: ', comparisons))
        
        for(i in comparisons){
          
          print(i)
          
          top_table <- toptable(DE_fit(), number = Inf,
                                coef = i,
                                genelist = DE_fit()$genes)
          
          interest <- top_table[abs(top_table$logFC) > input$fc &
                                top_table[,input$fdr_pvalue] < input$alpha,]
          
          
          write.xlsx(interest, file, row.names = FALSE,
                     sheetName = i, append = T)
          }
        
        
      }else{
        
        top_table <- toptable(DE_fit(), number = Inf,
                              coef = input$top_contrast,
                              genelist = DE_fit()$genes)
        
        interest <- top_table[ abs(top_table$logFC) > input$fc &
                                 top_table[,input$fdr_pvalue] < input$alpha,]
        
        write.xlsx(interest, file, row.names = FALSE,
                   sheetName = input$top_contrast)
      }
      
      
      
    }
    
  
  )
  
  # FASTA annotation section -----
  
  
  
  
   ###### Table with the groups -----
  
  annot <- reactive({
    
    annot <- read.table(file = input$annotation$datapath,
                        sep = input$del,
                        skip = input$skip,
                        header = T)
    
  })
  
  output$columns <- renderTable({


   data.frame(Columns = colnames(annot()))


  })
  
  ### FASTA dowload ----
  
  
  output$fasta_down <- downloadHandler(
    
    filename = 'DEGs_seqs.fasta',
    
    content = function(file){
      
      colseq <- input$col_seq
      
      print(colseq)
      
      colid <- input$col_id
      
      print (colid)
      
      
      top_table <- toptable(DE_fit(), number = Inf,
                            coef = input$group_fasta,
                            genelist = DE_fit()$genes)
      
      print(summary(annot()[,colid] %in% top_table[, input$top_id]))
      
      top_table <- top_table[ abs(top_table$logFC) > input$fc &
                               top_table[,input$fdr_pvalue] < input$alpha,]
      
      degs <- top_table[, input$top_id]
      
      print(degs)
      
      
      
      
      
      deg_seq <- annot()[annot()[,colid] %in% degs,]
      
      print(summary(annot()[,colid] %in% degs))
      
      print(head(deg_seq, n= 20))
      
      
      
      deg_seq <- deg_seq[!(deg_seq[,colseq] == ''),]
      
      
      print(summary(deg_seq[,colseq] == ''))
      
      return(seqinr::write.fasta(sequences = as.list(as.character(deg_seq[,colseq])),
                  names = as.character(deg_seq[,colid]),
                  file.out = file,
                  open = 'w'))
  
      
    }
    
    
  )
}

# Run the app
shinyApp(ui, server)
