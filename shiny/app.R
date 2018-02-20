library(shiny)
library(shinyRGL)
library(rgl)
library(data.table)

setwd("/Users/zhenyu/github/tsne/shiny")
tsne = read.table("combined.tsne3.PC200.dim3.3_3_1_1.perm1000.txt", h=T, stringsAsFactors=T)
projects = unique(tsne$project)
cols = rainbow(length(projects))
tsne$col = cols[tsne$project]
tsne$alpha = 0.8
tsne$radius = 2


ui <- fluidPage(
  headerPanel("Molecular Profiling of PanCancer Gene Expression/ miRNA Expression/ Methylation/ Copy Number Variation"),
  sidebarPanel(
    checkboxGroupInput( "project_selected", 
                        "Cancer Type", 
                        choices = projects, 
                        selected = projects)),
  mainPanel(
#    tableOutput("data"))
    plotOutput("rglPlot"))
)

server <- function(input, output) {
#  data <- reactive({
#    tsne3[which(tsne3$project %in% input$project_selected), ]
#    new.tsne = tsne3
#    new.tsne[!new.tsne$project %in% input$project_selected, "col"] = "#55555555"
#    new.tsne
#  })
  output$rglPlot <- renderPlot({
    try(rgl.close())
    project_selected = input$project_selected
    data = tsne
    if(length(project_selected) != length(projects)) {
      w = which(!data$project %in% project_selected)
      data$col[w] = "#333333FF"
      data$radius[w] = data$radius[w]/2
      data$alpha[w] = 0.2
    }
    plot3d(data$feature_1, data$feature_2, data$feature_3, col=data$col, alpha=data$alpha, add=F,type="s",radius=data$radius)
    box3d(labels=FALSE, tick=FALSE, box=FALSE)
    rglwidget()
  })
}

shinyApp(ui = ui, server = server)

library(shiny)
library(shinyRGL)
library(rgl)

ui <- fluidPage(
    
    titlePanel("GDC Shiny app :: Load and Render Previously Calculated PCoA with Corresponding Metadata"), 
  
    fluidRow(
        splitLayout(cellWidths = c("4in", "4in", "12in"),
                    
                    sidebarPanel(
                        width=12,
                        selectInput(inputId="pcoa_in", label="select PCoA file (*.PCoA)", choices=dir(pattern=".PCoA$")),
                        actionButton(inputId="load_PCoA.button", label="Load PCoA"),
                        selectInput(inputId="metadata_in", label="select metadata file (*.metadata.txt)", choices=dir(pattern=".txt$")),
                        actionButton(inputId="load_metadata.button", label="Load metadata"),
                        uiOutput('PCO_1'),
                        uiOutput('PCO_2'),
                        uiOutput('PCO_3'),
                        uiOutput('metadata_select'),
                        actionButton(inputId="render_pcoa.button", label="render PCoA"),
                        textOutput('pcoa_in'),
                        textOutput('metadata_in'),
                        #actionButton(inputId="reset.button", label="Reset"), # coming soon
                        textOutput('blank_1'),
                        textOutput('blank_2'),
                        textOutput('blank_3'),
                        textOutput('blank_4'),
                        textOutput('blank_5'),
                        textOutput('blank_6')
                    ),
                    
                    
                    plotOutput("my_legend", height="36in"),
                    
                    webGLOutput("myWebGL", width="12in", height="12in"))
    )
        
)




server <- function(input, output){
    
    import_data <- function(file_name){
        data.matrix(read.table(file_name, row.names=1, header=TRUE, sep="\t", comment.char="", quote="", check.names=FALSE))
    }  
    
    import_metadata <- function(group_table){ #, group_column, sample_names){
        metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
            read.table(
                file=group_table,row.names=1,header=TRUE,sep="\t",
                colClasses = "character", check.names=FALSE,
                comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
            )
        )
    }
    
    load_pcoa_data.v2 <- function(PCoA_in){
        con_1 <- file(PCoA_in)
        con_2 <- file(PCoA_in)
                                        # read through the first time to get the number of samples
        open(con_1);
        num_values <- 0
        data_type = "NA"
        while ( length(my_line <- readLines(con_1,n = 1, warn = FALSE)) > 0) {
            if ( length( grep("PCO", my_line) ) == 1  ){
                num_values <- num_values + 1
            }
        }
        close(con_1)
                                        # create object for values
        eigen_values <- matrix("", num_values, 1)
        dimnames(eigen_values)[[1]] <- 1:num_values
        eigen_vectors <- matrix("", num_values, num_values)
        dimnames(eigen_vectors)[[1]] <- 1:num_values
                                        # read through a second time to populate the R objects
        value_index <- 1
        vector_index <- 1
        open(con_2)
        current.line <- 1
        data_type = "NA"
        while ( length(my_line <- readLines(con_2,n = 1, warn = FALSE)) > 0) {
            if ( length( grep("#", my_line) ) == 1  ){
                if ( length( grep("EIGEN VALUES", my_line) ) == 1  ){
                    data_type="eigen_values"
                } else if ( length( grep("EIGEN VECTORS", my_line) ) == 1 ){
                    data_type="eigen_vectors"
                }
            }else{
                split_line <- noquote(strsplit(my_line, split="\t"))
                if ( identical(data_type, "eigen_values")==TRUE ){
                                        #dimnames(eigen_values)[[1]][value_index] <- noquote(split_line[[1]][1])
                    dimnames(eigen_values)[[1]][value_index] <- gsub("\"", "", noquote(split_line[[1]][1]))
                    eigen_values[value_index,1] <- noquote(split_line[[1]][2])       
                    value_index <- value_index + 1
                }
                if ( identical(data_type, "eigen_vectors")==TRUE ){
                                        #dimnames(eigen_vectors)[[1]][vector_index] <- noquote(split_line[[1]][1])
                    dimnames(eigen_vectors)[[1]][vector_index] <- gsub("\"", "", noquote(split_line[[1]][1]))
                    for (i in 2:(num_values+1)){
                        eigen_vectors[vector_index, (i-1)] <- as.numeric(noquote(split_line[[1]][i]))
                    }
                    vector_index <- vector_index + 1
                }
            }
        }
        close(con_2)
                                        # finish labeling of data objects
        dimnames(eigen_values)[[2]] <- "EigenValues"
        dimnames(eigen_vectors)[[2]] <- dimnames(eigen_values)[[1]]
        class(eigen_values) <- "numeric"
        class(eigen_vectors) <- "numeric"
                                        # write imported data to global objects
        
        {return(list(eigen_values=eigen_values, eigen_vectors=eigen_vectors))}
        
    }
    
                                        # import *.PCoA (matR formatted)
    observeEvent(input$load_PCoA.button, {
        
                                        #output$pcoa_in <- renderText( paste(input$pcoa_in, "LOADING", sep ="  ::  ") )
        my_pcoa <<- load_pcoa_data.v2(input$pcoa_in)
        output$pcoa_in <- renderText( paste(input$pcoa_in, "LOADED", sep ="  ::  ") )  

        output$blank_1 <- renderText( "" )
        output$blank_2 <- renderText( "" )
        output$blank_3 <- renderText( "" )
        output$blank_4 <- renderText( "" )
        output$blank_5 <- renderText( "" )
        output$blank_6 <- renderText( "" )
        
        
                                        # load PCS to the interface (just do 3 to keep it simple)
        output$PCO_1 = renderUI({
            selectInput(inputId='PCO_1', label='PC0_1', choices=(colnames(my_pcoa$eigen_vectors)), selected=colnames(my_pcoa$eigen_vectors)[1])
        })
        
        output$PCO_2 = renderUI({
            selectInput(inputId='PCO_2', label='PC0_2', choices=(colnames(my_pcoa$eigen_vectors)), selected=colnames(my_pcoa$eigen_vectors)[2])
        })
        
        output$PCO_3 = renderUI({
            selectInput(inputId='PCO_3', label='PC0_3', choices=(colnames(my_pcoa$eigen_vectors)), selected=colnames(my_pcoa$eigen_vectors)[3])
        })
        
    })
    
                                        # import metadata (matR formatted) and
    observeEvent(input$load_metadata.button, {
        
                                        #output$metadata_in <- renderText( paste(input$metadata_in, "LOADING", sep="  ::  ") )
        my_metadata <<- import_metadata( input$metadata_in )
        output$metadata_in <- renderText( paste(input$metadata_in, "LOADED", sep="  ::  ") )
        
        output$metadata_select = renderUI({
            selectInput(inputId='metadata_select', label='Select Metadata', choices=(colnames(my_metadata)), selected=colnames(my_metadata)[1])
        })

        output$blank_1 = renderText({""})
        output$blank_2 = renderText({""})
        output$blank_3 = renderText({""})
        output$blank_4 = renderText({""})
        output$blank_5 = renderText({""})
        output$blank_6 = renderText({""})
        
    })
    
    
    observeEvent(input$render_pcoa.button, {

        weight_distances = TRUE # if TRUE , x y and z will have the same range regardless of weight
        
        ######################
        # SUB( ): Create optimal contrast color selection using a color wheel
        # adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html 
        ######################
        col.wheel <- function(num_col, my_cex=0.75) {
            cols <- rainbow(num_col)
            col_names <- vector(mode="list", length=num_col)
            for (i in 1:num_col){
                col_names[i] <- getColorTable(cols[i])
            }
            cols
        }
        ######################
        ######################
        
        ######################
        # SUB( ): Automtically generate colors from metadata with identical text or values
        ######################
        create_colors <- function(metadata_column, color_mode = "auto"){ # function to     
            my_data.color <- data.frame(metadata_column)

            #my_data.color <- my_data.color[ order(as.numeric(my_data.color)), ]
            
            column_factors <- as.factor(metadata_column[,1])
            column_levels <- levels(as.factor(metadata_column[,1]))
            num_levels <- length(column_levels)
            color_levels <- col.wheel(num_levels)
            levels(column_factors) <- color_levels
            my_data.color[,1]<-as.character(column_factors)
            
            return(my_data.color)
        }
        ######################
        ######################

        ######################
        # SUB(9): The inverse function to col2rgb()
        # adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html
        ######################
        rgb2col <- function(rgb) {
            rgb <- as.integer(rgb)
            class(rgb) <- "hexmode"
            rgb <- as.character(rgb)
            rgb <- matrix(rgb, nrow=3)
            paste("#", apply(rgb, MARGIN=2, FUN=paste, collapse=""), sep="")
        }
        ######################
        ######################

        ######################
        # SUB(10): Convert all colors into format "#rrggbb"
        # adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html
        ######################
        getColorTable <- function(col) {
            rgb <- col2rgb(col);
            col <- rgb2col(rgb);
            sort(unique(col))
        }
        ######################
        ######################

        ######################
        # SUB(4): Sub to provide scaling for title and legened cex
        ######################
        calculate_cex <- function(my_labels, my_pin, my_mai, reduce_by=0.30, debug){
            
                                        # get figure width and height from pin
            my_width <- my_pin[1]
            my_height <- my_pin[2]
            
                                        # get margine from mai
            my_margin_bottom <- my_mai[1]
            my_margin_left <- my_mai[2]
            my_margin_top <- my_mai[3]
            my_margin_right <- my_mai[4]
            
                                        # find the longest label (in inches), and figure out the maximum amount of length scaling that is possible
            label_width_max <- 0
            for (i in 1:length(my_labels)){  
                label_width <- strwidth(my_labels[i],'inches')
                if ( label_width > label_width_max){ label_width_max<-label_width  }
            }
            label_width_scale_max <- ( my_width - ( my_margin_right + my_margin_left ) )/label_width_max
            
                                        # find the number of labels, and figure out the maximum height scaling that is possible
            label_height_max <- 0
            for (i in 1:length(my_labels)){  
                label_height <- strheight(my_labels[i],'inches')
                if ( label_height > label_height_max){ label_height_max<-label_height  }
            }
            adjusted.label_height_max <- ( label_height_max + label_height_max*0.4 ) # fudge factor for vertical space between legend entries
            label_height_scale_max <- ( my_height - ( my_margin_top + my_margin_bottom ) ) / ( adjusted.label_height_max*length(my_labels) )
            
                                        # max possible scale is the smaller of the two 
            scale_max <- min(label_width_scale_max, label_height_scale_max)
                                        # adjust by buffer
        ###scale_max <- scale_max*(100-buffer/100) 
            adjusted_scale_max <- ( scale_max * (1-reduce_by) )
                                        #if(debug==TRUE){ print(cat("\n", "adjusted_scale_max: ", adjusted_scale_max, "\n", sep=""))  }
            return(adjusted_scale_max)   
        }
        ######################
        ######################
        
        ######################
        # SUB(5): Fetch par values of the current frame - use to scale cex
        ######################
        par_fetch <- function(){
            my_pin<-par('pin')
            my_mai<-par('mai')
            my_mar<-par('mar')
            return(list("my_pin"=my_pin, "my_mai"=my_mai, "my_mar"=my_mar))    
        }
        ######################
        ######################


        ## # generate key and colors from seleced metadata
        metadata_column <- my_metadata[ , input$metadata_select , drop=FALSE ] # get column i from the metadata matrix
        
        sample_names <- rownames(my_pcoa$eigen_vectors)
        
        suppressWarnings( numericCheck <- as.numeric(metadata_column) ) # check to see if metadata are numeric, and sort accordingly # Check this block - not sure what it is doing : Kevin 11-10-16
        if( is.na(numericCheck[1])==FALSE ){
            column_name = colnames(metadata_column)[1]
            row_names = rownames(metadata_column)
            metadata_column <- matrix(numericCheck, ncol=1)
            colnames(metadata_column) <- column_name
            rownames(metadata_column) <- row_names
        }
        
        color_column <- create_colors(metadata_column, color_mode = "auto") # set parameters for plotting
        
        ########## ncol.color_matrix <- 1 
        column_factors <- as.factor(metadata_column) 
        column_levels <- levels(as.factor(metadata_column))
        num_levels <- length(column_levels)
        color_levels <- col.wheel(num_levels)
        suppressWarnings( column_levels <- column_levels[ order(as.numeric(column_levels)) ] ) # NEW (order by levels values) # was using order, got an alphanumeric sort
        suppressWarnings( color_levels <- color_levels[ order(as.numeric(column_levels)) ] ) # NEW (order by levels values)
        
        ## # make sure that vectors are sorted identically to the colors
        rownames(my_pcoa$eigen_vectors) <- gsub("\"", "", rownames(my_pcoa$eigen_vectors)) 
        my_pcoa$eigen_vectors <- my_pcoa$eigen_vectors[ rownames(color_column), ]
        
        ## # get the colors list form the column
        pcoa_colors <- as.character(color_column[,1])
        
        pcx <- my_pcoa$eigen_vectors[,input$PCO_1]
        pcy <- my_pcoa$eigen_vectors[,input$PCO_2]
        pcz <- my_pcoa$eigen_vectors[,input$PCO_3]
        cols <- pcoa_colors

                                        # get (normalized) eigen values for the axist labels
        names(my_pcoa$eigen_values) <- colnames(my_pcoa$eigen_vectors)
        x_eigen <- paste( "PC_", input$PCO_1, " :: ", round( my_pcoa$eigen_values[input$PCO_1], 4 ), sep="" )
        x_eigen.TEST <<- paste( "PC_", input$PCO_1, " :: ", round( my_pcoa$eigen_values[input$PCO_1], 4 ), sep="" )
        y_eigen <- paste( "PC_", input$PCO_2, " :: ", round( my_pcoa$eigen_values[input$PCO_2], 4 ), sep="" )
        z_eigen <- paste( "PC_", input$PCO_3, " :: ", round( my_pcoa$eigen_values[input$PCO_3], 4 ), sep="" )

        eigen_vectors.TEST <<- my_pcoa$eigen_vectors
        eigen_values.TEST <<- my_pcoa$eigen_values

        output$myWebGL <- renderWebGL({
            plot3d(pcx, pcy, pcz, col=cols, alpha=0.6, add=T,type="s",radius=.000025)
            #legend3d("topright",legend=c(column_levels), col=c(color_levels), pch=16, cex=0.8, inset=c(0.05), border=NA, bty="n")
            if ( weight_distances == TRUE ){
                aspect3d(1,1,1)
            }else{
            }
            box3d(labels=FALSE, tick=FALSE, box=FALSE)
            decorate3d(xlab=x_eigen, ylab=y_eigen, zlab=z_eigen, axes=FALSE, cex=1)

        })


        #output$my_legend <- renderPlot ({ legend("center", legend=column_levels, pch=16, col=color_levels, cex=0.8) })
        output$my_legend <- renderPlot ({
            plot.new()
            legend_par <- par_fetch()
            legend_cex <- calculate_cex(column_levels, legend_par$my_pin, legend_par$my_mai, reduce_by=0.20)
            legend_cex.TEST <<- legend_cex
            legend("top", legend=column_levels, pch=16, col=color_levels, cex=legend_cex)
        })
        
    })

    observeEvent(input$reset.button, {})
        
}
                              





