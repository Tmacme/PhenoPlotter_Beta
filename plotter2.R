# plotting tools
################ FUNCTIONS

library(ggplot2)
library(readxl)
library(stringr)
library(dplyr)
library(reshape2)
library(lubridate)
library(gdata)
library(readr)
library(tidyr)
library(Hmisc)
library(plyr)
library(zoo)
library(signal)


getProcedureTable <- function(procedure_code, directory){
  dataDir <- paste(directory,procedure_code,"\\ScrapedProcedureTable",sep="")
  setwd(dataDir)
  # will get ONLY THE FIRST procedure table in the directory <- Needs to be improved to take into account multiple procedure tables
  files <- list.files(path=dataDir, pattern = "*.csv",full.names=TRUE, recursive = FALSE)
  procedure_table <- read_delim(files[1], delim = ',')
  colnames(procedure_table)[1] <- "all_titles"
  procedure_table$all_titles <- make.names(procedure_table$all_titles)
  return(procedure_table)
}

getIMPCDictionary <- function(directory){
  setwd(directory)
  IMPC_dictionary <- read_delim("IMPC_DataDictionary.csv", delim = ',')
  IMPC_dictionary$Parameter <- make.names(IMPC_dictionary$Parameter)
  return(IMPC_dictionary)
}

getWideCSVFiles <- function(procedure_code, directory){
  dataDir <- paste(directory,procedure_code,"\\WideData",sep="")
  setwd(dataDir)
  files <- list.files(path=dataDir, pattern = "*.csv",full.names=TRUE, recursive = FALSE)
  return(files)
}

getCenterData <- function(filename, procedure_code = ""){
  file_extension <- strsplit(filename, "/")[[1]][2]
  center_data <- read_delim(filename, delim = ',')
  # center_data <- read.csv(filename, stringsAsFactors = F)
  if('Start_hms_new' %in% colnames(center_data)){
    center_data <- read.csv(filename,colClasses = c("Start_hms_new" = "character"))
  }
  colnames(center_data) <- make.names(colnames(center_data))
  if(procedure_code == "OFD" && is.element("Center.distance.travelled", colnames(center_data)) && is.element("Periphery.distance.travelled", colnames(center_data))){
    center_data$Total.distance.travelled <- center_data$Center.distance.travelled + center_data$Periphery.distance.travelled
  }
  return(center_data)
}

getCenterName <- function(filename){
  file_extension <- strsplit(filename, "/")[[1]][2]
  center_name <- strsplit(file_extension, "_")[[1]][1]
}

getPhenotypes <- function(procedure_table, IMPC_dictionary, center_data, procedure_code){
  
  
  # IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol == procedure_code,]
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Phenotype == 'yes',]
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Use == 'yes',]
  # gets list of desired subset of phenotypes from the DataDictionary
  index = 1
  
  desired_phenotypes <- unique(IMPC_dictionary$Parameter)
  
  # gets list of phenotypes with data
  index = 1
  
  phenotype_list <- list()
  
  if(procedure_code == "OFD"){
    phenotype_list <- list("Total.distance.travelled")
  }
  for(i in desired_phenotypes){
    if(is.element(i,colnames(center_data))){
      unique_vals <- length(unique(center_data[[i]]))
      
      if(any(is.na(center_data[[i]]))){
        unique_vals <- unique_vals - 1
      }
      
      if(unique_vals>1 && !is.null(center_data[[i]])){
        phenotype_list[[i]] <- i
      }
    }
  }
  return(phenotype_list)
}

getMetadata <- function(procedure_table, IMPC_dictionary, center_data, procedure_code, pipeline = "IMPC"){
  full_metadata_list <- list()
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol==procedure_code,]
  # extract phenotype and metadata lists from procedure_table
  index = 1
  on_metadata <- FALSE
  for(i in 1:(nrow(procedure_table)-1)){
    if(procedure_table$all_titles[i] == 'NA.'){
      index = 1
      on_metadata<-TRUE
    }
    if(on_metadata){
      full_metadata_list[[index]] <- procedure_table$all_titles[i+1]
    }      
    index <- index + 1 
  }
  full_metadata_list <- lapply(full_metadata_list, function(x) strsplit(x, paste(".",pipeline,sep=""))[[1]][1])

  # gets list of desired subset of metadata from the DataDictionary
  index = 2
  desired_metadata <- list('production_phenotype')
  for(i in full_metadata_list){
    if(is.element(i,IMPC_dictionary$Parameter)){
      row_num <- which(IMPC_dictionary$Parameter==i)
      # check if "use" is marked in the data dictionary
      if(!is.na(IMPC_dictionary$Use[row_num]) && IMPC_dictionary$Use[row_num]=="yes" && IMPC_dictionary$Field[row_num]=="factor"){
        desired_metadata[index] <- i
        index <- index + 1
      }
    }
  }
  
  # gets list of plottable metadata
  index = 1
  metadata_list <- list()
  for(i in desired_metadata){
    if(is.element(i,colnames(center_data))){
      unique_vals <- length(unique(center_data[[i]]))
      if(any(is.na(center_data[[i]]))){
        unique_vals <- unique_vals - 1
      }
      
      if(unique_vals>1 && unique_vals<100 && i!="Start.Time"){
        # print(paste("metadata = ",i," ; unique_vals = ", unique(center_data[[i]])[1]," ", unique(center_data[[i]])[2], sep = ""))
        metadata_list[index] <- i
        index <- index + 1
      }
    }
  }
  return(metadata_list)
}

getTimeSeries <- function(IMPC_dictionary, center_data, procedure_code){
  timeseries <- list("date_of_experiment")
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol==procedure_code,]
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Field=="time",]
  timeseries <- c(timeseries, IMPC_dictionary$Parameter)
  return(timeseries)
}

produceBoxPlot <-function(center_data, phenotype, factors, title = ""){
  x <- factors
  y <- phenotype
  center_data <- center_data[complete.cases(center_data[[x]]),]
  center_data <- center_data[complete.cases(center_data[[y]]),]
  center_mean <- mean(center_data[[y]], na.rm = TRUE)
  center_std <- sd(center_data[[y]], na.rm = TRUE)
  if(length(unique(center_data[[x]]))>1){
    print(paste(phenotype,factors,sep=" "))
    plot <- ggplot() + 
      geom_point(aes(x = as.factor(center_data[[x]]),y = center_data[[y]] ,colour = sex),data=center_data,shape = 16,size = 1.0, position = position_jitter(w = 0.2, h = 0)) +  
      scale_colour_brewer(guide = guide_legend(),palette = 'Set1') +
      geom_crossbar(aes(y = center_data[[y]],x = reorder(as.factor(center_data[[x]]),center_data[[y]], FUN = mean, na.rm=TRUE), na.rm=TRUE),data=subset(center_data, !is.na(center_data[[y]])),na.rm = TRUE, colour = "Blue", width = 0.5, fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
      theme_bw() + labs(x=x, y=y) + ggtitle(title) +
      geom_hline(yintercept = center_mean, color = "grey", size = 1) + geom_hline(yintercept = center_mean + center_std, color = "grey", size = 1) + geom_hline(yintercept = center_mean - center_std, color = "grey", size = 1) +
      theme(axis.text=element_text(size=6),axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
    return(plot)
  }
}

produceAllBoxPlot <- function(procedure_code, directory, isControl = TRUE, pipeline = "IMPC"){
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  IMPC_dictionary <- IMPC_dictionary[IMPC_dictionary$Protocol == procedure_code,]
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  if(isControl){
    control_or_not <- "Controls"
  }else{
    control_or_not <- "All"
  }
  pdf(paste(date,control_or_not,"_BOXPLOT_",procedure_code,".pdf", sep = ""), 6,6)
  for(filename in files){
    center_name <- getCenterName(filename)
    print(paste("STARTED ", center_name))
    center_data <- getCenterData(filename,procedure_code)
    phenotype_list <- getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    metadata_list <- getMetadata(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code, pipeline = pipeline) 
    if(isControl){
      center_data <- center_data[center_data$biological_sample_group == 'control',]
    }
    for(metadata in metadata_list){
      if(is.element(metadata, colnames(center_data))){
        # print(paste("METADATA2: ", metadata,sep=""))
        
        for(phenotype in phenotype_list){
          if(is.element(phenotype, colnames(center_data))){
            # print(paste("METADATA2x: ", metadata,sep=""))
            
            if(length(unique(center_data[[metadata]]))>1){
              # print(paste("METADATA3: ", metadata,sep=""))
              
              row_num <- which(IMPC_dictionary$Parameter==metadata)
              
              if(!is.na(IMPC_dictionary$MetadataCategory[row_num])){
                metadata_split <- IMPC_dictionary$MetadataCategory[row_num]
                title <- paste(center_name,metadata,";",phenotype,"\n","(",control_or_not,";",metadata_split,")",sep=" ")
              }else{
                title <- paste(center_name,metadata,";",phenotype,control_or_not,sep=" ")
              }
              # print(paste("HERE IS METADATA: ", metadata))
              boxplot <- produceBoxPlot(center_data = center_data, phenotype = phenotype, factors = metadata, title = title)
              if(!is.null(boxplot)){
                print(boxplot)
              }
            }
          }
        }
      }
    }
    print(paste("COMPLETED BOXPLOT: ",control_or_not, " ",center_name,sep=""))
  }
  dev.off()
}

produceTimeSeries <- function(center_data, phenotype, timeplot, title = "", moving_average = 1){
  x <- timeplot
  y <- phenotype
  center_data <- center_data[complete.cases(center_data[[x]]),]
  center_data <- center_data[complete.cases(center_data[[y]]),]
  center_mean <- mean(center_data[[y]], na.rm = TRUE)
  center_std <- sd(center_data[[y]], na.rm = TRUE)
  center_data[[x]] <- parse_date_time(center_data[[x]],orders = c("Ymd", "Ymd H:M:S", "Y/m/d", "m/d/Y", "mdy", "mdy H:M","mdY"))
  center_data[[x]] <- as.Date(center_data[[x]])
  plot <- ggplot() +
    geom_point(aes(x = center_data[[x]],y = center_data[[y]],colour = sex),data=center_data,shape = 16,size = 1.0) +
    geom_line(aes(x = center_data[[x]],y = center_data[[y]]),data=center_data,size = 1.0,fun.data = mean_sdl,stat = 'summary') +
    scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + labs(x=x, y=y) + ggtitle(title) +
    stat_smooth(data=center_data, geom="smooth", position="") +
    geom_hline(yintercept = center_mean, color = "grey", size = 1) + geom_hline(yintercept = center_mean + center_std, color = "grey", size = 1) + geom_hline(yintercept = center_mean - center_std, color = "grey", size = 1) +
    theme(axis.text=element_text(size=6),axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
  return(plot)
}

produceAllTimeSeries <- function(procedure_code, directory, isControl = TRUE, moving_average = 1){
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  if(isControl){
    control_or_not <- "Controls"
  }else{
    control_or_not <- "All"
  }
  pdf(paste(date, control_or_not,"_TIMESERIES_",procedure_code,".pdf", sep = ""), 14,4)
  for(filename in files){
    center_data <- getCenterData(filename, procedure_code)
    center_name <- getCenterName(filename)
    phenotype_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    timeplot_list <- getTimeSeries(IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    
    if(isControl){
      center_data <- center_data[center_data$biological_sample_group == 'control',]
    }
    if(nrow(center_data)>1){
      for(phenotype in phenotype_list){
        if(is.element(phenotype, colnames(center_data))){
          for(timeplot in timeplot_list){
            if(is.element(timeplot, colnames(center_data))){
              title <- paste(center_name,phenotype,";", timeplot,"\n","(",control_or_not,")",sep=" ")
              time_series_plot <- produceTimeSeries(center_data = center_data, phenotype = phenotype, timeplot = timeplot, title = title, moving_average = moving_average)
              # tryCatch(print(time_series_plot),error=function(x){print(paste(center_name, phenotype, timeplot,sep=" "))})
              print(time_series_plot)
            }
          }
        }
      }
    }
  }
  dev.off()
}

produceMultiTimeSeries <- function(center_data, phenotype, timeplot, factors, title = "", moving_average = 1, isLine = TRUE){
  x <- timeplot
  y <- phenotype
  z <- factors
  center_data <- center_data[complete.cases(center_data[[x]]),]
  center_data <- center_data[complete.cases(center_data[[y]]),]
  center_data <- center_data[complete.cases(center_data[[z]]),]
  center_mean <- mean(center_data[[y]], na.rm = TRUE)
  center_std <- sd(center_data[[y]], na.rm = TRUE)
  center_data[[x]] <- parse_date_time(center_data[[x]],orders = c("Ymd", "Ymd H:M:S", "Y/m/d", "m/d/Y", "mdy", "mdy H:M","mdY"))
  
  if(length(unique(center_data[[z]]))>1){
    geom <- list()
    
    temp <- data.frame()
    for(val in unique(center_data[[z]])){
      if(sum(center_data[[z]]==val)>1){
        val_data <- center_data[center_data[[z]] == val,]
        temp <- rbind(temp, data.frame(x = val_data[[x]],y = val_data[[y]],color=toString(val)))
      }
    }
    
    
    if(isLine){
      g <- geom_line(data=temp,aes(x=x,y=y,col=color,group=color),na.rm=T)
    }else{
      g <- geom_point(data=temp,aes(x=x,y=y,col=color),na.rm=T, alpha=0.5)
    }
    geom <- c(geom, g)
    
    # PREVIOUS ATTEMPT AT ROLLING MEAN
    # center_zoo <- zoo(val_data[[y]], val_data[[x]])
    # m_av <- rollmean(center_zoo, moving_average, fill=list(NA,NULL,NA))
    # temp <- data.frame(x = val_data[[x]],y = coredata(m_av), color = toString(val))
    # if(isLine){
    #   g <- geom_line(data=temp,aes(x=x,y=y,col=color),na.rm=T)
    # }else{
    #   g <- geom_point(data=temp,aes(x=x,y=y,col=color),na.rm=T, alpha=0.5)
    # }  , group = combined_df$phenotyping_center
    
    
    plot <- ggplot() + unlist(geom) + labs(x=x,y=y,colour=z) + ggtitle(title) +
      geom_hline(yintercept = center_mean, color = "grey", size = 1) + geom_hline(yintercept = center_mean + center_std, color = "grey", size = 1) + geom_hline(yintercept = center_mean - center_std, color = "grey", size = 1) +
      theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
    return(plot)
  }else{
    return("")
  }
}

produceAllMultiTimeSeries <- function(procedure_code, directory, isControl = TRUE, moving_average = 1, isLine = TRUE, pipeline = "IMPC"){
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  if(isControl){
    control_or_not <- "Controls"
  }else{
    control_or_not <- "All"
  }
  if(!isLine){
    line_or_not <- "ASPOINTS"
    pdf_name <- paste(date, control_or_not,"_ASPOINTS_MULTITIMESERIES_",procedure_code,".pdf", sep = "")
  }else{
    line_or_not <- "ASLINE"
    pdf_name <- paste(date, control_or_not,"_MULTITIMESERIES_",procedure_code,".pdf", sep = "")
  }
  
  pdf(pdf_name, 14,4)
  for(filename in files){
    center_data <- getCenterData(filename, procedure_code)
    center_name <- getCenterName(filename)
    
    phenotype_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    timeplot_list <- getTimeSeries(IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    factor_list <- getMetadata(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code, pipeline = pipeline) 
    if(isControl){
      center_data <- center_data[center_data$biological_sample_group == 'control',]
    }
    for(phenotype in phenotype_list){
      if(is.element(phenotype, colnames(center_data))){
        for(timeplot in timeplot_list){
          if(is.element(timeplot, colnames(center_data))){
            for(factors in factor_list){
              if(is.element(factors, colnames(center_data))){
                row_num <- which(IMPC_dictionary$Parameter==factors)
                if(length(row_num)>1){
                  if(!is.na(IMPC_dictionary$MetadataCategory[row_num])){
                    metadata_split <- IMPC_dictionary$MetadataCategory[row_num]
                  }
                  title <- paste(center_name,phenotype,";", timeplot,";",factors,"\n","(",control_or_not,";",metadata_split,")",sep=" ")
                  multitime_series_plot <- produceMultiTimeSeries(center_data = center_data, phenotype = phenotype, timeplot = timeplot, factors = factors, title = title, moving_average = moving_average, isLine = isLine)
                  if(!multitime_series_plot == ""){
                    print(multitime_series_plot)
                  }
                }
              }
            }
          }
        }
      }
    }
    
    print(paste("COMPLETED MULTITIMESERIES: ",control_or_not," ",line_or_not," ",center_name,sep=""))
    
  }
  dev.off()
}

produceAllelePlot <- function(center_data, phenotype, title = "", isSorted = FALSE, center_name = "", z_score_plots = FALSE){
  x <- "allele_symbol_zygosity"
  y <- phenotype
  
  
  center_data$allele_symbol <- as.character(center_data$allele_symbol)
  center_data$strain_name <- as.character(center_data$strain_name)
  center_data$allele_symbol[(center_data$allele_symbol == "" | is.na(center_data$allele_symbol)) & center_data$biological_sample_group == 'control'] <- center_data$strain_name[(center_data$allele_symbol == "" | is.na(center_data$allele_symbol)) & center_data$biological_sample_group == 'control']
  
  # center_data$zygosity <- " "
  center_data$allele_symbol_zygosity <- paste(center_data$allele_symbol,center_data$zygosity,sep=" ")
  controls_only <- center_data[center_data$biological_sample_group == 'control',]
  controlMean <- mean(controls_only[[y]], na.rm = T)
  controlSD <- sd(controls_only[[y]], na.rm = T)
  center_data$z_score <- (center_data[[y]]-controlMean) / controlSD
  if(z_score_plots){
    plot <- ggplot() +
      geom_pointrange(aes(x = reorder(center_data[[x]], center_data$z_score, FUN = mean,na.rm=TRUE),y = center_data$z_score, na.rm = TRUE),data= center_data,colour = '#0000ff',fill = '#ffffff',size = 0.8, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
      geom_point(aes(x = center_data[[x]],y = center_data$z_score, na.rm = TRUE),data=center_data,shape = 20,colour = '#999999',size = 1.0, alpha = 0.3) +
      theme_bw(base_size = 4.0) + labs(title = title) +
      geom_hline(yintercept = -1, color  = "red") + 
      geom_hline(yintercept = 0, color  = "red") +
      geom_hline(yintercept = 1, color  = "red") +
      labs(x=x,y=y) +
      coord_flip() +
      theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
  }else{
    if(isSorted){
      plot <- ggplot() +
        geom_pointrange(aes(x = reorder(center_data[[x]], center_data[[y]], FUN = mean,na.rm=TRUE),y = center_data[[y]], na.rm = TRUE),data= center_data,colour = '#0000ff',fill = '#ffffff',size = 0.8, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
        geom_point(aes(x = center_data[[x]],y = center_data[[y]], na.rm = TRUE),data=center_data,shape = 20,colour = '#999999',size = 1.0, alpha = 0.3) +
        scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + theme_bw(base_size = 4.0) + labs(title = title) +
        geom_hline(yintercept = controlMean, color  = "red") + 
        geom_hline(yintercept = controlMean + controlSD, color  = "red") +
        geom_hline(yintercept = controlMean - controlSD, color  = "red") +
        labs(x=x,y=y) +
        coord_flip() +
        theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
      return(plot)
    }else{
      plot <- ggplot() +
        geom_point(aes(x = center_data[[x]],y = center_data[[y]], na.rm = TRUE),data= center_data,shape = 20,colour = '#999999',size = 1.0, alpha = 0.3) +
        geom_pointrange(aes(x = center_data[[x]],y = center_data[[y]], na.rm = TRUE),data= center_data,colour = '#0000ff',fill = '#ffffff',size = 0.8, shape = 20,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary',na.rm=T) + 
        scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + theme_bw(base_size = 4.0) + labs(title = title) + 
        geom_hline(yintercept = controlMean, color  = "red") + 
        geom_hline(yintercept = controlMean + controlSD, color  = "red") +
        geom_hline(yintercept = controlMean - controlSD, color  = "red") + 
        labs(x=x,y=y) +
        coord_flip() +
        theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
      return(plot)
    }
  }
}

produceAllAllelePlot <- function(procedure_code, directory){
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  pdf(paste(date,"_ALLELEPLOT_",procedure_code,".pdf", sep = ""), 10,10)
  for(filename in files){
    center_name <- getCenterName(filename)
    center_data <- getCenterData(filename, procedure_code)
    if(center_name == 'HMGU'){
      phenotype_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
      for(phenotype in phenotype_list){
        if(is.element(phenotype, colnames(center_data))){
          title <- paste(center_name,"allele_symbol",";",phenotype,sep=" ")
          alleleplot_unsorted <- produceAllelePlot(center_data = center_data, phenotype = phenotype, title = title, isSorted = FALSE, center_name = center_name)
          print(alleleplot_unsorted)
          alleleplot_sorted <- produceAllelePlot(center_data = center_data, phenotype = phenotype, title = title, isSorted = TRUE, center_name = center_name)
          print(alleleplot_sorted)
        }  
      }
    }
    
  }
  dev.off()
}

produceAllZScorePlots <- function(procedure_code, directory, plotAllFile = T){
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  finished_phenotypes <- list()
  if(plotAllFile){
    pdf(paste(date,"_ALL_ZPlots_",procedure_code,".pdf", sep = ""), onefile = T, 10,50)
  }
  
  for(filename in files){
    center_data <- getCenterData(filename, procedure_code)
    # View(center_data[,c("Total.distance.travelled", "biological_sample_id")])
    phenotype_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    for(phenotype in phenotype_list){
      if(!is.element(phenotype, finished_phenotypes)){
        combined_df <- data.frame()
        #combined_df <- data.frame(biological_sample_id = integer(0), phenotyping_center = character(0), allele_symbol_zygosity = character(0), z_score = numeric(0))
        
        if(!plotAllFile){
          pdf(paste(date,"_",phenotype,"_ZPlots_",procedure_code,".pdf", sep = ""), 10,20)
        }
        for(current in files){
          current_data <- getCenterData(current, procedure_code)
          current_name <- getCenterName(current)
          current_list <-getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = current_data, procedure_code = procedure_code)
          if(is.element(phenotype, current_list) & is.element(phenotype,colnames(current_data))){
            
            
            # produce graph for single center
            if(!plotAllFile){
              title <- paste(current_name,"allele_symbol",";",phenotype,("z-scores"),sep=" ")
              alleleplot_sorted <- produceAllelePlot(center_data = current_data, phenotype = phenotype, title = title, isSorted = TRUE, center_name = current_name,z_score_plots = TRUE)
              print(alleleplot_sorted)
            }
            
            
            
            ######CODE TO MODIFY/MERGE DATA FRAMES TO CREATE FINAL COMBINED ALLELE PLOT#################
            
            # replace empty cells in allele column with controls
            modified_current <- current_data
            modified_current$allele_symbol <- as.character(modified_current$allele_symbol)
            modified_current$strain_name <- as.character(modified_current$strain_name)
            modified_current$phenotyping_center <- as.character(modified_current$phenotyping_center)
            # replacing blanks w/ "control"
            modified_current$allele_symbol[(modified_current$allele_symbol == "" | is.na(modified_current$allele_symbol)) & modified_current$biological_sample_group == 'control'] <- modified_current$strain_name[(modified_current$allele_symbol == "" | is.na(modified_current$allele_symbol)) & modified_current$biological_sample_group == 'control']
            modified_current$allele_symbol_zygosity <- paste(modified_current$phenotyping_center,modified_current$allele_symbol,modified_current$zygosity,sep=" ")
            
            # get z-scores
            modified_current_controls <- modified_current[modified_current$biological_sample_group == 'control',]
            controlMean <- mean(modified_current_controls[[phenotype]], na.rm = T)
            controlSD <- sd(modified_current_controls[[phenotype]], na.rm = T)
            modified_current$z_score <- (modified_current[[phenotype]]-controlMean) / controlSD
            
            
            # subset modified dataframe
            modified_current <- modified_current[,c("biological_sample_id", "phenotyping_center", "allele_symbol_zygosity", "z_score")]
            
            # add dataframe to combined dataframe
            combined_df <- rbind(combined_df, modified_current)
            
          }
        }
        x <- "allele_symbol_zygosity"
        y <- "z_score"
        # View(modified_current)
        # View(combined_df)
        # plot <- ggplot() +
        #   geom_pointrange(aes(x = reorder(combined_df$allele_symbol_zygosity, combined_df$z_score, FUN = mean,na.rm=TRUE),y = combined_df$z_score, na.rm = TRUE),data= combined_df,colour = '#0000ff',fill = '#ffffff',size = 0.4, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
        #   geom_point(aes(x = combined_df$allele_symbol_zygosity,y = combined_df$z_score, na.rm = TRUE),data=combined_df,shape = 20,colour = '#999999',size = 1.0, alpha = 0.3) +
        #   theme_bw(base_size = 4.0) + labs(title = title) +
        #   geom_hline(yintercept = -1, color  = "red") + 
        #   geom_hline(yintercept = 0, color  = "red") +
        #   geom_hline(yintercept = 1, color  = "red") +
        #   labs(x=x,y=y) +
        #   coord_flip() +
        #   theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
        if(plotAllFile){
          plot <- ggplot() +
            geom_pointrange(aes(x = reorder(combined_df$allele_symbol_zygosity, combined_df$z_score, FUN = mean,na.rm=TRUE),y = combined_df$z_score, color = combined_df$phenotyping_center, group = combined_df$phenotyping_center, na.rm = TRUE),data= combined_df,size = 0.4, shape = 20 ,fun.data = mean_sdl,fun.args = list(mult = 1),stat = 'summary') + 
            scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + theme_bw(base_size = 4.0) + labs(title = paste(phenotype," z-scores")) +
            geom_hline(yintercept = -1, color  = "red") + 
            geom_hline(yintercept = 0, color  = "red") +
            geom_hline(yintercept = 1, color  = "red") +
            labs(x=x,y=y,color="phenotyping centers") +
            coord_flip() +
            theme(axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
          print(plot)
          print(paste("COMPLETED FULL ZPLOT: ", phenotype,sep=""))
          
        }
        
        if(!plotAllFile){
          dev.off()
        }
        finished_phenotypes <- c(finished_phenotypes, phenotype)
      }
      print(paste("COMPLETED INDIVIDUAL ZPLOT: ", phenotype, sep=""))
    }
  }
  if(plotAllFile){
    dev.off()
  }
}

produceTimeDaySeries <- function(center_data, phenotype, timeplot, title = "", moving_average = 1, isLine = TRUE){
  x <- timeplot
  y <- phenotype
  center_data <- center_data[complete.cases(center_data[[x]]),]
  center_data <- center_data[complete.cases(center_data[[y]]),]
  center_mean <- mean(center_data[[y]], na.rm = TRUE)
  center_std <- sd(center_data[[y]], na.rm = TRUE)
  center_data[[x]] <- parse_date_time(center_data[[x]],orders = c("H:M:S"))
  # center_zoo <- zoo(center_data[[y]], center_data[[x]])
  # m.av <- rollmean(center_zoo, moving_average, fill=list(NA,NULL,NA))
  # center_data$m_av <- coredata(m.av)
  if(isLine){
    plot <- ggplot() +
      geom_point(aes(x = center_data[[x]],y = center_data[[y]],colour = sex),data=center_data,shape = 16,size = 0.5) +
      geom_line(aes(x = center_data[[x]],y = center_data[[y]]),data=center_data,size = 0.5,fun.data = mean_sdl,stat = 'summary') +
      scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + labs(x=x, y=y) + ggtitle(title) +
      stat_smooth(data=center_data, geom="smooth", position="") +
      geom_hline(yintercept = center_mean, color = "grey", size = 1) + geom_hline(yintercept = center_mean + center_std, color = "grey", size = 1) + geom_hline(yintercept = center_mean - center_std, color = "grey", size = 1) +
      theme(axis.text=element_text(size=6),axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
  }else{
    plot <- ggplot() +
      geom_point(aes(x = center_data[[x]],y = center_data[[y]],colour = sex),data=center_data,shape = 16,size = 0.5) +
      scale_colour_brewer(guide = guide_legend(),palette = 'Set1') + labs(x=x, y=y) + ggtitle(title) +
      stat_smooth(data=center_data, geom="smooth", position="") +
      geom_hline(yintercept = center_mean, color = "grey", size = 1) + geom_hline(yintercept = center_mean + center_std, color = "grey", size = 1) + geom_hline(yintercept = center_mean - center_std, color = "grey", size = 1) +
      theme(axis.text=element_text(size=6),axis.title=element_text(size=8,face="bold"), plot.title = element_text(size = 9,face="bold"), legend.text = element_text(size=7), legend.title=element_text(size=8,face="bold"))
    
  }
  
  return(plot)
}

produceAllTimeDaySeries <- function(procedure_code, directory, isControl = TRUE, moving_average = 1,isLine=TRUE){
  procedure_table <- getProcedureTable(procedure_code, directory)
  IMPC_dictionary <- getIMPCDictionary(directory)
  files <- getWideCSVFiles(procedure_code, directory)
  output_directory <- paste(directory, procedure_code, sep="")
  setwd(output_directory)
  date <- format(Sys.time(),"%Y%m%d_%HH%MM_")
  if(isControl){
    control_or_not <- "Controls"
  }else{
    control_or_not <- "All"
  }
  if(!isLine){
    line_or_not <- "ASPOINTS"
    pdf_name <- paste(date, control_or_not,"_ASPOINTS_TIMEDAYSERIES_",procedure_code,".pdf", sep = "")
  }else{
    line_or_not <- "ASLINE"
    pdf_name <- paste(date, control_or_not,"_TIMEDAYSERIES_",procedure_code,".pdf", sep = "")
  }
  pdf(pdf_name, 14,4)
  for(filename in files){
    center_data <- getCenterData(filename, procedure_code)
    center_name <- getCenterName(filename)
    phenotype_list <- getPhenotypes(procedure_table = procedure_table, IMPC_dictionary = IMPC_dictionary, center_data = center_data, procedure_code = procedure_code)
    timeplot_list <- list('Start_hms_new')
    if(isControl){
      center_data <- center_data[center_data$biological_sample_group == 'control',]
    }
    for(phenotype in phenotype_list){
      if(is.element(phenotype, colnames(center_data))){
        for(timeplot in timeplot_list){
          if(is.element(timeplot, colnames(center_data))){
            title <- paste(center_name,phenotype,";", timeplot,"\n","(",control_or_not,")",sep=" ")
            time_series_plot <- produceTimeDaySeries(center_data = center_data, phenotype = phenotype, timeplot = timeplot, title = title, moving_average = moving_average,isLine = isLine)
            tryCatch(print(time_series_plot), error = function(e){print(paste(center_name,"threw an error",sep=" "))})
          }
        }
      }
    }
    print(paste("COMPLETED TIMEDAYSERIES: ",control_or_not," ",line_or_not," ",center_name,sep=""))
  }
  dev.off()
}

CV <- function(mean, sd){
  return((sd/mean)*100)
}

produceCVPlot <- function(procedure_code, directory, phenotype, title = "", isSorted = T){
  x <- "allele_symbol_zygosity"
  y <- phenotype
  
  files <- getWideCSVFiles(procedure_code, directory)
  index = 0
  for(file in files && index == 0){
    
    index = 1
    
    center_data <- getCenterData(file, procedure_code)
    
    center_data$allele_symbol <- as.character(center_data$allele_symbol)
    center_data$strain_name <- as.character(center_data$strain_name)
    
    # set blanks to 'control'
    center_data <- center_data[which((center_data$allele_symbol == "" | is.na(center_data$allele_symbol)) & center_data$biological_sample_group == 'control'), ]
    center_data$allele_symbol <- center_data$strain_name
    View(center_data)    
  }
  
}

produceAllPlots <- function(procedure_code, directory, pipeline = "IMPC"){
  
  start_time <- Sys.time()
  print(paste("START TIME: ", start_time, sep=""))
  # BOXPLOT
  print(paste("STARTED BOX PLOTS: ", procedure_code, sep=""))
  produceAllBoxPlot(procedure_code, directory, isControl = F, pipeline)
  produceAllBoxPlot(procedure_code, directory, isControl = T, pipeline)
  print(paste("FINISHED BOX PLOTS: ", procedure_code, sep=""))
  
  # ALLELE
  print(paste("STARTED ZSCORE PLOTS: ", procedure_code, sep=""))
  produceAllZScorePlots(procedure_code, directory, plotAllFile = F)
  produceAllZScorePlots(procedure_code, directory, plotAllFile = T)
  print(paste("FINISHED ZSCORE PLOTS: ", procedure_code, sep=""))
  
  # TIMESERIES
  print(paste("STARTED TIMESERIES PLOTS: ", procedure_code, sep=""))
  produceAllTimeSeries(procedure_code, directory, isControl = F,moving_average = 1)
  produceAllTimeSeries(procedure_code, directory, isControl = T,moving_average = 1)
  print(paste("FINISHED TIMESERIES PLOTS: ", procedure_code, sep=""))
  
  # MULTITIMESERIES
  print(paste("STARTED MULTITIMESERIES PLOTS: ", procedure_code, sep=""))
  produceAllMultiTimeSeries(procedure_code,directory,isControl=T,isLine=F, pipeline=pipeline)
  produceAllMultiTimeSeries(procedure_code,directory,isControl=F,isLine=F, pipeline=pipeline)
  print(paste("FINISHED MULTITIMESERIES PLOTS: ", procedure_code, sep=""))
  
  # TIMEDAYSERIES
  # print(paste("STARTED MULTITIMESERIES PLOTS: ", procedure_code, sep=""))
  # produceAllTimeDaySeries(procedure_code,directory,isControl=F,moving_average=1,isLine=T)
  # produceAllTimeDaySeries(procedure_code,directory,isControl=T,moving_average=1,isLine=T)
  # produceAllTimeDaySeries(procedure_code,directory,isControl=F,moving_average=1,isLine=F)
  # produceAllTimeDaySeries(procedure_code,directory,isControl=T,moving_average=1,isLine=F)
  # print(paste("FINISHED MULTITIMESERIES PLOTS: ", procedure_code, sep=""))
  
  end_time <- Sys.time()
  print(paste("END TIME: ", end_time, sep=""))
  
  elapsed_time <- end_time - start_time
  print(paste("ELAPSED TIME: ", elapsed_time, sep=""))
  
  
}

#############MAIN####################
directory <- "C:\\Users\\svprice\\Dropbox\\IMPC_Data\\"
procedure_code <- "ECG"
suppressWarnings(produceAllPlots(procedure_code, directory,pipeline="IMPC"))

procedure_code <- "HEM"
pipeline <- "IMPC"
produceAllBoxPlot(procedure_code, directory, isControl = F, pipeline)
produceAllBoxPlot(procedure_code, directory, isControl = T, pipeline)
produceAllMultiTimeSeries(procedure_code,directory,isControl=T,isLine=F, pipeline=pipeline)
produceAllMultiTimeSeries(procedure_code,directory,isControl=F,isLine=F, pipeline=pipeline)

