#'
#' Called within the plotManhattan() function that generates a PheWAS Manhattan plot with phenotypes on the x-axis and p-values on a -log10 scale on the y-axis. This function enables some additional plotting options compared to the plotManhattan() function.
#'
#' @param d The data frame passed into the plotManhattan() function, which is typically the data 	frame of PheWAS results returned by the phewas_ext() function.
#' @param suggestive.line The p-value at which to draw a blue line to show a significance threshold. The default value is 0.05. An NA value results in no line.
#' @param significant.line The p-value at which to draw a red line to show an adjusted significance threshold. The default applies a Bonferroni correction. An NA value results in no line.
#' @param size.x.labels The size of the x-axis labels (default = 9).
#' @param size.y.labels The size of the y-axis labels (default = 9).
#' @param point.size The size of the points on the plot (default = 3).
#' @param annotate.phenotype Determines whether phenotype descriptions (corresponding to phecodes with p < annotate.level or phecodes specified in annotate.list) are annotated on the Manhattan plot. The default is TRUE.
#' @param annotate.angle The angle used for the annotation text (default = 0).
#' @param annotate.size The size of the annotation text (default = 3).
#' @param annotate.level The p-value threshold for annotating phenotypes. The default value is that provided in the significant.line argument.
#' @param annotate.list A character vector of phecodes to be annotated on the plot regardless of significance.
#' @param lc.labels Determines if the labels should be displayed in lower case (default is FALSE).
#' @param y.axis.interval The interval used for the y-axis of the Manhattan plot (default value = 5).
#' @return phenotype-significance Manhattan plot to its parent function
#' @source https://github.com/monikagrabowska/PedsPheWAS/blob/main/R/phenotypePlot.R
#' @source https://github.com/PheWAS/PheWAS/blob/master/R/phenotypePlot.R
#' @export

phenotypePlot <-
  function(d, suggestive.line, significant.line,
           size.x.labels=9, size.y.labels=9,
           point.size = 3,
           annotate.phenotype=T,
           annotate.angle=0, annotate.size=3, annotate.level,
           annotate.list,
           lc.labels=F,
           y.axis.interval=y.axis.interval) {

    d=merge(d,ProcWAS::ccsr_phewas_plot_annotations_cleaned,by.x="phenotype",by.y="proc_code")

    d=d[!is.na(d$groupnum),]

    d$size=3

    #Remove lc.labels flag if no annotations
    if(!annotate.phenotype) lc.labels=F

    #Sort by the phenotype
    d=d[order(d$phenotype),]

    #Set the maximum x value to fit all phenotypes
    max.x = length(unique(d$phenotype)) + 1

    #Create the list of phenotypes, finding the best values for each phenotype
    phenotypes=aggregate(value ~ phenotype + groupnum, d,FUN=max)
    #Remove the least significant phenotypes; only has an effect if max.x was specified.
    phenotypes=phenotypes[order(phenotypes$value, decreasing=T),][1:min(nrow(phenotypes),max.x),]

    phenotypes=phenotypes[order(phenotypes$groupnum,phenotypes$phenotype),]

    phenotypes$seq = 1:nrow(phenotypes)

    #Limit to phenotype and seq, as they are the only relevant columns
    #Include value as min.value for annotation purposes
    phenotypes=phenotypes[,c("phenotype","seq","value", "groupnum")]
    names(phenotypes)[3]="min.value"

    #Add sequence information
    d=inner_join(phenotypes,d,by=c("phenotype", "groupnum"))
    d=d[order(d$seq),]

    #Define the max y axis value if not provided
    max.y=ceiling(max(d$value))+1

    labels= summarize(group_by(d, groupnum), tick=mean(unique(seq)),label=as.character(clinical_domain[1]))
    labels=labels[order(labels$tick),]

    color.palette = unique(d[order(d$seq),]$colors)
    names(color.palette)=color.palette

    y.axis.label=expression(-log[10](italic(p)))
    x.axis.label=""

    #Generate the inital plot
    plot=ggplot(d,ylab=y.axis.label,xlab=x.axis.label) + theme(axis.title = element_text(face="bold",size=15))

    #Set the Y scale and labels
    plot=plot+scale_y_continuous(y.axis.label, limits=c(0,max.y), breaks=seq(0,max.y,y.axis.interval), expand=c(0,0))

    #Include lines for significance thresholds
    if (suggestive.line<=max.y && !missing(suggestive.line) && !is.na(suggestive.line)) {
      plot=plot+geom_hline(yintercept=suggestive.line,colour="blue", alpha=I(1/3),size=0.3)
    }
    if (significant.line<=max.y && !missing(significant.line) && !is.na(significant.line)) {
      plot=plot+geom_hline(yintercept=significant.line,colour="red",alpha=I(1/3),size=0.3)
    }

    plot=plot+aes(seq,value,size=size,colour=colors)
    plot=plot+scale_size(range=c(point.size,point.size),guide="none")
    #Add points
    plot=plot+geom_point()

    #Color as defined
    plot = plot + scale_colour_manual(values= color.palette, guide="none")

    #Label the X axis with the groups
    plot=plot+scale_x_continuous(name=x.axis.label, limits=c(1,max.x), breaks=labels$tick, labels=labels$label, expand=c(0,0))

    #Set the default theme
    plot=plot+theme(
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      axis.text.x=element_text(size=size.x.labels, colour="black", angle=-50, hjust=0, vjust=1),
      axis.text.y=element_text(size=size.y.labels, colour="black"),
      axis.line =element_line(colour="black"),
      axis.ticks=element_line(colour="black")
    )

    #Hide the legend by default
    plot = plot+theme(legend.position = "none")

    #Add OR information
    plot=plot+aes(shape = factor(direction), fill=colors) +
      scale_shape_manual(values=c(25,24)) +
      scale_fill_manual(values= color.palette, guide="none")

    #If annotation present, start the definitions, otherwise skip it
    if (annotate.phenotype) 	{
      d$annotate=F
      #If provided with a list of phenotypes to annotate, select those.
      if(!missing(annotate.list)) d[d$phenotype %in% annotate.list, ]$annotate=T
      #Include those above the given threshold
      if((!missing(annotate.level) && !is.na(annotate.level)) && (sum(d$value>=annotate.level)>0)) d[d$value>=annotate.level, ]$annotate=T

      #Cap annotation length
      d$proc_code_description = substr(d$proc_code_description,1,100)
      #Add leading space
      d$proc_code_description = paste0("  ",d$proc_code_description)

      #If lower case labels are requested, lower case them.
      if(lc.labels) d$proc_code_description=tolower(d$proc_code_description)

      if(sum(d$annotate)==0) {
        message("No points met the annotation criteria.")
      }
      else {
        plot = plot + ggrepel::geom_text_repel(aes(label=proc_code_description),colour="black",data=d[d$annotate,],size=annotate.size,angle=annotate.angle,force = 5)
      }
    }
    print(plot)
  }