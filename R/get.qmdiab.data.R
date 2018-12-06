#' Calculate eigenmetabolites as module representatives
#'
#'
#' @param url the link to where QMDiab data is stored.
#' @references
#' \insertRef{donetwork2015}{MoDentify}
#' \insertRef{mookkanamori2013}{MoDentify}
#' \insertRef{suhreconnecting2017}{MoDentify}
#' \insertRef{yousri2015}{MoDentify}
#' @export
#' @importFrom readxl read_excel
#' @return list of three, data: \code{\link[data.table]{data.table}} of 
#' preprocssed metabolite concentrations, annos: 
#' \code{\link[data.table]{data.table}} of annotations for each metabolite, 
#' phenos: \code{\link[data.table]{data.table}} of age, sex, BMI, and T2D information
get.qmdiab.data <- function(filename){
  
  ###- Read data
  
  #- metabolites:
  pdat <- read_excel(filename, sheet="plasma") # plasma data
  udat <- read_excel(filename, sheet="urine") # urine data
  sdat <- read_excel(filename, sheet="saliva") # saliva data
  
  #- annotations:
  panno <- read_excel(filename, sheet="plasma annotations")  # plasma annotations
  uanno <- read_excel(filename, sheet="urine annotations")  # urine annotations
  sanno <- read_excel(filename, sheet="saliva annotations")  # saliva annotations
  
  
  ###- Get overlap of all three fluids
  
  inter <- intersect(intersect(pdat$`QMDiab-ID`, udat$`QMDiab-ID`), sdat$`QMDiab-ID`)
  
  pind <- match(inter, pdat$`QMDiab-ID`)
  uind <- match(inter, udat$`QMDiab-ID`)
  sind <- match(inter, sdat$`QMDiab-ID`)
  
  `QMDiab-ID` <-  t(pdat[pind,1])
  
  # phenotypes
  phenos <- pdat[pind,c("AGE", "GENDER", "BMI", "T2D")]
  rownames(phenos) <- `QMDiab-ID`
  
  # data
  pdat <- pdat[pind,2:(dim(pdat)[2]-5)]
  udat <- udat[uind,2:(dim(udat)[2]-5)]
  sdat <- sdat[sind,2:(dim(sdat)[2]-5)]
  
  data <- cbind(pdat,udat,sdat)
  rownames(data) <- `QMDiab-ID`
  
  
  ###- Assign fluid-specific annotations
  
  annos <- rbind(panno,uanno,sanno)
  annos$Fluid <- c(rep("P", length(panno$BIOCHEMICAL)),rep("U", length(uanno$BIOCHEMICAL)),rep("S", length(sanno$BIOCHEMICAL)))
  
  annos$name <- paste0(annos$Fluid, "::", annos$BIOCHEMICAL) # this column is needed by MoDentify
  annos$SUPER_PATHWAY <- paste0(annos$Fluid, "::", annos$SUPER_PATHWAY)
  annos$SUB_PATHWAY <- paste0(annos$Fluid, "::", annos$SUB_PATHWAY)
  
  colnames(data) <- annos$name

  ###- scoring
  data <- scale(data)
  
  ###- return
  return(list(annos=annos, data=data, phenos=phenos))
}