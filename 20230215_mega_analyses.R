########################################################################################################################################
##  ENIGMA - Panic Disorder - Linear Mixed Model Mega-analyses                                                                        ##
##  by Willem Benjamin Bruin                                                                                                          ##
##                                                                                                                                    ##
##  this script is adapted from the work of Max Laansma (ENIGMA-Parkinson's Disorder Working Group)                                   ##
##  and existing ENIGMA scripts for multiple linear models see: https://github.com/ENIGMA-git/ENIGMADiseaseWorkingGroupStats          ##
##                                                                                                                                    ##
##  Usage:                                                                                                                            ##
##  source('path/to/scripts/20230215_mega_analyses.R', encoding='ISO-8859-1')                                                          ##
##                                                                                                                                    ##
########################################################################################################################################
##                                                                                                                                    ##
##                                                  #### SECTIONS OUTLINE: ####                                                       ##
##                                                                                                                                    ##
## SECTION 0: Set up                                                                                                                  ##
##            set working directory and download required R packages                                                                  ##
##            functions needed for effect statistics (d and r), SE for each effect statistic, and CI calculations                     ##
##                                                                                                                                    ##
## SECTION 1: Group comparisons - Primary analyses of PD vs HC and clinical subgroup analyses                                         ##
## SECTION 2: Sub-analyses: interactions on full dataset                                                                              ##
## SECTION 3: Sub-analyses: lme for continuous variable models using Pearson's r        
## SECTION 4: Post-Hoc analyses: between patient comparisons   
##                                                                                                                                    ##
##            .1 Create filtered data sets/Set up variables for each model                                                            ##
##            .2 Setup analysis:                                                                                                      ##
##            .3 Loop through all FreeSurfer ROIs (volume) and run each LME model separately                                                 ##      
##            .4 Collect effect size data                                                                                             ##      
##            .5 Combine all output variables into table                                                                              ##
##                                                                                                                                    ##
########################################################################################################################################

options(error=recover)

###### SECTION 0: Set up ######

# Threshold for minimum number of PD patients required per site. 
N_min = 5 
# Threshold for minimum number of sites required. 
N_min_sites = 4

# Preferred Multiple Comparison Procedure (MCP) settings, 'fdr' uses Benjamini & Hochberg procedure
MCP_method = 'fdr'#
MCP_alpha = 0.05 # 

### Set working directory and create result directory for outputs 
mainDir <- "/X/X/X/ENIGMA_PD/data/"
mainResultsDir <- file.path(mainDir, paste("LH_mega_analysis_results_final",  MCP_method, MCP_alpha, sep = "_", 
                                           collapse = NULL))
dir.create(mainResultsDir, showWarnings = FALSE)
setwd(mainResultsDir)
cat("Working directory is set as: ", getwd(), "\nReading source file.\n")
sink("mega_analysis.log")

# Download required R packages 
# Comment any packages that have already been locally installed
library("data.table")
library("ppcor")
library("lme4")
library("nlme")
library("car")
library("dplyr")
library("emmeans")
library("xlsx")
library("ggplot2")
library("effects")
library("ggridges")

#### Read csv file containing the FreeSurfer (FS) data and covariates
full_data <- read.csv("/path/to/data/ENIGMA_Anxiety_DATA.csv", 
                      na.strings=c("","NA"), sep =',', dec=".")

# Check whether all subjects have required columns (Age, Sex & Dx [Diagnosis])
{missing_data<-full_data[!complete.cases(full_data[ , c('Age', 'Sex', 'Dx')]), ]
if (nrow(missing_data) > 0){
  stop("Data in provided CSV is incomplete. Ensure all partipants have Age, Sex and Dx!")
}}

# Import functions required for LME analyses
cat("Importing functions needed for effect size, se and CI calculations \n")
source("/path/to/scripts/20230215_stat_functions.R")
cat("Importing functions needed for exclusion of samples with insufficient size \n")
source("/path/to/scripts/20230215_filter_functions.R")

# Specify column ID depicting sample Site ID
multi_site_col_ID <- which(colnames(full_data) == 'MultiSiteID')

# Specify column indices for different FS modalities separately
subcort_volume_idx <- which(colnames(full_data) == 'LLatVent'):which(colnames(full_data) == 'Raccumb')
cort_surfarea_idx <- which(colnames(full_data) == 'L_bankssts_surfavg'):which(colnames(full_data) == 'R_insula_surfavg')
cort_thickness_idx <- which(colnames(full_data) == 'L_bankssts_thickavg'):which(colnames(full_data) == 'R_insula_thickavg')
total_thickness_idx <- which(colnames(full_data) == 'LThickness'):which(colnames(full_data) == 'RThickness')
total_surfarea_idx <- which(colnames(full_data) == 'LSurfArea'):which(colnames(full_data) == 'RSurfArea')
FS_indices<-list(subcort_volume_idx, c(cort_surfarea_idx, total_surfarea_idx), c(cort_thickness_idx, total_thickness_idx))
names(FS_indices)<-c("subcortical_volume", "cortical_surface_area", "cortical_thickness")

# Specify columns that depict if subject has data for particular FS modality
has_FS_cols = list('has_LandRvolumes.csv', 'has_CorticalMeasuresENIGMA_SurfAvg.csv', 
                   'has_CorticalMeasuresENIGMA_ThickAvg.csv', 'has_tot_cort_thick', 'has_tot_cort_surf')
names(has_FS_cols)<-c("subcortical_volume", "cortical_surface_area", "cortical_thickness", 
                     "total_cortical_thickness", "total_cortical_surface_area")

### Define age groups to perform analyses separately (optional)
data_age_grouped=list(full_data)
names(data_age_grouped)<-c("combined")
# adolescent_PDHC <- data.frame(filter(full_data, (AgeGroup=='adolescent')))
# adult_PDHC <- data.frame(filter(full_data, (AgeGroup=='adult')))
# data_age_grouped=list(full_data, adolescent_PDHC, adult_PDHC)
# names(data_age_grouped)<-c("combined", adolescent", "adult")

# Ensure FreeSurfer columns contain no strings (other than NaN)
for (f in 1:length(FS_indices)){
  feature_col_indices <- FS_indices[[f]]
  for (i in feature_col_indices){
    FreeSurfer_ROI = full_data[, i]
    if (is.character(FreeSurfer_ROI)){
      stop("String found in ",names(FS_indices[f]),", col ID:",i,", ",colnames(full_data)[i],". Stopping analysis!")
    }
  }
} 

# Loop through age groups (i.e. adolescents, adult or combined samples)
for (agegroup_i in 1:length(data_age_grouped)){
  
  agegroup_data <- data_age_grouped[[agegroup_i]]
  agegroup_label <- names(data_age_grouped[agegroup_i])
  agegroupResultsDir <- file.path(mainResultsDir, agegroup_label)
  dir.create(agegroupResultsDir, showWarnings = FALSE)
  setwd(agegroupResultsDir)
  plot.age_dist(agegroup_data, N_min) # Plot age distributions per site
  
  # Loop through different FS modalities and perform analyses separately
  for (f in 1:length(FS_indices)){
    feature_col_indices <- FS_indices[[f]]
    feature_label <- names(FS_indices[f])
    subResultsDir <- file.path(agegroupResultsDir, feature_label)
    cat("Running mega-analyses for", agegroup_label, "samples and ", feature_label, "features \n\n")
    dir.create(subResultsDir, showWarnings = FALSE)
    setwd(subResultsDir)
    
    ###### SECTION 1: Group comparisons - Primary analyses of PD vs HC and clinical subgroup analyses ######
    #### 1.1 Create filtered data sets                                                                ######
    #### Only use samples that have specified FS modality (using "has_" columns) and are in age_group ######
    agegroup_feature_data <- agegroup_data[which(agegroup_data[has_FS_cols[[f]]] == "True"), ]
    N_sites <- length(unique(agegroup_feature_data$MultiSiteID))
    N_ROIs <- length(feature_col_indices)

    # Define group comparisons for main case-control analyses
    m.PDHC<-data.frame(filter(agegroup_feature_data, (PD>0) | (Dx==0 )))
    m1.Cur<-data.frame(filter(agegroup_feature_data, (PD==2) | (Dx==0 )))
    m2a.ANX<-data.frame(filter(agegroup_feature_data, (PD>0 & ComANXD==1) | (Dx==0 )))
    m2b.woANX<-data.frame(filter(agegroup_feature_data, (PD>0 & ComANXD==0) | (Dx==0 )))
    m3a.MDD<-data.frame(filter(agegroup_feature_data, (PD>0 & MDD>0) | (Dx==0 )))
    m3b.woMDD<-data.frame(filter(agegroup_feature_data, (PD>0 & MDD==0) | (Dx==0 )))
    m4a.Meds<-data.frame(filter(agegroup_feature_data, (PD>0 & Med==2) | (Dx==0 )))
    m4b.SSRI<-data.frame(filter(agegroup_feature_data, (PD>0 & SSRI_SNRI==1) | (Dx==0 )))
    m4c.woMeds<-data.frame(filter(agegroup_feature_data, (PD>0 & Med==1) | (Dx==0 )))
    m5a.EarlyOns<-data.frame(filter(agegroup_feature_data, (PD>0 & AgeO<=21) | (Dx==0 )))
    m5b.LateOns<-data.frame(filter(agegroup_feature_data, (PD>0 & AgeO>21) | (Dx==0 )))
    
    datasets<-list(m.PDHC,m1.Cur,m2a.ANX,m2b.woANX,m3a.MDD,m3b.woMDD,m4a.Meds,m4b.SSRI,m4c.woMeds,m5a.EarlyOns, m5b.LateOns)
    names(datasets)<-c("m.PDHC","m1.Cur","m2a.ANX","m2b.woANX","m3a.MDD","m3b.woMDD","m4a.Meds","m4b.SSRI","m4c.woMeds", "m5a.EarlyOns", "m5b.LateOns")
    cat("Beginning SECTION 1: Primary analysis and sub-analyses with filtered data sets\n\nModels in section 1 are:\n", names(datasets), "\n")
    
    ### Filter out sites with < N_min patients per (sub)-analysis across ALL regions within a given FreeSurfer modality ###
    dataList_1 = filter_datasets_section_1(datasets, N_sites, N_ROIs, N_min, N_min_sites)
    if (length(dataList_1) > 0){
      cat("Finished filtering datasets for:", feature_label, "\nFollowing models survived thresholding:", names(dataList_1), "\n\n")
    } else {
      cat("Finished filtering datasets for:", feature_label, "\nNO models survived thresholding! Skipping analyses... \n\n")
      next
    }
    
    #### 1.2 Setup analysis: effect of diagnostic group ####
    for(l in 1:length(dataList_1)){
      
      ## Run LME models for each (filtered) dataset separately ##
      cat("working on model: ", names(dataList_1)[l], "\n") 
      
      ### Initialize arrays to store LME fit results
      n.controls=rep(0,N_ROIs)
      n.PD=rep(0,N_ROIs)
      n.sites=rep(0,N_ROIs)
      p.value<-rep(0,N_ROIs)
      db.PDHC=rep(0,N_ROIs)
      se.db.PDHC=rep(0,N_ROIs)
      perc.diff.PDHC=rep(0,N_ROIs)
      L.ci.db.PDHC=rep(0,N_ROIs)
      U.ci.db.PDHC=rep(0,N_ROIs)
      est.dxPDHC=rep(NA,N_ROIs)
      se.beta.PDHC=rep(0,N_ROIs)

      #### 1.3 Loop through all FreeSurfer features in set and run each LME model separately ####
      sink(paste0(names(dataList_1)[l],"_ModelSummaryStats.txt"))
      
      for (i in seq_along(feature_col_indices)){
        
        # FreeSurfer column index of interest
        x = feature_col_indices[[i]]
        
        ### Remove sites that have less than N_min examples per class
        FreeSurfer_ROI_HC=data.frame(filter(dataList_1[[l]],(Dx==0)))[c(multi_site_col_ID, x)]
        FreeSurfer_ROI_PD=data.frame(filter(dataList_1[[l]],(Dx>0)))[c(multi_site_col_ID, x)]
        names(FreeSurfer_ROI_HC)[2]<-"feature_x"; names(FreeSurfer_ROI_PD)[2]<-"feature_x"
        non_na_HC<-setDT(FreeSurfer_ROI_HC)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
        non_na_PD<-setDT(FreeSurfer_ROI_PD)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
        
        # Check if remaining number of sites is >=N_min_sites, and at least one site has both classes included, else skip model
        filtered_dataList = copy(dataList_1[[l]])
        sites_to_exclude_HC<-data.frame(filter(non_na_HC,(non_na<N_min)))$MultiSiteID
        sites_to_exclude_PD<-data.frame(filter(non_na_PD,(non_na<N_min)))$MultiSiteID
        if (length((sites_to_exclude_HC)) > 0){
          filtered_dataList=subset(filtered_dataList, !(MultiSiteID%in%sites_to_exclude_HC & Dx==0))
        }
        if (length((sites_to_exclude_PD)) > 0){
          filtered_dataList=subset(filtered_dataList, !(MultiSiteID%in%sites_to_exclude_PD & Dx>0))
        }        
        filtered_tbl<-table(filtered_dataList$MultiSiteID, filtered_dataList$Dx)
        if (dim(filtered_tbl)[1] < N_min_sites | dim(filtered_tbl)[2] < 2 | !any(apply(filtered_tbl, 1, function(x) all(x>0)))){
          next
        }
        
        # FreeSurfer label of interest
        cat(names(filtered_dataList)[x], "\n")
        FreeSurfer_ROI=filtered_dataList[,x]

        # Define base model variables
        Dx<-as.factor(filtered_dataList$Dx)
        Sex<-as.factor(filtered_dataList$Sex)
        Age<-as.numeric(filtered_dataList$Age)
        AgeC<-Age-mean(Age)
        AgeC2<-AgeC^2
        ICV<-as.numeric(filtered_dataList$ICV)
        MultiSiteID<-as.factor(filtered_dataList$MultiSiteID)
        
        ctrl <- lmeControl(opt='optim');
        # LME models for subcortical volumes and cortical surface area include TIV
        if (feature_label %in% list("subcortical_volume", "cortical_surface_area")){
          lme1_tmp_data <- data.frame(FreeSurfer_ROI, Dx, Sex, AgeC, AgeC2, ICV, MultiSiteID)
          lme1=lme(FreeSurfer_ROI ~ Dx + Sex + AgeC + AgeC2 + ICV, 
                   random= ~1 | MultiSiteID, na.action="na.exclude", control=ctrl, data=lme1_tmp_data)
        } else {
          lme1_tmp_data <- data.frame(FreeSurfer_ROI, Dx, Sex, AgeC, AgeC2, MultiSiteID)
          lme1=lme(FreeSurfer_ROI ~ Dx + Sex + AgeC + AgeC2, 
                   random= ~1 | MultiSiteID, na.action="na.exclude", control=ctrl, data=lme1_tmp_data)
        }
        fit1=summary(lme1); lsm1=summary(lsmeans(lme1, "Dx"))
        lme1_data=retrieve.LME.sample(lme1)
        cat("Model summary statistics are: \n"); print(fit1)
        coefs1 <- data.frame(coefficients(fit1))
        p.value[[i]]=coefs1[which(rownames(coefs1) == "Dx1"), ]$p.value
        n.controls[[i]]=length(which(lme1_data$Dx==0)); n.PD[[i]]=length(which(lme1_data$Dx==1))
        n.obs=fit1[["dims"]][["N"]]; n.groups=fit1[["dims"]][["ngrps"]][["MultiSiteID"]]
        n.sites[[i]]=n.groups; 
        param=length(fit1$coefficients$fixed) # (k=num parameters including intercept)
        perc.diff.PDHC[[i]]=lsm(lsm1[2,2],lsm1[1,2]) # First PD, then HC
        varcorr=VarCorr(fit1) # Extract intercept (random effect-stat) variance and residual variance 
        R.btwn=as.numeric(varcorr[which(rownames(varcorr) == "(Intercept)"), which(colnames(varcorr) == "Variance")])
        R.wthn=as.numeric(varcorr[which(rownames(varcorr) == "Residual"), which(colnames(varcorr) == "Variance")])
        R.calc=as.numeric(mixeff.R(R.btwn,R.wthn))
        # Extract T value for Diagnosis
        tval.Dx=fit1$tTable[which(rownames(fit1$tTable) == "Dx1"), which(colnames(fit1$tTable) == "t-value")]
        # Calculate effect size, standard error and 95 CI
        db.PDHC[[i]]=mixeff.d(tval.Dx, n.groups, n.obs, n.controls[[i]], n.PD[[i]], R.calc, param)
        se.db.PDHC[[i]]=se.db(db.PDHC[[i]],n.controls[[i]],n.PD[[i]])
        bound.db.PDHC=CI1(db.PDHC[[i]],se.db.PDHC[[i]])
        L.ci.db.PDHC[[i]]=bound.db.PDHC[1]; U.ci.db.PDHC[[i]]=bound.db.PDHC[2]
        est.dxPDHC[[i]]=coefs1[which(rownames(coefs1) == "Dx1"), ]$Value
        se.beta.PDHC[[i]]=coefs1[which(rownames(coefs1) == "Dx1"), ]$Std.Error
      }  
      sink()
      
      ## 1.4.1 calculate corrected p-values for overall model ##           
      tmp.MCP.p.overall=p.adjust(p.value,method=MCP_method) # Default using FDR with alpha p_corrected < 0.05
      MCP.p.overall=rep(0,N_ROIs)
      for (w in 1:N_ROIs){MCP.p.overall[w]=tmp.MCP.p.overall[w]}
      significant_MCP<-as.character(MCP.p.overall<MCP_alpha)
      #### 1.5 Combine all output variables into table ####
      outmatoverall=cbind(significant_MCP,db.PDHC,se.db.PDHC,perc.diff.PDHC,L.ci.db.PDHC,U.ci.db.PDHC,
                          MCP.p.overall,p.value,n.PD,n.controls,n.sites,est.dxPDHC,se.beta.PDHC)
      rownames(outmatoverall)=names(dataList_1[[l]])[feature_col_indices]
      cat("Generating final output for model: ", names(dataList_1)[l],"\n")
      ## 1.5.1 Save output
      write.xlsx(outmatoverall, file = paste0("megaResults.xlsx"),sheetName=(names(dataList_1)[l]), append=TRUE)
      
    }
    cat("\n")
    
    ###### SECTION 2: Sub-analyses: interactions (models: 6a,6b,6c) on full dataset ######
    cat("\nSection 1 completed. \n\nBeginning SECTION 2: Sub-analyses for interactions 6a,6b,6c \n")
    
    #### 2.1 Create datasets and set up variables for each model                                       ####
    #### No need to filter these models again because same data from main model (m.PDHC) is used here. ####
    m6a.DxAgeC<-dataList_1[[which(names(dataList_1) == "m.PDHC")]]
    m6b.DxSex<-dataList_1[[which(names(dataList_1) == "m.PDHC")]]
    m6c.DxAgeC2<-dataList_1[[which(names(dataList_1) == "m.PDHC")]]
    dataList_2<-list(m6a.DxAgeC,m6b.DxSex,m6c.DxAgeC2)
    names(dataList_2)<-c("m6a.DxAgeC","m6b.DxSex","m6c.DxAgeC2")
    sample_sizes_before_QC <- as.data.frame(read.xlsx("site_sample_sizes_BEFORE_filtering.xlsx", 
                                        sheetIndex=which(names(dataList_1) == "m.PDHC"), header=TRUE))
    sample_sizes_after_QC <- as.data.frame(read.xlsx("site_sample_sizes_AFTER_filtering.xlsx", 
                                       sheetIndex=which(names(dataList_1) == "m.PDHC"), header=TRUE))
    colnames(sample_sizes_before_QC)[1] = ""; colnames(sample_sizes_after_QC)[1] = ""
    for (d in 1:length(dataList_2)){
      write.xlsx(sample_sizes_before_QC, file = paste0("site_sample_sizes_BEFORE_filtering.xlsx"), row.names = FALSE, 
                 sheetName=names(dataList_2)[[d]], append=TRUE)
      write.xlsx(sample_sizes_after_QC, file = paste0("site_sample_sizes_AFTER_filtering.xlsx"), row.names = FALSE,
                 sheetName=names(dataList_2)[[d]], append=TRUE)
    }
  
    # Define model variables as strings to paste into lme (with or without ICV) ##
    if (feature_label %in% list("subcortical_volume", "cortical_surface_area")){
      lm.DxAgeC<-list("Dx","Sex","AgeC","AgeC2","ICV","Dx:AgeC")
      lm.DxSex<-list("Dx","Sex","AgeC","AgeC2","ICV","Dx:Sex")
      lm.DxAgeC2<-list("Dx","Sex","AgeC","AgeC2","ICV","Dx:AgeC2")
    } else {
      lm.DxAgeC<-list("Dx","Sex","AgeC","AgeC2","Dx:AgeC")
      lm.DxSex<-list("Dx","Sex","AgeC","AgeC2","Dx:Sex")
      lm.DxAgeC2<-list("Dx","Sex","AgeC","AgeC2","Dx:AgeC2")
    }
    lm.List1<-list(lm.DxAgeC,lm.DxSex,lm.DxAgeC2)
    names(lm.List1)<-c("lm.DxAgeC","lm.DxSex","lm.DxAgeC2") 
    
    #### 2.2 Setup analysis: interaction effects ####
    ## Loop through models ##
    for(m in 1:length(dataList_2)){
      
      cat("working on model: ", names(dataList_2)[m], "\n")

      ## 2.2.1 allocate empty vectors to store adjusted Cohen's D effect sizes, se, ci ##
      n.controls=rep(0,N_ROIs)
      n.PD=rep(0,N_ROIs)
      n.sites=rep(0,N_ROIs)
      p.value<-rep(0,N_ROIs)
      db.interaction=rep(0,N_ROIs)
      se.db.interaction=rep(0,N_ROIs)
      L.ci.db.interaction=rep(0,N_ROIs)
      U.ci.db.interaction=rep(0,N_ROIs)
      est.interaction=rep(NA,N_ROIs)
      se.beta.interaction=rep(0,N_ROIs)

      sink(paste0(names(dataList_2)[m],"_ModelSummaryStats.txt"))
      
      ## 2.2.2 create a named list to hold the fitted LME models for plotting
      fitlist <- as.list(1:N_ROIs)
      names(fitlist) <- names(agegroup_feature_data)[feature_col_indices]
      
      #### 2.3 Loop through all FreeSurfer features and run each LME model separately ####  
      for (i in seq_along(feature_col_indices)){
        
        # FreeSurfer column index of interest
        x = feature_col_indices[[i]]
        
        ### Remove sites that have less than N_min examples per class
        FreeSurfer_ROI_HC=data.frame(filter(dataList_2[[m]],(Dx==0)))[c(multi_site_col_ID, x)]
        FreeSurfer_ROI_PD=data.frame(filter(dataList_2[[m]],(Dx>0)))[c(multi_site_col_ID, x)]
        names(FreeSurfer_ROI_HC)[2]<-"feature_x"; names(FreeSurfer_ROI_PD)[2]<-"feature_x"
        non_na_HC<-setDT(FreeSurfer_ROI_HC)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
        non_na_PD<-setDT(FreeSurfer_ROI_PD)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
        
        filtered_dataList = copy(dataList_2[[m]])
        sites_to_exclude_HC<-data.frame(filter(non_na_HC,(non_na<N_min)))$MultiSiteID
        sites_to_exclude_PD<-data.frame(filter(non_na_PD,(non_na<N_min)))$MultiSiteID
        if (length((sites_to_exclude_HC)) > 0){
          filtered_dataList=subset(filtered_dataList, !(MultiSiteID%in%sites_to_exclude_HC & Dx==0))
        }
        if (length((sites_to_exclude_PD)) > 0){
          filtered_dataList=subset(filtered_dataList, !(MultiSiteID%in%sites_to_exclude_PD & Dx>0))
        }        
        # Check if remaining number of sites is >=N_min_sites, and at least one site has both classes included, else skip!
        filtered_tbl<-table(filtered_dataList$MultiSiteID,filtered_dataList$Dx)
        if (dim(filtered_tbl)[1] < N_min_sites | dim(filtered_tbl)[2]< 2 | !any(apply(filtered_tbl, 1, function(x) all(x>0)))){
          next
        }
        
        # FreeSurfer label of interest
        cat(names(filtered_dataList)[x], "\n")
        FreeSurfer_ROI=filtered_dataList[,x]
        
        # Define base model variables
        Dx<-as.factor(filtered_dataList$Dx)
        Sex<-as.factor(filtered_dataList$Sex)
        Age<-as.numeric(filtered_dataList$Age)
        AgeC<-Age-mean(Age)
        AgeC2<-AgeC^2
        ICV<-as.numeric(filtered_dataList$ICV)
        MultiSiteID<-as.factor(filtered_dataList$MultiSiteID)
        
        # Create new data frame that only includes variables used to fit LME model
        if (feature_label %in% list("subcortical_volume", "cortical_surface_area")){
          lme2_tmp_data <- data.frame(FreeSurfer_ROI, Dx, Sex, Age, AgeC, AgeC2, ICV, MultiSiteID)
        } else {
          lme2_tmp_data <- data.frame(FreeSurfer_ROI, Dx, Sex, Age, AgeC, AgeC2, MultiSiteID)
        }
        
        # Run LME model for interaction effect
        ctrl <- lmeControl(opt='optim');
        fml=as.formula(paste("FreeSurfer_ROI ~ ",paste(lm.List1[[m]], collapse = "+"),sep = ""))
        lme2=lme(fml, random= ~1 | MultiSiteID, na.action="na.exclude", control=ctrl, data=lme2_tmp_data)
        # This will ensure the full formula is stored inside the lme object so it can be accessed later
        lme2$call$fixed<-eval(lme2$call$fixed) 
        fit2=summary(lme2)
        lme2$call$data<-eval(lme2$call$data) # Add the actual data used to fit our model AFTER printing its summary
        lme2_data=retrieve.LME.sample(lme2)
        fitlist[[i]]<-lme2
        
        cat("Model summary statistics are: \n"); print(fit2); cat("\n")
        coefs2 <- data.frame(coefficients(fit2))
        n.controls[[i]]=length(which(lme2_data$Dx==0)); n.PD[[i]]=length(which(lme2_data$Dx==1))
        n.obs=fit2[["dims"]][["N"]]; n.groups=fit2[["dims"]][["ngrps"]][["MultiSiteID"]]
        n.sites[[i]]=n.groups
        param=length(fit2$coefficients$fixed)
        varcorr=VarCorr(fit2) # Extract intercept (random effect-stat) variance and residual variance 
        R.btwn=as.numeric(varcorr[which(rownames(varcorr) == "(Intercept)"), which(colnames(varcorr) == "Variance")])
        R.wthn=as.numeric(varcorr[which(rownames(varcorr) == "Residual"), which(colnames(varcorr) == "Variance")])
        R.calc=as.numeric(mixeff.R(R.btwn, R.wthn))
        interaction_label=tail(row.names(coefs2), 1) # Extract interaction effect label and T-value
        tval.interaction=coefs2[which(rownames(coefs2)==interaction_label), which(colnames(coefs2) == "t.value")]
        p.value[[i]]=coefs2[which(rownames(coefs2)==interaction_label), ]$p.value # pval for interaction term
        
        # # Calculate effect sizes using mixeff.d 
        db.interaction[[i]]=mixeff.d(tval.interaction, n.groups, n.obs, n.controls[[i]], n.PD[[i]], R.calc, 
                                     param) ## (k=num parameters incl intercept)
        se.db.interaction[[i]]=se.db(db.interaction[[i]], n.controls[[i]], n.PD[[i]])
        bound.db.interaction=CI1(db.interaction[[i]], se.db.interaction[[i]])
        L.ci.db.interaction[[i]]=bound.db.interaction[1]; U.ci.db.interaction[[i]]=bound.db.interaction[2]
        est.interaction[[i]]=coefs2[which(rownames(coefs2) == interaction_label), ]$Value
        se.beta.interaction[[i]]=coefs2[which(rownames(coefs2) == interaction_label), ]$Std.Error
        
      }
      sink()
      
      ## 2.4.1 calculate corrected p-values for overall model ##
      tmp.MCP.p.overall=p.adjust(p.value,method=MCP_method)
      MCP.p.overall=rep(0,N_ROIs)
      for (w in 1:N_ROIs){MCP.p.overall[w]=tmp.MCP.p.overall[w]}
      significant_MCP<-as.character(MCP.p.overall<MCP_alpha)
      
      # #### 2.5 Combine all output variables using Cohen's D into table ####
      outmatoverall=cbind(significant_MCP,db.interaction,se.db.interaction,L.ci.db.interaction,U.ci.db.interaction,
                          MCP.p.overall,p.value,n.PD,n.controls,n.sites,est.interaction,se.beta.interaction)

      rownames(outmatoverall)=names(agegroup_feature_data)[feature_col_indices]
    
      cat("Generating final output for model: ", names(dataList_2)[m],"\n")
      ## 2.5.1 Save output
      write.xlsx(outmatoverall, file = paste0("megaResults.xlsx"),sheetName=(names(dataList_2)[m]), append=TRUE)
      #### 2.5.2 Create plots for significant interaction effects
      sign_interaction_fits = fitlist[as.logical(significant_MCP)]
      if (length(sign_interaction_fits) > 0){
        for (i in seq_along(sign_interaction_fits)){
          FS_label<-names(sign_interaction_fits[i])
          interaction_term<-unlist(tail(lm.List1[[m]], 1)); 
          interaction_fit<-sign_interaction_fits[[i]]
          cat("Creating plot for significant interaction. Term:", interaction_term, ", Feature:", FS_label, "\n")
          plot.lme_interaction(interaction_fit, interaction_term, FS_label)
        }
      }
      rm(fitlist) # Remove LME models list to clear up memory
    }
    cat("\n")
    
    ###### SECTION 3: Sub-analyses: lme for continuous variable models using Pearson's r (models: 6d, 7abcde) ######
    cat("\nSection 2 completed. \n\nBeginning SECTION 3: Sub-analyses, Pearson's r for continous var models 6d, 7abcde\n\n")
    
    #### 3.1 Create datasets and set up variables for each model ####
    #### Select only patients with Age of Onset or questionnaire data (BAI, STAI, ASI)
    m6d.AgeO<-data.frame(filter(agegroup_feature_data,PD>0 & complete.cases(agegroup_feature_data$AgeO)))
    m7a.BAI<-data.frame(filter(agegroup_feature_data,PD>0 & complete.cases(agegroup_feature_data$BAI)))
    m7b.STAI<-data.frame(filter(agegroup_feature_data,PD>0 & complete.cases(agegroup_feature_data$STAI_T)))
    m7c.ASI<-data.frame(filter(agegroup_feature_data,PD>0 & complete.cases(agegroup_feature_data$ASI)))
    corr_dataSets<-list(m6d.AgeO,m7a.BAI,m7b.STAI,m7c.ASI)
    names(corr_dataSets)<-c("m6d.AgeO","m7a.BAI","m7b.STAI","m7c.ASI")
    
    # Define model variables as strings to paste into LME (with or without ICV) ##
    if (feature_label %in% list("subcortical_volume", "cortical_surface_area")){
      lm.AgeO<-list("AgeO","Sex","AgeC","AgeC2","ICV")
      lm.BAI<-list("BAI","Sex","AgeC","AgeC2","ICV")
      lm.STAI<-list("STAI_T","Sex","AgeC","AgeC2","ICV")
      lm.ASI<-list("ASI","Sex","AgeC","AgeC2","ICV")
    } else {
      lm.AgeO<-list("AgeO","Sex","AgeC","AgeC2")
      lm.BAI<-list("BAI","Sex","AgeC","AgeC2")
      lm.STAI<-list("STAI_T","Sex","AgeC","AgeC2")
      lm.ASI<-list("ASI","Sex","AgeC","AgeC2")
    }
    lm_corr.List<-list(lm.AgeO,lm.BAI,lm.STAI,lm.ASI)
    names(lm_corr.List)<-c("lm.AgeO","lm.BAI","lm.STAI","lm.ASI")
    
    ### Filter out sites with < N_min patients per (sub)-analysis across ALL regions within a given FreeSurfer modality ###
    dataList_3 = filter_datasets_section_3(corr_dataSets, N_sites, N_ROIs, N_min, N_min_sites)
    # Remove models with insufficient examples from lm_corr.List
    corr_dataList_to_remove = which(!names(corr_dataSets) %in% names(dataList_3))
    lm_corr.List[unlist(corr_dataList_to_remove)] <- NULL
    if (length(dataList_3) > 0){
      cat("Finished filtering corr models for:", feature_label, "\nFollowing models survived thresholding:", names(dataList_3), "\n\n")
    } else {
      cat("Finished filtering corr models for:", feature_label, "\nNO models survived thresholding! Skipping analyses... \n\n")
      next
    }
    
    #### 3.2 Setup analyses: lme for Continuous Variables ####
    ## Loop through correlation models ##
    for(c in 1:length(dataList_3)){
      
      cat("Working on model: ", names(dataList_3)[c], "\n")  
      
      # ####3.2.1 allocate empty vectors to store adjust effect sizes, se, ci ####
      n.PD=rep(0,N_ROIs)
      n.sites=rep(0,N_ROIs) 
      p.value<-rep(0,N_ROIs) 
      rME.PD=rep(0,N_ROIs) 
      se.rME.PD=rep(0,N_ROIs) 
      L.ci.rME.PD=rep(0,N_ROIs) 
      U.ci.rME.PD=rep(0,N_ROIs) 
      est.beta.corr=rep(NA,N_ROIs) 
      se.beta.corr=rep(0,N_ROIs) 
      
      sink(paste0(names(dataList_3)[c],"_ModelSummaryStats.txt"))
      
      #### 3.3 Loop through all FreeSurfer features and perform each regression ####  
      for (i in seq_along(feature_col_indices)){
        
        x = feature_col_indices[[i]]
        
        ### Remove sites that have less than (N_min) PD examples after removing NaNs
        FreeSurfer_ROI_PD=data.frame(filter(dataList_3[[c]]))[c(multi_site_col_ID, x)]
        names(FreeSurfer_ROI_PD)[2]<-"feature_x"
        non_na_PD<-setDT(FreeSurfer_ROI_PD)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
        filtered_dataList = copy(dataList_3[[c]])
        sites_to_exclude_PD<-data.frame(filter(non_na_PD,(non_na <  N_min)))$MultiSiteID
        if (length((sites_to_exclude_PD)) > 0){
          filtered_dataList=subset(filtered_dataList, !(MultiSiteID%in%sites_to_exclude_PD))
        }        
        # Check if remaining number of sites is at least N_min_sites, else skip!
        filtered_tbl<-table(filtered_dataList$MultiSiteID,filtered_dataList$Dx)
        if (dim(filtered_tbl)[1] < N_min_sites){
          next
        }
        
        cat(names(filtered_dataList)[x], "\n")
        FreeSurfer_ROI=filtered_dataList[,x]
        
        # Define base model variables
        Dx<-as.factor(filtered_dataList$Dx)
        Sex<-as.factor(filtered_dataList$Sex)
        Age<-as.numeric(filtered_dataList$Age)
        AgeC<-Age-mean(Age)
        AgeC2<-AgeC^2
        ICV<-as.numeric(filtered_dataList$ICV)
        MultiSiteID<-as.factor(filtered_dataList$MultiSiteID)
        AgeO<-as.numeric(filtered_dataList$AgeO)
        BAI<-as.numeric(filtered_dataList$BAI)
        STAI_T<-as.numeric(filtered_dataList$STAI_T)
        ASI<-as.numeric(filtered_dataList$ASI)
        
        # Combine variables in temporary dataframe used to fit LME model
        lme_vars<-unlist(c("FreeSurfer_ROI", lm_corr.List[[c]], "MultiSiteID"))
        lme3_tmp_data<-data.frame(mget(lme_vars))
        
        # Run LME model
        ctrl <- lmeControl(opt='optim');
        fml = as.formula(paste("FreeSurfer_ROI ~ ",paste(lm_corr.List[[c]],collapse="+"),sep=""))
        lme3=lme(fml, random= ~ 1 | MultiSiteID, na.action="na.exclude", control=ctrl, data=lme3_tmp_data)
        lme3$call$fixed<-eval(lme3$call$fixed) # This will ensure the full formula is stored in lme object
        lme_corr_fit=summary(lme3)
        
        cat("Model summary statistics are: \n"); print(lme_corr_fit); cat("\n")
        print(lme_corr_fit)
        coefs3 <- data.frame(coefficients(lme_corr_fit))
        corr_label=lm_corr.List[[c]][[1]]
        tval.corr=coefs3[which(rownames(coefs3) == corr_label), which(colnames(coefs3) == "t.value")]
        n.obs=lme_corr_fit[["dims"]][["N"]]; n.PD[[i]]=n.obs
        n.groups=lme_corr_fit[["dims"]][["ngrps"]][["MultiSiteID"]]; n.sites[[i]]=n.groups
        param=length(lme_corr_fit$coefficients$fixed)
        
        # Extract intercept (random effect-stat) variance and residual variance 
        varcorr=VarCorr(lme_corr_fit)
        R.btwn=as.numeric(varcorr[which(rownames(varcorr) == "(Intercept)"), which(colnames(varcorr) == "Variance")])
        R.wthn=as.numeric(varcorr[which(rownames(varcorr) == "Residual"), which(colnames(varcorr) == "Variance")])
        R.calc=as.numeric(mixeff.R(R.btwn,R.wthn))
        
        # Extract p-value for correlation of interest
        p.value[[i]]=coefs3[which(rownames(coefs3) == corr_label), ]$p.value
        
        #### 3.4 collect effect size data, standard error and 95 CI####
        rME.PD[[i]]=mePears.r(tval.corr,n.groups,n.obs,R.calc,param) ##AA: (k=num parameters incl intercept)
        se.rME.PD[[i]]=Zr.and.se2(rME.PD[[i]],n.obs)[2]
        bound.rME.PD=CI1(rME.PD[[i]],se.rME.PD[[i]])
        L.ci.rME.PD[[i]]=bound.rME.PD[1]
        U.ci.rME.PD[[i]]=bound.rME.PD[2]
        est.beta.corr[[i]]=lme_corr_fit$tTable[2,1]  
        se.beta.corr[[i]]=lme_corr_fit$tTable[2,2]
        cat("\n")
      }

      sink()
      
      ## 3.4.1 calculate corrected p-values for overall model ##
      tmp.MCP.p.overall=p.adjust(p.value,method=MCP_method)
      MCP.p.overall=rep(0,N_ROIs)
      for (w in 1:N_ROIs){MCP.p.overall[w]=tmp.MCP.p.overall[w]}
      ## Couple p-values to FreeSurfer features ##
      significant_MCP_corr<-as.character(MCP.p.overall<MCP_alpha)
      #### 3.5 Combine all output variables into table ####
      outmatoverall=cbind(significant_MCP_corr,rME.PD,se.rME.PD,L.ci.rME.PD,U.ci.rME.PD,MCP.p.overall,
                          p.value,n.PD,n.sites,est.beta.corr,se.beta.corr)
      rownames(outmatoverall)=names(agegroup_feature_data)[feature_col_indices]
      cat("Generating final output for model: ", names(dataList_3)[c],"\n")
      ## 3.5.1 Save output ##
      write.xlsx(outmatoverall, file = paste0("megaResults.xlsx"),sheetName=(names(dataList_3)[c]), append=TRUE)
    }
    
    cat("\n")

    ###### SECTION 4: Post-hoc between patient group sub-analyses ######
    cat("\nSection 3 completed. \n\nBeginning SECTION 4: Post-hoc between patient group sub-analyses")
    
    s1.Cur<-data.frame(filter(agegroup_feature_data, (PD==2) | (PD==1)))
    s1.Cur$PD[s1.Cur$PD==1]=0; s1.Cur$PD[s1.Cur$PD==2]=1
    s2.ANX<-data.frame(filter(agegroup_feature_data, (PD>0 & ComANXD==1) | (PD>0 & ComANXD==0)))
    s3.MDD<-data.frame(filter(agegroup_feature_data, (PD>0 & MDD>0) | (PD>0 & MDD==0)))
    s3.MDD$MDD.dich[s3.MDD$MDD==0]=0; s3.MDD$MDD.dich[s3.MDD$MDD>0]=1
    s4.Meds<-data.frame(filter(agegroup_feature_data, (PD>0 & Med==2) | (PD>0 & Med==1)))
    s4.Meds$Med[s4.Meds$Med==1]=0; s4.Meds$Med[s4.Meds$Med==2]=1
    s5.SSRI<-data.frame(filter(agegroup_feature_data, (PD>0 & SSRI_SNRI==1) | (PD>0 & SSRI_SNRI==0)))
    s6.EarlyOns<-data.frame(filter(agegroup_feature_data, (PD>0 & AgeO<=21) | (PD>0 & AgeO>21)))
    s6.EarlyOns$AO[s6.EarlyOns$AgeO<=21]=0; s6.EarlyOns$AO[s6.EarlyOns$AgeO>21]=1
    
    datasets<-list(s1.Cur, s2.ANX, s3.MDD, s4.Meds, s5.SSRI, s6.EarlyOns)
    names(datasets)<- c("s1.Cur", "s2.ANX", "s3.MDD", "s4.Meds", "s5.SSRI","s6.EarlyOns")
      
    cat("\nModels in section 1 are:\n", names(datasets), "\n")
    
    # Define model variables as strings to paste into LME (with or without ICV) ##
    if (feature_label %in% list("subcortical_volume", "cortical_surface_area")){
      lm.Cur<-list("PD","Sex","AgeC","AgeC2","ICV")
      lm.ANX<-list("ComANXD","Sex","AgeC","AgeC2","ICV")
      lm.MDD<-list("MDD.dich","Sex","AgeC","AgeC2","ICV")
      lm.Meds<-list("Med","Sex","AgeC","AgeC2","ICV")
      lm.SSRI<-list("SSRI_SNRI","Sex","AgeC","AgeC2","ICV")
      lm.EarlyOns<-list("AO","Sex","AgeC","AgeC2","ICV")
    } else {
      lm.Cur<-list("PD","Sex","AgeC","AgeC2")
      lm.ANX<-list("ComANXD","Sex","AgeC","AgeC2")
      lm.MDD<-list("MDD.dich","Sex","AgeC","AgeC2")
      lm.Meds<-list("Med","Sex","AgeC","AgeC2")
      lm.SSRI<-list("SSRI_SNRI","Sex","AgeC","AgeC2")
      lm.EarlyOns<-list("AO","Sex","AgeC","AgeC2")
    }
    
    lm_sub.List<-list(lm.Cur,lm.ANX,lm.MDD,lm.Meds,lm.SSRI,lm.EarlyOns)
    names(lm_sub.List)<-c("lm.Cur","lm.ANX","lm.MDD","lm.Meds","lm.SSRI","lm.EarlyOns")

    ## Filter out sites with < N_min patients per (sub)-analysis across ALL regions within a given FreeSurfer modality ###
    dataList_4 = filter_datasets_section_3(datasets, N_sites, N_ROIs, N_min, N_min_sites)
    ## Remove models with insufficient examples from lm_sub.List
    dataList_4_to_remove = which(!names(datasets) %in% names(dataList_4))
    lm_sub.List[unlist(dataList_4_to_remove)] <- NULL
    if (length(dataList_4) > 0){
      cat("Finished filtering datasets for:", feature_label, "\nFollowing models survived thresholding:", names(dataList_4), "\n\n")
    } else {
      cat("Finished filtering datasets for:", feature_label, "\nNO models survived thresholding! Skipping analyses... \n\n")
      next
    }
    
    for(l in 1:length(dataList_4)){
      
      ## Run LME models for each (filtered) dataset separately ##
      cat("working on model: ", names(dataList_4)[l], "\n") 
      
      ## Initialize arrays to store LME fit results
      n.p1=rep(0,N_ROIs)
      n.p2=rep(0,N_ROIs)
      n.sites=rep(0,N_ROIs)
      p.value<-rep(0,N_ROIs)
      db.sub=rep(0,N_ROIs)
      se.db.sub=rep(0,N_ROIs)
      perc.diff.sub=rep(0,N_ROIs)
      L.ci.db.sub=rep(0,N_ROIs)
      U.ci.db.sub=rep(0,N_ROIs)
      est.Effsub=rep(NA,N_ROIs)
      se.beta.sub=rep(0,N_ROIs)
      
      ## Parse effect of interest
      effect_of_interest = unlist(lm_sub.List[[l]][1])
      
      ## Loop through all FreeSurfer features in set and run each LME model separately
      sink(paste0(names(dataList_4)[l],"_ModelSummaryStats.txt"))
      
      for (i in seq_along(feature_col_indices)){
        
        # FreeSurfer column index of interest
        x = feature_col_indices[[i]]
        
        # Remove sites that have less than N_min examples per class
        FreeSurfer_ROI_P1=data.frame(filter(dataList_4[[l]], (get(effect_of_interest)==0)))[c(multi_site_col_ID, x)]
        FreeSurfer_ROI_P2=data.frame(filter(dataList_4[[l]], (get(effect_of_interest)==1)))[c(multi_site_col_ID, x)]
        names(FreeSurfer_ROI_P1)[2]<-"feature_x"; names(FreeSurfer_ROI_P2)[2]<-"feature_x"
        non_na_P1<-setDT(FreeSurfer_ROI_P1)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
        non_na_P2<-setDT(FreeSurfer_ROI_P2)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
        
        # Check if remaining number of sites is >=N_min_sites, and at least one site has both classes included, else skip model
        filtered_dataList = copy(dataList_4[[l]])
        sites_to_exclude_P1<-data.frame(filter(non_na_P1,(non_na<N_min)))$MultiSiteID
        sites_to_exclude_P2<-data.frame(filter(non_na_P2,(non_na<N_min)))$MultiSiteID
        if (length((sites_to_exclude_P1)) > 0){
          filtered_dataList=subset(filtered_dataList, !(MultiSiteID%in%sites_to_exclude_P1 & (get(effect_of_interest)==0)))
        }
        if (length((sites_to_exclude_P2)) > 0){
          filtered_dataList=subset(filtered_dataList, !(MultiSiteID%in%sites_to_exclude_P2 & (get(effect_of_interest)==0)))
        }        
        filtered_tbl<-table(filtered_dataList$MultiSiteID, filtered_dataList[c(effect_of_interest)][[1]])
        if (dim(filtered_tbl)[1] < N_min_sites | dim(filtered_tbl)[2] < 2 | !any(apply(filtered_tbl, 1, function(x) all(x>0)))){
          cat("Skipping model for :", names(filtered_dataList)[x], "\nN remaining number of sites is not >=N_min_sites, or no site that has both classes included \n\n")
          next
        }
        
        # FreeSurfer label of interest
        cat(names(filtered_dataList)[x], "\n")
        FreeSurfer_ROI=filtered_dataList[,x]
        
        # Define base model variables
        PD<-as.factor(filtered_dataList$PD) 
        ComANXD<-as.factor(filtered_dataList$ComANXD)
        MDD.dich<-as.factor(filtered_dataList$MDD.dich)
        SSRI_SNRI<-as.factor(filtered_dataList$SSRI_SNRI)
        Med<-as.factor(filtered_dataList$Med)
        AO<-as.factor(filtered_dataList$AO)
        Sex<-as.factor(filtered_dataList$Sex)
        Age<-as.numeric(filtered_dataList$Age)
        AgeC<-Age-mean(Age)
        AgeC2<-AgeC^2
        ICV<-as.numeric(filtered_dataList$ICV)
        MultiSiteID<-as.factor(filtered_dataList$MultiSiteID)
        
        ctrl <- lmeControl(opt='optim');
        
        # Combine variables in temporary dataframe used to fit LME model
        lme_vars<-unlist(c("FreeSurfer_ROI", lm_sub.List[[l]], "MultiSiteID"))
    
        # LME models for subcortical volumes and cortical surface area include TIV
        if (feature_label %in% list("subcortical_volume", "cortical_surface_area")){
          lme4_tmp_data<-data.frame(mget(lme_vars),ICV)
          fml = as.formula(paste("FreeSurfer_ROI ~ ",paste(lm_sub.List[[l]],collapse="+"),sep="")) 
          lme4=lme(fml, random= ~ 1 | MultiSiteID, na.action="na.exclude", control=ctrl, data=lme4_tmp_data)
        } else {
          lme4_tmp_data<-data.frame(mget(lme_vars))
          fml = as.formula(paste("FreeSurfer_ROI ~ ",paste(lm_sub.List[[l]],collapse="+"),sep="")) 
          lme4=lme(fml, random= ~ 1 | MultiSiteID, na.action="na.exclude", control=ctrl, data=lme4_tmp_data)
        }
        
        fit4=summary(lme4); lsm4=summary(lsmeans(lme4, effect_of_interest))
        
        lme4_data=retrieve.LME.sample(lme4)
        cat("Model summary statistics are: \n"); print(fit4)
        coefs4 <- data.frame(coefficients(fit4))
        p.value[[i]]=coefs4[2,5] 
        n.p1[[i]]=length(which(lme4_data[,2]==0)); n.p2[[i]]=length(which(lme4_data[,2]==1))
        n.obs=fit4[["dims"]][["N"]]; n.groups=fit4[["dims"]][["ngrps"]][["MultiSiteID"]]
        n.sites[[i]]=n.groups; 
        param=length(fit4$coefficients$fixed) # (k=num parameters incl intercept)
        perc.diff.sub[[i]]=lsm(lsm4[2,2],lsm4[1,2]) 
        varcorr=VarCorr(fit4) # Extract intercept (random effect-stat) variance and residual variance 
        R.btwn=as.numeric(varcorr[which(rownames(varcorr) == "(Intercept)"), which(colnames(varcorr) == "Variance")])
        R.wthn=as.numeric(varcorr[which(rownames(varcorr) == "Residual"), which(colnames(varcorr) == "Variance")])
        R.calc=as.numeric(mixeff.R(R.btwn,R.wthn))
       
         # Extract T value for effect of interest
        tval.Eff=fit4$tTable[2, which(colnames(fit4$tTable)=="t-value")]
        
        # Calculate effect size, standard error and 95 CI
        db.sub[[i]]=mixeff.d(tval.Eff, n.groups, n.obs, n.p1[[i]], n.p2[[i]], R.calc, param)
        se.db.sub[[i]]=se.db(db.sub[[i]],n.p1[[i]],n.p2[[i]])
        bound.db.sub=CI1(db.sub[[i]],se.db.sub[[i]])
        L.ci.db.sub[[i]]=bound.db.sub[1]; U.ci.db.sub[[i]]=bound.db.sub[2]
        est.Effsub[[i]]=coefs4[2, ]$Value
        se.beta.sub[[i]]=coefs4[2, ]$Std.Error
        cat("\n")
      }  
      
      sink()
      
      ## Calculate corrected p-values for overall model
      tmp.MCP.p.overall=p.adjust(p.value,method=MCP_method)
      MCP.p.overall=rep(0,N_ROIs)
      for (w in 1:N_ROIs){MCP.p.overall[w]=tmp.MCP.p.overall[w]}
      ## Couple p-values to FreeSurfer features
      significant_MCP_corr<-as.character(MCP.p.overall<MCP_alpha)
      ## Combine all output variables into table
      outmatoverall=cbind(significant_MCP_corr, db.sub, se.db.sub, perc.diff.sub, 
                          L.ci.db.sub, U.ci.db.sub, MCP.p.overall,
                          p.value, n.p1, n.p2, n.sites, est.Effsub, se.beta.sub)
      rownames(outmatoverall)=names(agegroup_feature_data)[feature_col_indices]
      cat("Generating final output for model: ", names(dataList_4)[l],"\n")
      ## Save output
      write.xlsx(outmatoverall, file = paste0("megaResults.xlsx"),sheetName=(names(dataList_4)[l]), append=TRUE)
    }
    cat("\n")
  }
}
sink()
setwd(mainDir)
    