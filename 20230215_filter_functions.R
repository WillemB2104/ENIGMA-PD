filter_datasets_section_1<-function(datasets, N_sites, N_ROIs, N_min, N_min_sites){
  dataList<-list()
  dataList_to_remove<-list()
  for (d in 1:length(datasets)){
    cat("Filtering dataset:", names(datasets[d]), "\n")
    dataList[[d]]<-datasets[[d]]
    tbl<-data.frame(matrix(ncol=0,nrow=N_sites))
    tbl<-table(datasets[[d]]$MultiSiteID, datasets[[d]]$Dx)
    # Check if at least one site has BOTH groups (patients and controls)
    if(ncol(tbl) != 2){
      cat("There are NO sites that have data available for both patients and controls (N>=1). Removing model from analyses \n\n")
      dataList_to_remove<-append(dataList_to_remove, d)
      next
    }
    # Check if sites have at least ONE patient example (no requirement for controls)
    for (t in 1:nrow(tbl)){
      if(tbl[,2][t] < 1){  
        dataList[[d]]<-subset(dataList[[d]],MultiSiteID!=dimnames(tbl)[[1]][t])
      } 
    }
    tbl<-table(dataList[[d]]$MultiSiteID, dataList[[d]]$Dx)
    # Check if any site remains AFTER removing those with N<1 (step above)
    if(nrow(tbl) == 0){
      cat("There are NO sites that have data available for both patients and controls (N>=1). Removing model from analyses \n\n")
      dataList_to_remove<-append(dataList_to_remove, d)
      next
    }
    # Here we calculate the number of missing values per FreeSurfer feature, per class separately
    cgoodpd<-data.frame(matrix(ncol=0,nrow=nrow(tbl)))
    cgoodhcs<-data.frame(matrix(ncol=0,nrow=nrow(tbl)))
    tmppd<-data.frame(filter(dataList[[d]],(Dx>0)))
    tmphcs<-data.frame(filter(dataList[[d]],(Dx==0)))
    for (i in seq_along(feature_col_indices)){
      v = feature_col_indices[[i]]
      gpdvols<-data.frame(matrix(ncol=0,nrow=0))
      pdvols<-data.frame(tmppd)[,c(multi_site_col_ID,v)]
      names(pdvols)[2] <- "feature_x"
      gpdvols<-setDT(pdvols)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
      cgoodpd[,1]<-gpdvols[,c(1)]
      cgoodpd[,i+1]<-gpdvols[,c(2)]
    }
    for (i in seq_along(feature_col_indices)){
      v = feature_col_indices[[i]]
      ghcsvols<-data.frame(matrix(ncol=0,nrow=0))
      hcsvols<-data.frame(tmphcs)[,c(multi_site_col_ID,v)]
      names(hcsvols)[2] <- "feature_x"
      ghcsvols<-setDT(hcsvols)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
      # Add sites missing in ghcsvols overview (those that do not have control examples) with 0
      empty_sites=rownames(tbl)[!(rownames(tbl) %in% ghcsvols$MultiSiteID)]
      if(length(empty_sites) > 0){
        tmp=cbind(empty_sites, rep(0, length(empty_sites)))
        colnames(tmp)<-colnames(ghcsvols)
        ghcsvols<-rbind(ghcsvols, tmp)
        ghcsvols <- ghcsvols[order(ghcsvols$MultiSiteID),]
      }
      cgoodhcs[,1]<-ghcsvols[,c(1)]
      cgoodhcs[,i+1]<-ghcsvols[,c(2)]
    }
    # Sort cgoodhcs & cgoodpd on multisite ID to match with tbl
    cgoodhcs <- cgoodhcs[order(cgoodhcs$MultiSiteID),]
    cgoodpd <- cgoodpd[order(cgoodpd$MultiSiteID),]
    dx_nsite<-data.frame(matrix(ncol=0,nrow=nrow(tbl)))
    dx_nsite[,1]<-cgoodpd[,1]
    dx_nsite[,2]<-tbl[,2]
    dx_nsite[,3]<-apply(cgoodpd[,2:(1+N_ROIs)], 1, min)
    dx_nsite[,4]<-sapply(transpose(cgoodpd[,2:(1+N_ROIs)]), function(x) !all(x==0))
    dx_nsite[,5]<-cgoodhcs[,1]
    dx_nsite[,6]<-tbl[,1]
    dx_nsite[,7]<-apply(cgoodhcs[,2:(1+N_ROIs)], 1, min)
    dx_nsite[,8]<-sapply(transpose(cgoodhcs[,2:(1+N_ROIs)]), function(x) !all(x==0))
    colnames(dx_nsite)=c("SiteID", "N_patients", "min_PD_across_ROIs", "any_PD_with_non_NaN_features",
                         "SiteID", "N_controls", "min_HC_across_ROIs", "any_HC_with_non_NaN_features")
    write.xlsx(dx_nsite, file = paste0("site_sample_sizes_BEFORE_filtering.xlsx"),
               sheetName=names(datasets)[[d]], append=TRUE)
    # Filter out sites that have 0 non-NaN examples across ALL features within modality for PD patients
    for (t in 1:nrow(dx_nsite)){
      if(!dx_nsite[, which(colnames(dx_nsite) == 'any_PD_with_non_NaN_features')][t]){
        # Only filter for PD, not HC. These are sample sizes after QC
        dataList[[d]]<-subset(dataList[[d]],MultiSiteID !=dimnames(tbl)[[1]][t])
      } 
    }
    # Filter out sites that have less than N_min patients
    for (t in 1:nrow(dx_nsite)){
      if(dx_nsite[, which(colnames(dx_nsite) == 'N_patients')][t] < N_min){
        # Sample sizes after QC
        dataList[[d]]<-subset(dataList[[d]],MultiSiteID !=dimnames(tbl)[[1]][t])
      } 
    }
    tbl<-table(dataList[[d]]$MultiSiteID, dataList[[d]]$Dx)
    if(nrow(tbl) == 0){
      cat("Not enough data available after filtering (need at least", N_min, "PD per site). Removing model from analyses\n\n")
      dataList_to_remove<-append(dataList_to_remove, d)
    } else if (length(unique(dataList[[d]]$MultiSiteID)) < N_min_sites) {
      cat("Number of sites surviving thresholding is <", N_min_sites, ". Removing model from analyses\n\n")
      dataList_to_remove<-append(dataList_to_remove, d)
    } else {
      cat(names(datasets)[d], "\nSite count is: ",nrow(tbl), "\n Total Controls is: ", sum(tbl[,1]), "\nControls per site is: ", tbl[,1],"\n Total Patients is: ", sum(tbl[,2]), "\nPatients per site is: ", tbl[,2], "\n")
      #Save a separate overview for site's that are included after filtering
      dx_nsite_filtered<-subset(dx_nsite, SiteID%in%dimnames(tbl)[[1]])
      write.xlsx(dx_nsite_filtered, file = paste0("site_sample_sizes_AFTER_filtering.xlsx"), sheetName=names(datasets)[[d]], append=TRUE)
    }
    cat("\n")
  }
  names(dataList)<-names(datasets)
  
  # Remove models with insufficient examples from analyses specified in dataList
  dataList[unlist(dataList_to_remove)] <- NULL
  return(dataList)
}


filter_datasets_section_3<-function(corr_dataSets, N_sites, N_ROIs, N_min, N_min_sites){
  corr_dataList<-list()
  corr_dataList_to_remove<-list()
  for (w in 1:length(corr_dataSets)){
    cat("Filtering dataset:", names(corr_dataSets[w]), "\n")
    corr_dataList[[w]]<-corr_dataSets[[w]]
    tbl<-data.frame(matrix(ncol=0,nrow=N_sites))
    tbl<-table(corr_dataSets[[w]]$MultiSiteID, corr_dataSets[[w]]$Dx)
    # Check if sites have at least ONE patient example, otherwise exclude
    for (t in 1:nrow(tbl)){
      if(tbl[,1][t] < 1){  
        corr_dataList[[w]]<-subset(corr_dataList[[w]],MultiSiteID !=dimnames(tbl)[[1]][t])
      }
    }
    # Here we check if any data is available for analysis AFTER removing those with N<1 (step above)
    tbl<-table(corr_dataList[[w]]$MultiSiteID, corr_dataList[[w]]$Dx)
    if(nrow(tbl) == 0){
      cat("There are NO sites that have data available for patients (N>=1). Removing model from analyses \n\n")
      corr_dataList_to_remove<-append(corr_dataList_to_remove, w)
      next
    }
    # Calculate the number of missing values per FreeSurfer feature (for PD only)
    cgoodpd<-data.frame(matrix(ncol=0,nrow=nrow(tbl)))
    tmppd<-data.frame(filter(corr_dataList[[w]],(Dx>0)))
    for (i in seq_along(feature_col_indices)){
      v = feature_col_indices[[i]]
      gpdvols<-data.frame(matrix(ncol=0,nrow=0))
      pdvols<-data.frame(tmppd)[,c(multi_site_col_ID,v)]
      names(pdvols)[2] <- "feature_x"
      gpdvols<-setDT(pdvols)[, .(non_na = sum(!is.na(feature_x))), MultiSiteID]
      cgoodpd[,1]<-gpdvols[,c(1)]
      cgoodpd[,i+1]<-gpdvols[,c(2)]
    }
    cgoodpd <- cgoodpd[order(cgoodpd$MultiSiteID),]
    corr_nsite<-data.frame(matrix(ncol=0,nrow=nrow(tbl)))
    corr_nsite[,1]<-cgoodpd[,1]
    corr_nsite[,2]<-tbl[,1]
    corr_nsite[,3]<-apply(cgoodpd[,2:(1+N_ROIs)], 1, min)
    corr_nsite[,4]<-sapply(transpose(cgoodpd[,2:(1+N_ROIs)]), function(x) !all(x==0))
    colnames(corr_nsite)=c("SiteID", "N_patients", "min_N_across_ROIs", "any_PD_with_non_NaN_features")
    write.xlsx(corr_nsite, file = paste0("site_sample_sizes_BEFORE_filtering.xlsx"), 
               sheetName=names(corr_dataSets)[[w]], append=TRUE)
    
    # Filter out sites that have 0 non-NaN examples across ALL features for PD patients
    for (t in 1:nrow(corr_nsite)){
      if(!corr_nsite[, which(colnames(corr_nsite) == 'any_PD_with_non_NaN_features')][t]){
        # Only filter for PD, not HC. These are sample sizes after QC
        corr_dataList[[w]]<-subset(corr_dataList[[w]],MultiSiteID !=dimnames(tbl)[[1]][t])
      }}
    # Filter out sites that have less than N_min PD examples
    for (t in 1:nrow(corr_nsite)){
      if(corr_nsite[, which(colnames(corr_nsite) == 'N_patients')][t] < N_min){
        # Sample sizes after QC
        corr_dataList[[w]]<-subset(corr_dataList[[w]],MultiSiteID !=dimnames(tbl)[[1]][t])
      }}
    tbl<-table(corr_dataList[[w]]$MultiSiteID, corr_dataList[[w]]$Dx)
    if(nrow(tbl) == 0){
      cat("Not enough data available after filtering (need at least", N_min, "PD per site). Removing model from analyses\n\n")
      corr_dataList_to_remove<-append(corr_dataList_to_remove, w)
    } else if (length(unique(corr_dataList[[w]]$MultiSiteID)) < N_min_sites) {
      cat("Number of sites surviving thresholding is <", N_min_sites, ". Removing model from analyses\n\n")
      corr_dataList_to_remove<-append(corr_dataList_to_remove, w)
    } else {
      cat(names(corr_dataSets)[w], "\nSite count is: ",nrow(tbl), "\n Total Patients is: ", sum(tbl[,1]), "\nPatients per site is: ", tbl[,1], "\n")
      #Save a separate overview for site's that are included after filtering
      corr_nsite_filtered<-subset(corr_nsite, SiteID%in%dimnames(tbl)[[1]])
      write.xlsx(corr_nsite_filtered, file = paste0("site_sample_sizes_AFTER_filtering.xlsx"), 
                 sheetName=names(corr_dataSets)[[w]], append=TRUE)
    }
    cat("\n")
  }
  names(corr_dataList)<-names(corr_dataSets)
  
  # Remove models with insufficient examples from analyses specified in corr_dataList
  corr_dataList[unlist(corr_dataList_to_remove)] <- NULL
  return(corr_dataList)
}