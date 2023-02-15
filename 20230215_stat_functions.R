#### Functions needed for effect statistics (d and r), SE for each effect statistic, and CI calculations ####
## effect size calculations for linear mixed effects models based on Nakagawa & Cuthill (2007) ##
## equation numbers refer to corresponding equation in above paper ##

## repeatability/ICC value required for effect stat ##
mixeff.R<-function(s2.bt,s2.wt){
  R<-s2.bt/(s2.bt+s2.wt)
  names(R)<-"R for d and r (eq23)"
  return(R)
}

## lme group comparisons ##
mixeff.d<-function(t.val,n.i,n.o,n.o1,n.o2,R,k){
  d<-(t.val*(1+(n.i/n.o)*R))*((sqrt(1-R))*n.o)/
    (sqrt(n.o1*n.o2)*sqrt(n.o-k))
  names(d)<-"effect size d (specific for lm with random effects (eq22)"
  return(d)
}

se.db<-function(d,n1,n2){
  seCo<-sqrt(((n1+n2-1)/(n1+n2-3))*(4/(n1+n2))*(1+((d^2)/8)))
  names(seCo)<-"se for Cohens d (eq16)"
  return(seCo)
}

## lme interactions and linear regressions ## 
mePears.r<-function(t.val,n.i.,n.o.,R,k){
  r<-(t.val*(1+(n.i./n.o.)*R)*sqrt(1-R))/
    (sqrt(t.val^2*(1+(n.i./n.o.)*R)^2*sqrt(1-R)+n.o.-k))
  names(r)<-"effect size r (specific for lm with random effects (eq24)"
  return(r)
}

Zr.and.se2<-function(r,n){
  Zr<-0.5*log((1+r)/(1-r))
  se<-(1/sqrt(n-3))
  names(Zr)<-"Fisher's z"
  names(se)<-"se for Zr"
  return(c(Zr,se))
}

## confidence intervals ##
CI1<-function(ES,se){
  ci<-c((ES-(1.96)*se),(ES+(1.96)*se))
  names(ci)<-c("95% CI lower","95% CI upper")
  return(ci)
}

## calculate percentage difference
lsm<-function(mpd,mctrl){
  perc<-(mpd-mctrl)/mctrl*100
  names(perc)<-"perc diff"
  return(perc)
}

## retrieve sample data from lme model where NAs are omitted
retrieve.LME.sample<-function(lme_model){
  lme_data_df = copy(lme_model$data)
  row_indices = 1:nrow(lme_data_df)
  indices_to_drop = as.vector(unlist(lme_model$na.action)) # These are row indices that contain NAs
  indices_to_keep = setdiff(row_indices, indices_to_drop)
  used_data_df = lme_data_df %>% slice(indices_to_keep)
  return(used_data_df)
}

## create plots for LME interaction effects
plot.lme_interaction<-function(interaction_fit, interaction_term, FS_label){
  interaction_eff<-effect(term=interaction_term, mod=interaction_fit, xlevels=list(Dx=c(0,1)))
  # The computed effect absorbs the lower-order terms marginal to the term in question, and averages over other terms in the model
  interaction_df<-as.data.frame(interaction_eff) #convert the effects list to a data frame
  lme_data_df<-retrieve.LME.sample(interaction_fit)
  lme_data_1<-lme_data_df[,!(names(lme_data_df) %in% c("MultiSiteID"))]
  lme_data_2<-interaction_eff$data
  lme_data_2$Age<-lme_data_1$Age
  # Check if data used to fit LME retrieved from LME model object and effects package are identical
  if(all_equal(lme_data_1, lme_data_2)){
    interaction_df$Dx<-as.factor(interaction_df$Dx)
    x_axis_term=tail(unlist(strsplit(interaction_term, ":")), 1)
    
  if(x_axis_term == 'AgeC2'){
      mean_Age = mean(lme_data_df$Age)
      interaction_df_pos_Age <- interaction_df
      interaction_df_pos_Age['Age'] = mean_Age + sqrt(interaction_df[x_axis_term])
      interaction_df_neg_Age <- interaction_df
      interaction_df_neg_Age['Age'] = mean_Age - sqrt(interaction_df[x_axis_term])
      interaction_df = rbind(interaction_df_neg_Age, interaction_df_pos_Age)
      x_axis_term = "Age"
    }
    
    # Plot interaction effect using LME fits
    ggplot(interaction_df, aes_string(x=x_axis_term, y="fit", color="Dx", group="Dx")) + geom_point() + 
      geom_line(size=1, method=loess) + geom_ribbon(aes(ymin=lower, ymax=upper, fill=Dx), alpha=0.1, linetype=2) + 
      labs(y = FS_label) + labs(title = paste("Interaction for", interaction_term)) + 
      theme(plot.title = element_text(hjust = 0.5))
    filename_1 = sub(":", "-", paste("./interactions/interaction", interaction_term, FS_label, "fit", sep="_"))
    dir.create(dirname(filename_1), showWarnings = FALSE)
    ggsave(paste(filename_1, 'tiff', sep='.'), dpi=300)
    ggsave(paste(filename_1, 'pdf', sep='.'), dpi=300)
    # Store data used to create scatter plot
    write.csv(interaction_df, paste(filename_1, 'csv', sep='.'), row.names = FALSE) 
    
    # Create scatter plot for interaction effect
    lme_data_df$Dx<-as.factor(lme_data_df$Dx)
    filename_2 = sub(":", "-", paste("./interactions/interaction", interaction_term, FS_label, "scatter", sep="_"))
    ggplot(lme_data_df, aes_string(x=x_axis_term, y="FreeSurfer_ROI", color="Dx")) + 
      geom_jitter(width=2, alpha=0.5) + geom_smooth(method=loess, formula=y~x, fullrange=TRUE) + labs(y = FS_label) + 
      labs(title = paste("Interaction for", interaction_term)) + theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste(filename_2, 'tiff', sep='.'), dpi=300)
    ggsave(paste(filename_2, 'pdf', sep='.'), dpi=300)
    # Store data used to create scatter plot
    write.csv(lme_data_df, paste(filename_2, 'csv', sep='.'), row.names = FALSE) 
    
    # Plot interaction effect using LME fits and scatter of raw data # scale_x_sqrt(limits=c(-1000, 1000)) to scale
    filename_3 = sub(":", "-", paste("./interactions/interaction", interaction_term, FS_label, "fit_scatter", sep="_"))
    ggplot(interaction_df, aes_string(x=x_axis_term, y="fit", color="Dx", group="Dx")) + geom_point() +
      geom_line(size=1, method=loess) + geom_ribbon(aes(ymin=lower, ymax=upper, fill=Dx), alpha=0.1, linetype=2) +
      geom_jitter(data=lme_data_df, aes_string(x=x_axis_term, y="FreeSurfer_ROI", color="Dx"), alpha=0.2, width=2) + 
      labs(y = FS_label) + labs(title = paste("Interaction for", interaction_term)) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste(filename_3, 'tiff', sep='.'), dpi=300)
    ggsave(paste(filename_3, 'pdf', sep='.'), dpi=300)
    ##Store data used to create scatter fit plot
    write.csv(interaction_df, paste(filename_3, 'csv', sep='.'), row.names = FALSE) 
    
    # anova(interaction_fit, type="marginal")
  } else {
    cat("Warning: Could not create interaction effect plots. Check model data stored in LME object!")
  }
}

## Plot age distributions per site
plot.age_dist<-function(agegroup_data, N_min){
  
  agegroup_data$Dx<-as.factor(agegroup_data$Dx)
  agegroup_data$MultiSiteID<-as.factor(agegroup_data$MultiSiteID)
  
  # Exclude sites with <N_min patients as these are not included for analysis
  tbl<-table(agegroup_data$MultiSiteID, agegroup_data$Dx)
  excluded_sites = row.names(tbl)[tbl[,2] < N_min]
  filtered_df<-copy(agegroup_data)
  for (site in excluded_sites){filtered_df <- subset(filtered_df, MultiSiteID!=site)}
  
  # Plot age distributions
  ggplot(filtered_df, aes(x=Age, y=MultiSiteID)) + 
    geom_density_ridges(aes(fill = Dx, point_color = Dx), alpha=.3, scale=1, 
                        jittered_points = TRUE, point_alpha = 0.35, point_size = 0.25) + 
    scale_y_discrete(limits = rev)
  ggsave(paste('age_distribution_per_site', 'pdf', sep='.'), dpi=300, width=5.66, height=8.98, units='in')
  ggsave(paste('age_distribution_per_site', 'tiff', sep='.'), dpi=300, width=5.66, height=8.98, units='in')
}


