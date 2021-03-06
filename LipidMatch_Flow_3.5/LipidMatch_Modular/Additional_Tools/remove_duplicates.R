#########################################################
# Lipidomics Workflow: Remove Duplicates                #
# Author: Andrew C. Patt                                #
# Email: patt.14@buckeyemail.osu.edu                    #
#########################################################

# Algorithm steps:
# 1. Get top ranked lipid per row
# 2. Find "," replace with "" in all names
# 3. Remove rows that have 4_ or 5_ (ones IDed by exact mass or not IDed)
# 4. Remove first 2 characters in names (MS/MS status)
# 5. Remove adducts from names:
#   find "+*" replace with ""
#   find " *" replace with ""
#   find "-H*" replace with ""
#   find "-2H*" replace with ""
# 6. Calculate median peak area or peak height
# 7. Identify highest median row for each group
# 8. Get adducts for all rows in that group and store those that fell within the retention time window
# 9. Return a data frame with just those rows
# 10.Remove rows that have sodium adducts

rm( list = ls() )

remove_duplicates<-function(data,row_start,adduct_col,annotation_col,window){
  #remove any commas in annotations
  data[,annotation_col]<-gsub(",","",data[,annotation_col])
  # create a vector with just the ranking score (1_, 2_, 3_, or 4_)
  ranking_class<-substr(as.vector(data[,annotation_col]),start=1,stop=2)
  # Identify and remove rows with 4_ or 5_
  not_MSMS<-sapply(ranking_class,function(x){
    if(x=="4_"||x=="5_"){
      return(FALSE)
    }else{
      return(TRUE)
    }
  })
  data<-data[which(not_MSMS),]
  # Calculate median peak area or peak height
  peak_cols<-which(grepl("Peak.area|Peak.height",colnames(data)))
  if (length(peak_cols)>1) {
    median_intensities<-apply(data[,peak_cols],1,median)
  }

  # A specific case where there is only one sample
  if (length(peak_cols)==1){
    median_intensities<-data[,peak_cols]
  }

  # Take top ranked lipid for each row
  top_ranked_lipid<-sapply(as.vector(data[,annotation_col]), function(x){
    strsplit(x,split = " | ")[[1]][1]
  })
  top_ranked_lipid<-as.vector(top_ranked_lipid)

  # Remove first 2 digits from lipid names
  top_ranked_lipid<-substring(top_ranked_lipid,3)

  # Remove adducts from names
  top_ranked_lipid<-gsub("\\+.*| .*|\\-H.*|\\-2H.*","",top_ranked_lipid)
  names(median_intensities)<-top_ranked_lipid
  # Identify sodiated rows
  not_sodiated<-sapply(as.vector(data[,adduct_col]),function(x){
    # Yang: if x is missing return TRUE, need to verify with Jeremy
    if (is.na(x)){
      return(TRUE)
    }
    else{
      if(x=="[M+Na]+"){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }
  })

  # Identify row # of most abundant duplicate for each group. Make a string of adducts within the window.
  keepers<-c()
  adducts<-c()
  for(i in unique(top_ranked_lipid)){
    rows = which(top_ranked_lipid==i)
    if(length(rows)==1){
      keepers<-c(keepers,rows)
      adducts<-c(adducts,as.vector(data[,adduct_col])[rows])
    }else{
      candidates<-median_intensities[rows]
      ranks<-order(candidates, decreasing = TRUE)

      # Prefer negative adducts if they're there
      all_adducts<-as.vector(data[,adduct_col])[rows]
      ### EDITED JPK (added new adducts for negative mode)
      if(grepl("+", paste(all_adducts,sep = "",collapse = "|")) && ("[M-H]-" %in% all_adducts || "[M-2H]-" %in% all_adducts || "[M+C2H3O2]-" %in% all_adducts || "[M+HCO2]-" %in% all_adducts || "[M+Cl]-" %in% all_adducts)){
        negative_adduct_indices<-c(match("[M-H]-",all_adducts),match("[M-2H]-",all_adducts),match("[M-2H]-",all_adducts),match("[M+C2H3O2]-",all_adducts),match("[M+HCO2]-",all_adducts),match("[M+Cl]-",all_adducts))
        negative_adduct_indices<-negative_adduct_indices[!is.na(negative_adduct_indices)]
        candidates_neg<-candidates[negative_adduct_indices]
        best_abundance<-max(candidates_neg)
        best_index<-intersect(which(median_intensities==best_abundance),rows)[1]
      } else {
        best_abundance<-max(candidates)
        best_index<-intersect(which(median_intensities==best_abundance),rows)[1]
      }

      # Use the runner-up if the top ranked adduct is sodium
      if(not_sodiated[best_index]==FALSE){
        best_index<-intersect(which(median_intensities==candidates[match(2,ranks)]),rows)[1]
      }
      # Find rows within RT window
      rows<-rows[ranks]
      in_window<-intersect(which(as.vector(data[,"row.retention.time"])[rows] <= as.vector(data[,"row.retention.time"])[best_index] + window/2),
                           which(as.vector(data[,"row.retention.time"])[rows] >= as.vector(data[,"row.retention.time"])[best_index] - window/2))
      in_window_candidates<-candidates[in_window]

      # Get adducts within that window
      adduct_indices<-rows[in_window]
      adduct_names<-unique(as.vector(data[,adduct_col])[adduct_indices])

      # Move negative adducts to the front of the list and sort by intensity
      if(grepl("+", paste(all_adducts,sep = "",collapse = "|")) && ("[M-H]-" %in% adduct_names || "[M-2H]-" %in% adduct_names || "[M+C2H3O2]-" %in% all_adducts || "[M+HCO2]-" %in% all_adducts || "[M+Cl]-" %in% all_adducts)){
        negative_adduct_indices<-c(match("[M-H]-",adduct_names),match("[M-2H]-",adduct_names),match("[M+C2H3O2]-",adduct_names),match("[M+HCO2]-",adduct_names),match("[M+Cl]-",adduct_names))
        negative_adduct_indices<-negative_adduct_indices[!is.na(negative_adduct_indices)]
        negative_adduct_indices<-negative_adduct_indices[order(in_window_candidates[negative_adduct_indices],decreasing = TRUE)]
        adduct_names<-c(adduct_names[negative_adduct_indices],adduct_names[-negative_adduct_indices])
      }
      adducts<-c(adducts,paste(adduct_names,sep = "",collapse = "|"))
      keepers<-c(keepers,best_index)
    }
  }

  # Return a frame with no duplicates
  data<-data.frame(data[keepers,],top_ranked_lipid[keepers],adducts)
  data<-data[which(as.vector(data[,ncol(data)])!="[M+Na]+"),]
  colnames(data)[ncol(data)-1] = "Molecular"
  colnames(data)[ncol(data)] = "Adducts Confirmed by MS/MS"
  #write.csv(data,paste(substr(Directory,1,nchar(Directory)-4),"_Molecular.csv",sep=""), row.names=FALSE, quote=FALSE)
  return(data)
}

#################################################################################################
# Merge modes is a separate function that can be used to merge two data frames of different modes
# The output of this function can be used for remove_duplicates
#################################################################################################
merge_modes<-function(data_pos,data_neg){
  #remove _pos and _neg from colnames
  colnames(data_pos)<-gsub("_pos","",colnames(data_pos),ignore.case = TRUE)
  colnames(data_neg)<-gsub("_neg","",colnames(data_neg),ignore.case = TRUE)
  colnames(data_pos)<-gsub("_p","",colnames(data_pos),ignore.case = TRUE)
  colnames(data_neg)<-gsub("_n","",colnames(data_neg),ignore.case = TRUE)

  #Store original number of columns in positive and negative mode
  original_pos_col<-length(colnames(data_pos))
  original_neg_col<-length(colnames(data_neg))
  #take out columns found in one frame and not the other
  if(length(setdiff(colnames(data_pos),colnames(data_neg)))>0){
    data_pos<-data_pos[,-match(intersect(colnames(data_pos),setdiff(colnames(data_pos),colnames(data_neg))),colnames(data_pos))]
  }
  if(length(setdiff(colnames(data_pos),colnames(data_neg)))>0){
    data_neg<-data_neg[,-match(intersect(colnames(data_neg),setdiff(colnames(data_pos),colnames(data_neg))),colnames(data_neg))]
  }
  #reorder neg to match pos column order
  data_neg<-data_neg[,match(colnames(data_pos),colnames(data_neg))]
  #rbind pos and neg.
  data_merged<-rbind(data_pos,data_neg)
  peak_cols<-which(grepl("Peak.area|Peak.height",colnames(data_merged)))
  # Print warning and stop code if there are no matching file names in positive and negative mode
  if (length(peak_cols)<1) {
    stop("No samples found in common between modes")
  }
  # Print warning if there are files without matching names in both positive and negative mode
  if((original_pos_col-length(colnames(data_pos)))>0){
    print(paste("Warning: ",(original_pos_col-length(colnames(data_pos)))," files in positive mode did not have matching files in negative mode. Make sure names are exactly the same, except _neg and _pos, or _n and _p",sep=""))
  }
  if((original_neg_col-length(colnames(data_neg)))>0){
    print(paste("Warning: ",(original_neg_col-length(colnames(data_neg)))," files in negative mode did not have matching files in positive mode. Make sure names are exactly the same, except _neg and _pos, or _n and _p",sep=""))
  }
  return(data_merged)
}

# Running
row_start<-1
window<-0.2
## Input 1,2 and 3
PosMode<-TRUE
NegMode<-TRUE
PosNegCombined<-TRUE

if(PosMode==TRUE) {
Directory_pos <- "C:/DESKTOP/DATA/2019_07_02_MeeraOxidized/LM_Flow_AllOils/MeeraOxidizedFinal/PosIDed.csv"
  data_pos <- read.csv(Directory_pos, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE,stringsAsFactors=FALSE) # 2018-02-02: Yang add stringsAsFactors=FALSE flag to force variables (such as Intensity_Ranked) to be characters instead of factors
  adduct_col<-ncol(data_pos)-1
  annotation_col<-ncol(data_pos)-4
  pos_output<-remove_duplicates(data_pos,row_start,adduct_col,annotation_col,window)
  write.table(pos_output, file = paste(paste(substr(Directory_pos,1,nchar(Directory_pos)-11),"Pos_MolecularSpecies.csv", sep = ""), sep = "/"), sep = ",", col.names = TRUE, row.names = FALSE)
}

if(NegMode==TRUE) {
Directory_neg <- "C:/DESKTOP/DATA/2019_07_02_MeeraOxidized/LM_Flow_AllOils/MeeraOxidizedFinal/NegIDed.csv"
  data_neg <- read.csv(Directory_neg, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE,stringsAsFactors=FALSE) # 2018-02-02: Yang add stringsAsFactors=FALSE flag to force variables (such as Intensity_Ranked) to be characters instead of factors
  adduct_col<-ncol(data_neg)-1
  annotation_col<-ncol(data_neg)-4
  neg_output<-remove_duplicates(data_neg,row_start,adduct_col,annotation_col,window)
  write.table(neg_output, file = paste(paste(substr(Directory_neg,1,nchar(Directory_neg)-11),"Neg_MolecularSpecies.csv", sep = ""), sep = "/"), sep = ",", col.names = TRUE, row.names = FALSE)
}

if(PosNegCombined==TRUE) {
  data_merged<-merge_modes(data_pos,data_neg)
  adduct_col<-ncol(data_merged)-1
  annotation_col<-ncol(data_merged)-4
  merged_output<-remove_duplicates(data_merged,row_start,adduct_col,annotation_col,window)
  write.table(merged_output, file = paste(paste(substr(Directory_neg,1,nchar(Directory_neg)-11),"PosNeg_MolecularSpecies.csv", sep = ""), sep = "/"), sep = ",", col.names = TRUE, row.names = FALSE)
}
