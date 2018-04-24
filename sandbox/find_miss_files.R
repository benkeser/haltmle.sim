find_missing_files <- function(tag = "out",
                               parm = c("n","seed","r2"),
                               parm_values = list(n = c(200,1000),
                                                  seed = 1:1000,
                                                  r2 = c(0.01, seq(0.101, 0.901, by = 0.1))),
                               full_parm = NULL,
                               folder = "~/haltmle.sim/out"){
	all_files <- list.files(folder)
	tag_files <- all_files[grepl(tag, all_files)]
	# remove allOut files
	tag_files <- tag_files[!grepl("allOut", tag_files)]
	split_list <- sapply(tag_files,function(x){
		strsplit(tools::file_path_sans_ext(x),"_") 
		})
	l_x <- length(split_list[[1]])
	parm_list <- lapply(split_list, function(x){
		as.numeric(sapply(x[2:l_x], function(y){ strsplit(y, "=")[[1]][2] }))
	})
	done_parm <- data.frame(Reduce(rbind, parm_list))
	names(done_parm) <- parm
	if(is.null(full_parm)){
		full_parm <- expand.grid(parm_values)
	}
	rbinded_parm <- rbind(done_parm,full_parm)
	left_parm <- rbinded_parm[!(duplicated(rbinded_parm) | duplicated(rbinded_parm, fromLast = TRUE)), ]

	return(left_parm)
}


# # rename ks files
# ks_files <- all_files[grepl("ks_finalorig",all_files)]
# for(i in ks_files){
# 	new_name <- paste0("ksfinalorig", strsplit(i, "ks_finalorig")[[1]][2])
# 	file.rename(paste0("~/haltmle.sim/out/",i), paste0("~/haltmle.sim/out/",new_name))
# }
