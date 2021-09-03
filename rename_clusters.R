if (CL_FILTER != "") 
	FILTERED_CELLS <- FilterOutClusters(object = object, 
                                        id = CL_MODE, 
                                        cl_filter = CL_FILTER)

filtered_renamed_clusters <- NULL
if (RENAMING != "") {
	filtered_renamed_clusters <- unlist(strsplit(RENAMING, split=","))
	object <- RenameClusters(object = object, 
                             renaming = RENAMING, 
                             id = CL_MODE, 
                             cl_filter = CL_FILTER)
}
