# Author: Tülay Karakulak
# Created: 04.01.2024
# Description: The script takes two gff files which are the output of MAS-Seq Iso-Seq workflow and matches the same transcripts based on their exon coordinates with some flexibility.


library(rtracklayer)
library(dplyr)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)

# Check if the correct number of arguments (8) is provided
if(length(args) != 8) {
	  stop("Please provide the paths for two input GFF files, 2 Seurat genes.tsv, number of nucleotides flexibility for 5 prime and 3prime, output directory and output file")
}

# Assigning input and output file paths from the arguments
input_path1 <- args[1]
input_path2 <- args[2]
genes_path1 <- args[3]
genes_path2 <- args[4]
max5diff <-  as.numeric(args[5])
max3diff <- as.numeric(args[6])
output_dir <- args[7]
output_filename <- args[8]


# Importing GFF files
ccRCC_gff <- import(input_path1)
ccRCC_2_gff <- import(input_path2)

# Reading and processing gene name files
isoform_gene_names1 <- read.csv(genes_path1, header = FALSE, sep='\t')
isoform_gene_names2 <- read.csv(genes_path2, header = FALSE, sep='\t')

colnames(isoform_gene_names1) <- c('PBid', 'GeneName')
colnames(isoform_gene_names2) <- c('PBid', 'GeneName')

# Cleaning up the PBid columns by removing everything after ":"
isoform_gene_names1$PBid <- sub(":.*$", "", isoform_gene_names1$PBid)
isoform_gene_names2$PBid <- sub(":.*$", "", isoform_gene_names2$PBid)

# Extracting exon information from GFF data
exons1 <- ccRCC_gff[ccRCC_gff$type == "exon",]
exons2 <- ccRCC_2_gff[ccRCC_2_gff$type == "exon",]

df_exons1 <- as.data.frame(exons1)
df_exons2 <- as.data.frame(exons2)

# Find number of exons for each transcript id
df_exons1_grouped <- df_exons1 %>% 
	  group_by(transcript_id) %>%
	  mutate(concatenated_coord_single = paste(start, end, sep = "-")) %>% 
	  mutate(concatanated_coord = paste(concatenated_coord_single, collapse = ',')) %>%  
	  ungroup() %>%
	  distinct(transcript_id, concatanated_coord, .keep_all = TRUE)

df_exons2_grouped <- df_exons2 %>% 
	group_by(transcript_id) %>%
	mutate(concatenated_coord_single = paste(start, end, sep = "-")) %>% 
	mutate(concatanated_coord = paste(concatenated_coord_single, collapse = ',')) %>%  
	ungroup() %>%
	distinct(transcript_id, concatanated_coord, .keep_all = TRUE)

df_exons1_grouped$NumberOfExons <- sapply(strsplit(df_exons1_grouped$concatanated_coord, ","), length)
df_exons2_grouped$NumberOfExons <- sapply(strsplit(df_exons2_grouped$concatanated_coord, ","), length)

# Function to create a range from concatenated coordinates
make_a_range_all <- function(coordinates) {
    	coords <- strsplit(coordinates, ",")[[1]]
    	ranges <- lapply(coords, function(coord) {
        parts <- as.integer(unlist(strsplit(coord, "-")))
        return(list(start = parts[1], end = parts[2]))
    })
    return(ranges)
}



check_criteria <- function(ranges1, ranges2) {
	    # Check if the number of ranges in each list is the same and non-zero
	 if(length(ranges1) != length(ranges2) || length(ranges1) == 0) {
		            return(FALSE)
    }

    # Check first coordinate with flexibility
	    if(!(ranges2[[1]]$start >= (ranges1[[1]]$start - max5diff) && ranges2[[1]]$start <= (ranges1[[1]]$start + max5diff))) {
	            return(FALSE)
        }

        # Check last coordinate with flexibility
           if(!(ranges2[[length(ranges2)]]$end >= (ranges1[[length(ranges1)]]$end - max3diff) && ranges2[[length(ranges2)]]$end <= (ranges1[[length(ranges1)]]$end + max3diff))) {
		        return(FALSE)

	   }
	
	# check if ranges have sufficient length
	   if(length(ranges1) < 2 || length(ranges2) <2 ) {
	  	return(FALSE)
	   }

        # Check intermediate coordinates for exact match
        for(i in 2:(length(ranges1) - 1)) {
		        if(!length(ranges1[[i]]) || !length(ranges2[[i]]) || ranges1[[i]]$start != ranges2[[i]]$start || ranges1[[i]]$end != ranges2[[i]]$end) {
				            return(FALSE)
	        }
	    }

	    return(TRUE)
}

match_ids_local <-  data.frame(
			gene_id = I(vector("list", 0)),
			transcript1 = I(vector("list", 0)),
			transcript2 = I(vector("list", 0)),
			concatanated_coord1 = I(vector("list", 0)),
			concatanated_coord2 = I(vector("list", 0)),
			numberOfExons = I(vector("list", 0)))




# Open the file and write the header
output_file_path <- paste0(output_dir, output_filename)  # Adjust the path accordingly

# Specify the file name and extension in the output path
output_file <- file(output_file_path, open = "wt")
write.table(x = data.frame(), file = output_file, sep = '\t', row.names = FALSE, col.names = TRUE)

for(each_coordinate in df_exons1_grouped$concatanated_coord) {
    transcript_id1 <- df_exons1_grouped[df_exons1_grouped$concatanated_coord == each_coordinate, 'transcript_id']
    gene_id1 <- isoform_gene_names1[isoform_gene_names1$PBid == transcript_id1$transcript_id, 'GeneName']

    if(length(gene_id1) == 1) {
        ranges1 <- make_a_range_all(each_coordinate)

    	if(length(ranges1) >= 1) {

        number_of_exons <- df_exons1_grouped[df_exons1_grouped$concatanated_coord == each_coordinate, 'NumberOfExons']

        transcript_ids <- isoform_gene_names2[isoform_gene_names2$GeneName == gene_id1, 'PBid']

        if(length(transcript_ids) >= 1) {
            transcript2_positions <- df_exons2_grouped[df_exons2_grouped$NumberOfExons == number_of_exons$NumberOfExons & df_exons2_grouped$transcript_id %in% transcript_ids, 'concatanated_coord']
            
            if(length(transcript2_positions$concatanated_coord) >= 1) {
                for(each_positions in 1:length(transcript2_positions$concatanated_coord)) {
                    coordinates_transcript2 <- transcript2_positions$concatanated_coord[each_positions]
                    ranges2 <- make_a_range_all(coordinates_transcript2)
                    transcript_id2 <- df_exons2_grouped[df_exons2_grouped$concatanated_coord == coordinates_transcript2, 'transcript_id']

                    result <- check_criteria(ranges1, ranges2)

                    if(isTRUE(unique(result))) {
                        match_ids <- data.frame(
                            gene_id = I(list(unique(gene_id1))),
                            transcript1 = I(list(transcript_id1$transcript_id)),
                            transcript2 = I(list(transcript_id2$transcript_id)),
                            concatanated_coord1 = I(list(unique(each_coordinate))),
                            concatanated_coord2 = I(list(unique(coordinates_transcript2))),
                            numberOfExons = I(list(unique(number_of_exons$NumberOfExons)))
                        )

		    	# Append the match_ids data frame to the file
		         write.table(match_ids, file = output_file, sep = '\t', row.names = FALSE, col.names = FALSE, append = TRUE)

                        #match_ids_local <- rbind(match_ids_local, match_ids)
		}    }
                }
            }
        }
    }
}



# close the file
close(output_file)

