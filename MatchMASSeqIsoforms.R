# Author: Tülay Karakulak
# Created: 29.11.2023
# Description: The script takes two gff files which are the output of MAS-Seq Iso-Seq workflow and matches the same transcripts based on their exon coordinates with some flexibility.


library(rtracklayer)
library(dplyr)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)

# Check if the correct number of arguments (5) is provided
if(length(args) != 7) {
	  stop("Please provide the paths for two input GFF files, 2 Seurat genes.tsv, one output file, number of nucleotides flexibility for 5 prime and 3prime")
}

# Assigning input and output file paths from the arguments
input_path1 <- args[1]
input_path2 <- args[2]
genes_path1 <- args[3]
genes_path2 <- args[4]
output_path <- args[5]
max5diff <-  as.numeric(args[6]) # Max nucleotide length at the 5' end
max3diff <- as.numeric(args[7])  # Max nucleotide length at the 3' end


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
make_a_range <- function(coordinates) {
		
			start_pos <- strsplit(strsplit(coordinates,  ",")[[1]], '-')[[1]][1]
			end_pos <- strsplit(strsplit(coordinates,  ",")[[length(coordinates)]], '-')[[1]][2]        
			df <- data.frame(start_pos = as.integer(start_pos), end_pos = as.integer(end_pos))
			
			return(df)
		}

# Function to check if two sets of coordinates meet the defined criteria - matching transcript based on max3diff and max5diff
check_criteria <- function(df1, df2) {
			  
			start_pos_check <- df2$start_pos >= (df1$start_pos - max5diff) & df2$start_pos <= (df1$start_pos + max5diff)
			end_pos_check <- df2$end_pos >= (df1$end_pos - max3diff) & df2$end_pos <= (df1$end_pos + max3diff)
			    
			return(start_pos_check & end_pos_check)
		    }

# Initializing a data frame to store matched IDs
match_ids_local <- data.frame(
			gene_id = I(vector("list", 0)),
			transcript1 = I(vector("list", 0)),
			transcript2 = I(vector("list", 0)),
			concatanated_coord1 = I(vector("list", 0)),
			concatanated_coord2 = I(vector("list", 0)),
			numberOfExons = I(vector("list", 0)))


# Main loop to process and match coordinates between the two datasets
for(each_coordinate in df_exons1_grouped$concatanated_coord) {
			    	
	transcript_id1 <- df_exons1_grouped[df_exons1_grouped$concatanated_coord == each_coordinate, 'transcript_id']
	gene_id1 <- isoform_gene_names1[isoform_gene_names1$PBid == transcript_id1$transcript_id, 'GeneName']

	if(length(gene_id1) == 1) {
				  
			coordinates <- strsplit(each_coordinate,  ",")[[1]]
		    start_end_pos <- make_a_range(c(coordinates[1], coordinates[length(coordinates)]))
			number_of_exons <- df_exons1_grouped[df_exons1_grouped$concatanated_coord == each_coordinate, 'NumberOfExons']
				
			  # take the coordinates of the transcripts which has equal number of exons

			transcript_ids <- isoform_gene_names2[isoform_gene_names2$GeneName == gene_id1, 'PBid']

			if(length(transcript_ids) >= 1) {
			transcript2_positions <- df_exons2_grouped[df_exons2_grouped$NumberOfExons == number_of_exons$NumberOfExons & df_exons2_grouped$transcript_id %in% transcript_ids, 'concatanated_coord']
			    

			if(length(transcript2_positions$concatanated_coord) >= 1){

			for(each_positions in 1:length(transcript2_positions$concatanated_coord)) {
				      
				coordinates_transcript2 <- transcript2_positions$concatanated_coord[each_positions]
			    coordinates2 <- strsplit(coordinates_transcript2,  ",")[[1]]
				start_end_pos_transcript2 <- make_a_range(c(coordinates2[1], coordinates2[length(coordinates2)]))
				transcript_id2 <- df_exons2_grouped[df_exons2_grouped$concatanated_coord == coordinates_transcript2, 'transcript_id']

				result <- check_criteria(start_end_pos, start_end_pos_transcript2)

				if(isTRUE(unique(result))) {
					match_ids <- data.frame(
					gene_id = I(list(unique(gene_id1))),
					transcript1 = I(list(transcript_id1$transcript_id)),
					transcript2 = I(list(transcript_id2$transcript_id)),
					concatanated_coord1 = I(list(unique(each_coordinate))),
					concatanated_coord2 = I(list(unique(coordinates_transcript2))),
					numberOfExons = I(list(unique(number_of_exons$NumberOfExons))))
					match_ids_local <- rbind(match_ids_local, match_ids)
				}
				}
			    }
		    }}}

match_ids_final <- match_ids_local %>% dplyr::distinct()
# Combine all the data frames in the results list to a single data frame
write.table(match_ids_final, file = output_path, row.names=FALSE, sep='\t')
saveRDS(match_ids_final, file = paste0(output_path, ".RDS"))

