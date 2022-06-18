options(warn=-1)
suppressWarnings(suppressMessages(library("Seurat", quietly = TRUE)))
suppressWarnings(suppressMessages(library("Epi", quietly = TRUE)))
suppressWarnings(suppressMessages(library("dplyr", quietly = TRUE)))
suppressWarnings(suppressMessages(library("reshape2", quietly = TRUE)))
suppressWarnings(suppressMessages(library("pscl", quietly = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
path <- args[1] #"/Volumes/archive/brownlab/davidc/kennylab_copy/ACME.rds"
str_group_1 <- args[2] #"0,1,4,6,12,17,25,33"
str_group_2 <- args[3] #"16"
fileout <- args[4] #"ACME_DE_4.txt"

obj <- readRDS(path)

group_1 <- as.numeric(strsplit(str_group_1, ",")[[1]])
group_2 <- as.numeric(strsplit(str_group_2, ",")[[1]])
cells_1 <- WhichCells(obj,idents=group_1)
cells_2 <- WhichCells(obj,idents=group_2)

data_count <- obj[["RNA"]]@counts
genes <- rownames(data_count)
ngene <- length(genes)

count_table <- function(v, tag) {
	counts <- as.numeric(table(cut(v, c(0:201,Inf), right=F)))
	labels <- as.character(0:201)
	labels <- paste(tag, "_", labels, sep="")
	labels[202] <- paste(labels[202], "_Plus", sep="")
	df_count <- as.data.frame(as.list(counts))
	df_sum <- data.frame(X=sum(counts))
	names(df_count) <- labels
	names(df_sum)[1] <- paste(tag, "_count_total", sep="")
	df_count <- cbind(df_count, df_sum)
	return(df_count)
}

base_case <- function() {
	df_out <- data.frame(A=0,B=0,C=0,D=1)
	names(df_out) <- c("Effect_size", "CI95_low", "CI95_high", "P_value")
	df_stat <- data.frame(Scenario=0,
	                      Mean_case_count_gene=0,
						  Mean_control_count_gene=0)
	df_out <- cbind(df_stat,df_out)
	return(df_out)
}

zip_one_group <- function(data, v) {
	nn <- length(data)
	offs <- log(1/nn)
	mean_count <- mean(data)
	df <- data.frame(X=1,Y=data,W=offs)
	
	ml <- zeroinfl(Y~1+offset(W)|1,data=df)
	effs <- coefficients(ml)[1:1]
	df_out <- as.data.frame(ci.lin(ml))[1:1,c(1,5,6,4)]
	names(df_out) <- c("Effect_size", "CI95_low", "CI95_high", "P_value")
	df_out$Effect_size <- effs
	df_out$CI95_low[!is.finite(df_out$CI95_low)] <- -Inf
	df_out$CI95_high[!is.finite(df_out$CI95_high)] <- Inf
	df_out$P_value[!is.finite(df_out$P_value)] <- 1
	if (v<0) {
		df_out$Effect_size <- sign(v) * df_out$Effect_size
		df_out$CI95_low <- sign(v) * df_out$CI95_low
		df_out$CI95_high <- sign(v) * df_out$CI95_high
		tmp <- df_out$CI95_high
		df_out$CI95_high <- df_out$CI95_low
		df_out$CI95_low <- tmp
	}
	df_stat <- data.frame(Scenario=1,
	                      Mean_case_count_gene=ifelse(v>0,mean_count,0),
						  Mean_control_count_gene=ifelse(v>0,0,mean_count))
	df_out <- cbind(df_stat,df_out)
	return(df_out)
}

poisson_one_group <- function(data, v) {
	nn <- length(data)
	offs <- log(1/nn)
	mean_count <- mean(data)
	df <- data.frame(X=1,Y=data,W=offs)
	
	ml <- glm(Y~1+offset(W),data=df,family="poisson")
	effs <- coefficients(ml)[1:1]
	df_out <- as.data.frame(ci.lin(ml))[1:1,c(1,5,6,4)]
	df_out$Effect_size <- effs
	df_out$CI95_low[!is.finite(df_out$CI95_low)] <- -Inf
	df_out$CI95_high[!is.finite(df_out$CI95_high)] <- Inf
	df_out$P_value[!is.finite(df_out$P_value)] <- 1
	names(df_out) <- c("Effect_size", "CI95_low", "CI95_high", "P_value")
	if (v<0) {
		df_out$Effect_size <- sign(v) * df_out$Effect_size
		df_out$CI95_low <- sign(v) * df_out$CI95_low
		df_out$CI95_high <- sign(v) * df_out$CI95_high
		tmp <- df_out$CI95_high
		df_out$CI95_high <- df_out$CI95_low
		df_out$CI95_low <- tmp
	}
	df_stat <- data.frame(Scenario=2,
	                      Mean_case_count_gene=ifelse(v>0,mean_count,0),
						  Mean_control_count_gene=ifelse(v>0,0,mean_count))
	df_out <- cbind(df_stat,df_out)
	return(df_out)
}

poisson_two_group <- function(data1, data2) {
	gene_count_case <- data1
	gene_count_ctrl <- data2
	mean_count_gene_case <- mean(gene_count_case)
	mean_count_gene_ctrl <- mean(gene_count_ctrl)
	df_case <- data.frame(X=1,Y=gene_count_case)
	df_ctrl <- data.frame(X=0,Y=gene_count_ctrl)
	df_combine <- rbind(df_case, df_ctrl)
	
	ml <- glm(Y~X,data=df_combine,family="poisson")
	effs <- coefficients(ml)[2:2]
	df_out <- as.data.frame(ci.lin(ml))[2:2,c(1,5,6,4)]
	names(df_out) <- c("Effect_size", "CI95_low", "CI95_high", "P_value")
	df_out$Effect_size <- effs
	df_out$CI95_low[!is.finite(df_out$CI95_low)] <- -Inf
	df_out$CI95_high[!is.finite(df_out$CI95_high)] <- Inf
	df_out$P_value[!is.finite(df_out$P_value)] <- 1
	df_stat <- data.frame(Scenario=3,
	                      Mean_case_count_gene=mean_count_gene_case,
						  Mean_control_count_gene=mean_count_gene_ctrl)
	df_out <- cbind(df_stat,df_out)
	return(df_out)
}

zip_two_group <- function(data1, data2) {
	gene_count_case <- data1
	gene_count_ctrl <- data2
	mean_count_gene_case <- mean(gene_count_case)
	mean_count_gene_ctrl <- mean(gene_count_ctrl)
	df_case <- data.frame(X=1,Y=gene_count_case)
	df_ctrl <- data.frame(X=0,Y=gene_count_ctrl)
	df_combine <- rbind(df_case, df_ctrl)
	
	ml <- zeroinfl(Y~X|1,data=df_combine)
	effs <- coefficients(ml)[2:2]
	df_out <- as.data.frame(ci.lin(ml))[2:2,c(1,5,6,4)]
	names(df_out) <- c("Effect_size", "CI95_low", "CI95_high", "P_value")
	df_out$Effect_size <- effs
	df_out$CI95_low[!is.finite(df_out$CI95_low)] <- -Inf
	df_out$CI95_high[!is.finite(df_out$CI95_high)] <- Inf
	df_out$P_value[!is.finite(df_out$P_value)] <- 1
	df_stat <- data.frame(Scenario=4,
	                      Mean_case_count_gene=mean_count_gene_case,
						  Mean_control_count_gene=mean_count_gene_ctrl)
	df_out <- cbind(df_stat,df_out)
	return(df_out)
}

get_results <- function(x) {
	gene_count_case <- as.numeric(data_count[x,cells_1])
	gene_count_ctrl <- as.numeric(data_count[x,cells_2])
	b_case_all_zero <- all(gene_count_case == 0)
	b_ctrl_all_zero <- all(gene_count_ctrl == 0) 
	b_case_all_non_zero <- all(gene_count_case != 0)
	b_ctrl_all_non_zero <- all(gene_count_ctrl != 0)
	
	df_out <- data.frame()
	if (b_case_all_zero && b_ctrl_all_zero) {
		df_out <- base_case()
	} else if (b_case_all_zero && (! b_ctrl_all_zero)) {
		if (b_ctrl_all_non_zero) {
			df_out <- poisson_one_group(gene_count_ctrl,-1)
		} else {
			df_out <- zip_one_group(gene_count_ctrl,-1)
		}
	} else if ((! b_case_all_zero) && b_ctrl_all_zero) {
		if (b_case_all_non_zero) {
			df_out <- poisson_one_group(gene_count_case,1)
		} else {
			df_out <- zip_one_group(gene_count_case,1)
		}
	} else if (b_case_all_non_zero && b_ctrl_all_non_zero) {
		df_out <- poisson_two_group(gene_count_case,
		                            gene_count_ctrl)
	} else {
		df_out <- zip_two_group(gene_count_case,
		                        gene_count_ctrl)
	}
	
	wilcox <- wilcox.test(gene_count_case, gene_count_ctrl)
	df_wilcox <- data.frame(Wilcox_P_value=as.numeric(wilcox$p.value)[1])
	
	df_base <- data.frame(Gene=x)
	df_gene_count_tbl_case <- count_table(gene_count_case, "Group_1")
	df_gene_count_tbl_ctrl <- count_table(gene_count_ctrl, "Group_2")
	df_out <- cbind(df_base,df_gene_count_tbl_case,df_gene_count_tbl_ctrl,df_out,df_wilcox)
	return(df_out)
}

gene <- genes[1]
df <- get_results(gene)
write.table(df, file=fileout, 
                row.names=FALSE, 
                quote=FALSE, 
                sep="\t")

for (i in 2:ngene) {
	gene <- genes[i]
	df <- get_results(gene)
	write.table(df, file=fileout, 
	                row.names=FALSE, 
	                col.names=FALSE, 
	                quote=FALSE, 
	                sep="\t", 
	                append=TRUE)
}
