library(coloc)
library(reticulate)
library(data.table)

data = py_load_object(snakemake@input[[1]])
num_tissues = dim(data$betas)[1]
num_pairs = num_tissues * (num_tissues -1) / 2 + num_tissues
table<-data.frame(
	"t1"=1:num_pairs,
	"t2"=1:num_pairs,
	"PPH0"=1:num_pairs,
	"PPH1"=1:num_pairs,
	"PPH2"=1:num_pairs,
	"PPH3"=1:num_pairs,
	"PPH4"=1:num_pairs
)

i <- 0
for (t1 in c(1:num_tissues)){
	for (t2 in c(t1:num_tissues)){c()
		i <- i + 1
		beta1 = data$betas[t1,]
		beta2 = data$betas[t2,]
		v1 = data$s2[t1,]
		v2 = data$s2[t2,]
		res = coloc.abf(
			dataset1=list(beta=beta1, varbeta=v1, N=838, sdY=sd(data$Y[t1,]), type='quant'),
			dataset2=list(beta=beta2, varbeta=v2, N=838, sdY=sd(data$Y[t2,]), type='quant')
		)
		table$t1[i] <- t1
		table$t2[i] <- t2
		table$PPH0[i] <- res$summary[2]
		table$PPH1[i] <- res$summary[3]
		table$PPH2[i] <- res$summary[4]
		table$PPH3[i] <- res$summary[5]
		table$PPH4[i] <- res$summary[6]
	}
}

fwrite(table, snakemake@output[[1]], quote=F, row.names=F, sep="\t")