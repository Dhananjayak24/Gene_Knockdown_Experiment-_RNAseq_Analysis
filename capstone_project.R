#Description: You're given RNAseq-based FPKM gene expression matrix from gene knockdown 
#experiment. The experiment was designed as follows. 
#There are 4 conditions: un-treated cells (columns with prefix UT), 
#negative control cells (columns with prefix NC), 
#cells with knockdown method A (columns with prefix A), 
#and cells with knockdown method C (columns with prefix 
#C). There are also two time points: an early time point (columns with suffix .1) 
#and late time point (columns with suffix .2). #This means that gene expression was 
#measured at two time points following the start of the experiment. 
#UT cells have not been treated at all, and this means that those cells have just 
#been growing in their growth culture medium. On the other hand NC cells 
#have been treated just like in conditions A and C with the only exception that it had 
#missing the single component necessary to undertake the gene knockdown. Finally, 
#cells under A and C were treated for the depletion of the same gene with two different methods.

#Install packages
install.packages("rankProd")
install.packages("VennDiagram")

#Load packages
library(RankProd)
library(VennDiagram)

#setup working directory
setwd("D:/R/udemy_course")

#Preparing data
cap_data = read.csv("capstone.csv")
gene_name = cap_data$gene_name
cap_data$X = NULL
cap_data$gene_name = NULL
#Here I removed column "X" and gene neames because they are not required in the matrix

####################################################################

#Scenario 1: Comparing Untreated and Negetive control at time point 1

#getting column names of related columns
col_names1 = names(cap_data)[1:6]
#creating subset of original data
data_1 = cap_data[,col_names1]
cl = c(rep(0,3),rep(1,3)) #This is to create a variable to note untreated and control

#We will use rankProducts package to list DEG list

result_1 = RankProducts(data = data_1,
                        logged = F,
                        cl = cl,
                        rand = 999) #here cl is class and rand is used to give seed number so that result is reproducible

# now we will create top table to list the genes 
top_table1 = topGene(result_1,
                     cutoff = 0.05,
                     method = "pval",
                     gene.names = gene_name,
                     logged = F)
#INTERPRETATION: In table 1 see values of Fold Change (FC)
#They are less that 1 that means denominator(Negetive control) is larger than numerator(Untreated), hence UPREGULATED
#The same if you see table2 FC values are higher so NC is DOWN REGULATED
#RESULT:
#We have 4849 genes UPREGULATED (table 1)
#4996 genes DOWN REGULATED

######################################################################################
#Scenario 2: Compare NC vs A and NC vs C at time point 1: 
#can you identify differentially expressed genes between control and knockdown? 
#Are the results reproducible with the two different depletion methods? 
#If yes would you have more confidence in the impact of gene depletion for those genes?

col_names2A = c("NC1.1","NC2.1","NC3.1","A1.1","A2.1","A3.1")
col_names2C = c("NC1.1","NC2.1","NC3.1","C1.1","C2.1","C3.1")

data2A = cap_data[,col_names2A]
data2C = cap_data[,col_names2C]

#We will use rankProducts package to list DEG list
result_2A = RankProducts(data = data2A,
                        logged = F,
                        cl = cl,
                        rand = 999) #here cl is class and rand is used to give seed number so that result is reproducible

# now we will create top table to list the genes 
top_table2A = topGene(result_2A,
                     cutoff = 0.05,
                     method = "pval",
                     gene.names = gene_name,
                     logged = F)
#similarly for another set
#We will use rankProducts package to list DEG list
result_2C = RankProducts(data = data2C,
                         logged = F,
                         cl = cl,
                         rand = 999) #here cl is class and rand is used to give seed number so that result is reproducible

# now we will create top table to list the genes 
top_table2C = topGene(result_2C,
                     cutoff = 0.05,
                     method = "pval",
                     gene.names = gene_name,
                     logged = F)

diff_2A = c(rownames(top_table2A$Table1),rownames(top_table2A$Table2))
diff_2C = c(rownames(top_table2C$Table1),rownames(top_table2C$Table2))

#To answer the questions we need to draw Venn Diagram
#Either we can use package in Rstudio or use online tool https://bioinfogp.cnb.csic.es/tools/venny/

venn.diagram(
  x = list(diff_2A,diff_2C),
  category.names = c("differential_A","differential_C"),
  filename = "diff_A_C_venn.png",
  output = TRUE,
  fill = c("blue","red")
)
#Now the file in png format will be saved to your working directory
#ANSWERS:
#Yes, we can identify differentially expressed genes between control and knockdown
#Are the results reproducible with the two different depletion methods? 
#Parially yes because for two different knockdown methods there are almost equal number of genes differentially expressed(Venn diagram)

#If yes would you have more confidence in the impact of gene depletion for those genes?
#Answer : 
A_ratio = 3962/(4494+3962)
C_ratio = 3962/(3388+3962)
A_ratio
C_ratio
#We can see the ratio of similar genes expressed is almost same (0.46 and 0.54)
#So we can conclude the two methods expressions are almost same.

#####################################################################################
#QUESTION 3:
#Repeat the above for time point 2: are the effects more pronounced at the later time point? 
#Can you see higher reproducibility of gene expression patters in NC vs A and NC vs C  comparisons?

#We will copy same set of codes as above with some changes
col_names3A = c("NC1.2","NC2.2","NC3.2","A1.2","A2.2","A3.2")
col_names3C = c("NC1.1","NC2.2","NC3.2","C1.2","C2.2","C3.2")

data3A = cap_data[,col_names3A]
data3C = cap_data[,col_names3C]

#We will use rankProducts package to list DEG list
result_3A = RankProducts(data = data3A,
                         logged = F,
                         cl = cl,
                         rand = 999) #here cl is class and rand is used to give seed number so that result is reproducible

# now we will create top table to list the genes 
top_table3A = topGene(result_3A,
                      cutoff = 0.05,
                      method = "pval",
                      gene.names = gene_name,
                      logged = F)
#similarly for another set
#We will use rankProducts package to list DEG list
result_3C = RankProducts(data = data3C,
                         logged = F,
                         cl = cl,
                         rand = 999) #here cl is class and rand is used to give seed number so that result is reproducible

# now we will create top table to list the genes 
top_table3C = topGene(result_3C,
                      cutoff = 0.05,
                      method = "pval",
                      gene.names = gene_name,
                      logged = F)

diff_3A = c(rownames(top_table3A$Table1),rownames(top_table3A$Table2))
diff_3C = c(rownames(top_table3C$Table1),rownames(top_table3C$Table2))

#To answer the questions we need to draw Venn Diagram
#Either we can use package in Rstudio or use online tool https://bioinfogp.cnb.csic.es/tools/venny/

venn.diagram(
  x = list(diff_3A,diff_3C),
  category.names = c("differential_A","differential_C"),
  filename = "diff_A_C_venn_difftime.png",
  output = TRUE,
  fill = c("blue","red")
)
#calculating ratios
A_ratio_2 = 5107/(6264+5107)
C_ratio_2 = 5107/(3807+5107)
A_ratio_2
C_ratio_2

#Conclusion: There is higher intersection in second time data

#####################################################################################
#QUESTION 4:
#Is there an overlap in genes the appear differentially expressed by virtue of mere treatment and by virtue of actual knockdown? 
#If yes, do you think we can attribute changes in expression of those genes to the knockdown?

#Here we will compare for time point 1 data
#We will use output from question 1 and Question 2
#Que1 has control expression data and Que2 has expression at two methods
exp_genes = c(rownames(top_table1$Table1),rownames(top_table1$Table2))
knockdown_genes = c(rownames(top_table2A$Table1),rownames(top_table2A$Table2),
             rownames(top_table2C$Table1),rownames(top_table2C$Table2))
venn.diagram(
  x = list(exp_genes,knockdown_genes),
  category.names = c("UT vs NC","A&C vs NC"),
  filename = "UT_vs_A&C.png",
  output = TRUE,
  fill = c("blue","red")
)
#genes.

# There is a chunk of the genes that are differentially expressed in the negative control
# versus C experiment are also being differentially expressed in the untreated versus negative control.
#So that suggests there is some, genes that are possibly linked more to the treatment than the knockdown, 
#But the genes Differentially expressed by method of depletion are greater!!!

##############################################################################
#QUESTION : 5
#Knockdown has been undertaken for gene INTS12.
#Has this gene been successfully silenced? 
#Is there difference between time point 1 and time point 2?

#we will check whether gene INTS12 is present in different conditions

"INTS12" %in% diff_2A
"INTS12" %in% diff_2C
"INTS12" %in% diff_3A
"INTS12" %in% diff_3C

#So the answer is "true" for all the four consitions. hence it is succesfully silenced.

#Check difference in time point 1 and 2

View(top_table2A$Table2["INTS12",])
View(top_table3A$Table2["INTS12",])

#we check in Table2 because we are checking for DOWNREGULATION
#The FC value changed from 3.9 to 8.9 in time point 2 that means it is more silenced.

#####################################################################################
#QUESTON 6: Compare each condition with itself looking at differences between time point 1 and time point 2. 
#Can we say that some genes change in expression over time?

A_common = intersect(rownames(diff2A),rownames(diff3A))
#Here we are selecting only intersects for method A in both time point(after and before)


diff2A = rbind(top_table2A$Table1,top_table2A$Table2)
#we are combining upregulation and downregulation of time point 1 for method A
A_common_time1 = diff2A[A_common,]

diff3A = rbind(top_table3A$Table1,top_table3A$Table2)
#we are combining upregulation and downregulation of time point 2 for method A
A_common_time2 = diff3A[A_common,]


View(A_common_time1)
View(A_common_time2)

#Select down-regulated in KnockDown
select_down = intersect(rownames(A_common_time1[A_common_time1[,3]>1,]),
                   rownames(A_common_time2[A_common_time2[,3]>1,]))

#Number of genes that were more strongly downregulated in knockdown
sum(A_common_time1[select_down,3] < A_common_time2[select_down,3])

#Fraction of genes that were more strongly downregulated in KD (among downregulated both the times)
sum((A_common_time1[select_down,3] < A_common_time2[select_down,3])/length(select_down))
#So 81% of the genes are more strongly Down regulated in KD method

#Up-regulation

select_up = intersect(rownames(A_common_time1[A_common_time1[,3] < 1,]),
                        rownames(A_common_time2[A_common_time2[,3] < 1,]))

#Number of genes that were more strongly up-regulated in knockdown
sum(A_common_time1[select_up,3] > A_common_time2[select_up,3])

#Fraction of genes that were more strongly up-regulated in KD (among up-regulated both the times)
sum((A_common_time1[select_up,3] > A_common_time2[select_up,3])/length(select_up))
#So 77% of the genes are more strongly up regulated in KD method
