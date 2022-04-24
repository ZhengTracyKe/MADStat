library('R.matlab')
library('igraph')
library("ggplot2")
library("ggrepel")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Q_results <- readMat(paste("output/DiversityResults_coau.mat",sep = ''),
                     sparseMatrixClass=c("Matrix"),encoding = "utf-8")


name <- unlist(Q_results$authorNames)
p_value <- Q_results$SgnQpvals
N_list <- Q_results$N.coauthors
ranking <- Q_results$ranking[1,]

sum(p_value < 0.05)


hist(p_value,breaks = 20,freq = T)

######## histogram of p vlaues #########
barfill <- "#4271AE"
barlines <- "#1F3552"

df <- data.frame(p_value = p_value)
p <- ggplot(df, aes(x = p_value)) +
  geom_histogram(aes(y = ..count..), binwidth = 0.05,
                 colour = barlines, fill = barfill, alpha=0.6) +
  scale_y_continuous(name = "Count") +
  scale_x_continuous(name = "") +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.line = element_line(size=1, colour = "black"),
        axis.text.x=element_text(colour="black", size = 20),
        axis.text.y=element_text(colour="black", size = 20))
  ggtitle("")
p

dev.copy2pdf(file=paste("p_Q_coau_2.pdf",sep = ''),width = 8.5,height = 6.5)

######## histogram of N #########
barfill <- "coral"
barlines <- "#1F3552"

df <- data.frame(size = N_list)
p <- ggplot(df, aes(x = size)) +
  geom_histogram(aes(y = ..count..), binwidth = 5,
                 colour = barlines, fill = barfill, alpha=0.6) +
  scale_y_continuous(name = "Count") +
  scale_x_continuous(name = "") +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.line = element_line(size=1, colour = "black"),
        axis.text.x=element_text(colour="black", size = 20),
        axis.text.y=element_text(colour="black", size = 20))
ggtitle("")
p


dev.copy2pdf(file=paste("N_coau_1.pdf",sep = ''),width = 8.5,height = 6.5)







name_ranked[1:10]
name_ranked <- name[ranking[1:1000]]
output <- data.frame(name = name_ranked,p_value = p_value,size = N_list)
write.csv(output, file = "q_coauthor.csv",
          row.names = FALSE, quote = TRUE, fileEncoding = "utf-8")


hist(log(p_value + 1e-8,base = 10),breaks = 30)
dev.copy2pdf(file=paste("log_p_Q_coau_2.pdf",sep = ''))


plot(N_list ~ p_value)
dev.copy2pdf(file=paste("N_p_Q_coau_2.pdf",sep = ''))






############## directed networks ##################
rm(list=ls())
Q_results <- readMat(paste("output/Q_stat_citee_citer.mat",sep = ''),
                     sparseMatrixClass=c("Matrix"))

A <- Q_results$A
name <- unlist(Q_results$name)
ranking <- Q_results$ranking
name_interest <- name[ranking[1:1000]]

# p_citer<- Q_results$Q.pv.list.citer
# p_citee <- Q_results$Q.pv.list.citee
N_citer <- Q_results$N.net.list.citer
N_citee <- Q_results$N.net.list.citee

name_ranked <- unlist(Q_results$name.ranked)[1:1000]
q_citee <- Q_results$Q.st.list.citee[,1]
q_citer <- Q_results$Q.st.list.citer[,1]

sum(is.na(q_citer))
sum(is.na(q_citee))

id_outlier <- (is.na(q_citee) | is.na(q_citer))
sum(id_outlier)

q_citee <- q_citee[!id_outlier]
q_citer <- q_citer[!id_outlier]

name_ranked <- name_ranked[!id_outlier]


sum(q_citer > q_citee)
sum(q_citee == q_citer)
sum(q_citee > q_citer)
########## draw the scatter plot of q_citer vs. q_citee ###########
to_english <- function(s) {
  # 1 character substitutions
  old1 <- "ŚŠŞşśšŽŻŻŻżžźþÅȦàáăặắâãäąåḃČĆÇçćčĐèééėễėêëęȩİİìıíĭîïķľłŁðñņńØÖòóôõöőùůūúûůüÜýÿğĝǧḡeřţt"
  new1 <- "SSSsssZZZZzzzyAAaaaaaaaaaabCCCcccDeeeeeeeeeeIIiiiiiikllLdnnnOOoooooouuuuuuuUyyggggertt"
  s1 <- chartr(old1, new1, s)
  
  # 2 character substitutions
  old2 <- c("œ", "ß", "æ", "ø","a̧","g̃","t́")
  new2 <- c("oe", "ss", "ae", "oe","a","g","t")
  s2 <- s1
  for(i in seq_along(old2)) s2 <- gsub(old2[i], new2[i], s2, fixed = TRUE)
  
  return(s2)
}


c <- iconv(name_ranked, from = "", to = "utf-8")
c <- iconv(c, from = "utf-8", to = "")
c <- to_english(c)

name_ranked1 <- c
select <- which(q_citee < 600)
q_citer1 <- q_citer[select]
q_citee1 <- q_citee[select]
name_ranked1 <- name_ranked1[select]

Q <- data.frame(q_citer = q_citer1,q_citee = q_citee1)
line_lim <- min(max(q_citer1),max(q_citee1))
p <- ggplot(data = Q, aes(x = q_citer, y = q_citee)) + theme_bw() +
  geom_point(shape = 16, colour = "grey", size = 2) +
  geom_line(data = data.frame(q_citee = c(0,line_lim),q_citer = c(0,line_lim)),
            linetype = 2, colour = 'black',size = 1)+
  theme(text = element_text(size=20))

m <- 10

label_index <- c()
top <- 10
show_n <- 5
d <- rowSums(A)

name_ranked1[1:10]

tmp <- setdiff(order(q_citee1 - q_citer1,decreasing = T),1:m)
label_index <- append(label_index,tmp[1:(2*show_n)])
tmp <- setdiff(order(q_citer1 - q_citee1,decreasing = T),1:m)
label_index <- append(label_index,tmp[1:show_n])

anchor_point <- data.frame(q_citer = q_citer1[label_index],q_citee = q_citee1[label_index])
high_point <- data.frame(q_citer = q_citer1[1:m],q_citee = q_citee1[1:m])

label_size <- c(rep(6,dim(anchor_point)[1]),
                rep(6,dim(high_point)[1]))
p1 <- p + geom_point(data = rbind(anchor_point,high_point),
                     shape = 16,
                     colour = c(rep("darkorange",dim(anchor_point)[1]),
                                rep("red",dim(high_point)[1])),
                     size = 2) + 
  geom_text_repel(data = rbind(anchor_point,high_point),
                  aes(label = c(name_ranked1[label_index],name_ranked1[1:m])),
                  colour = c(rep("black",dim(anchor_point)[1]),
                             rep("red",dim(high_point)[1])),
                  size = label_size, 
                  box.padding = unit(0.2, "lines"))

p1
dev.copy2pdf(file=paste("q_citer_citee.pdf",sep = ''))