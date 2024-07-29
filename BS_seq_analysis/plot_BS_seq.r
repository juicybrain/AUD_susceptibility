
library("tidyverse")


dat = read.table("Rgs_9_CG_call.txt", header=T)
library("ComplexHeatmap")
library("circlize")


P1 = Heatmap(dat)

dat_M = as.matrix(dat[,2:8])
row.names(dat_M) = dat$Group

svg('D1_DMR.center.noscale_clr12_04.24.centroid.svg',width = 8,height = 8)
h1 <- Heatmap(dat_M,name="mCG",
            row_km = 3,
            row_gap = unit(1, "mm"),
            column_gap = unit(1, "mm"),
            #row_km_repeats = 10,
            col=colorRamp2(c(0,1), c("white", "black")),
            clustering_distance_rows =  "euclidean",
            clustering_method_rows = "centroid",
            show_row_names = T,
            row_title = "%s",
            border = TRUE,
            show_row_dend = F,
            use_raster = TRUE,
            #column_split=Group,
            cluster_column_slices = FALSE,
            column_order=c("CG1","CG2","CG3","CG4","CG5","CG6","CG7")      
)
dev.off()

h2 <- Heatmap(dat_M,name="mCG",
         #   row_km = 3,
         cluster_rows = F,
         cluster_columns = F,
            row_gap = unit(1, "mm"),
            column_gap = unit(1, "mm"),
            #row_km_repeats = 10,
            col=colorRamp2(c(0,1), c("white", "black")),
           # clustering_distance_rows =  "euclidean",
           # clustering_method_rows = "centroid",
            show_row_names = T,
            row_title = "%s",
            border = TRUE,
            show_row_dend = F,
            use_raster = TRUE,
            #column_split=Group,
            cluster_column_slices = FALSE,
            column_order=c("CG1","CG2","CG3","CG4","CG5","CG6","CG7")      
)

#### ggplot
library("reshape2")
dat$spl = paste(dat$Group,rownames(dat), sep="_")
dat_long = melt(dat[,c(9,2,3,4,5,6,7,8)], id.vars = "spl")

dat_long$mC = ifelse(dat_long$value >0.5, "mC", "C" )

p_bs <- ggplot(dat_long, aes(x=variable, y=spl, shape=value, size= 5 )) +
geom_point(aes(fill=mC), color="black", pch=21, stroke=0.5 )+ scale_fill_manual(values=c("mC"="black", "C"="white")) + theme(panel.background = element_rect(fill = "white", colour = "grey50"))

ggsave("BS.svg",plot=p_bs, device="svg",width=8, height=16, units = "cm")


ggsave("BS.svg",plot=p_bs, device="svg",width=8, height=8, units = "cm")

### separate by Groups
dat_long_2 = melt(dat, id.vars = c("Group", "spl"))
dat_long_2$mC = ifelse(dat_long$value >0.5, "mC", "C" )


# Assuming dat_long has a column named 'group' that contains the group information
p_bs <- ggplot(dat_long_2, aes(x = variable, y = spl, shape = value, size = 8)) +
  geom_point(aes(fill = mC), color = "black", pch = 21, stroke = 0.5) +
  scale_fill_manual(values = c("mC" = "black", "C" = "white")) +
  facet_wrap(~ Group, nrow = 4, scales = "free_y") +  # Faceting by group, with free y scales
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        strip.background = element_blank(),  # Optional: styling for facet labels
        #strip.text = element_text(size = 12, face = "bold")
        strip.text = element_blank()
        )  # Optional: styling for facet text

ggsave("separated_plots.svg", p_bs, width = 3, height = 7)  # Save the plot as a PNG file



a=c(1/3, 0.3)
b=c(5/9, 0.4)

success_a <- 30; total_a <- 100  # 30% of Group A
success_b <- 45; total_b <- 150  # 30% of Group B

# Chi-squared test
test <- prop.test(c(success_a, success_b), c(total_a, total_b))

contingency_table <- matrix(c(success_a, total_a - success_a, success_b, total_b - success_b), nrow = 2)
fisher_test <- fisher.test(contingency_table)

print(test)
print(fisher_test)