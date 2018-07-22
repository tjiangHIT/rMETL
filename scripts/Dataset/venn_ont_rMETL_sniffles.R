
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
area1 = 3287,
area2 = 63925,
cross.area = 2912,
category = c("A", "B"),
fill = c("blue", "red"),
cat.col = c("blue", "red"),
# lty = "blank",
cex = 2,
cat.cex = 2,
# cat.pos = c(285, 105),
# cat.dist = 0.09,
# cat.just = list(c(-1, -1), c(1, 1)),
ext.pos = c(30, 0),
# ext.dist = c(-0.05, 0),
# ext.length = 0.85,
# ext.line.lwd = 2,
# # ext.line.lty = "dashed"
);

# Writing to file
tiff(filename = "nanopore.tiff", compression = "lzw");
grid.draw(venn.plot);
dev.off();