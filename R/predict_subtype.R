
predict.molecular.subtype <- function(expression.matrix)
{
    c1.genes <- intersect(rownames(expression.matrix),sig.genes$C1$symbol)
    C1.exp <- expression.matrix[c1.genes,]
    C1 <- sig.genes$C1
    C1.sel <- C1[C1$symbol %in% c1.genes,]
    rownames(C1.sel) <- C1.sel$symbol
    C1.ord <- C1.sel[rownames(C1.exp),]

    c2.genes <- intersect(rownames(expression.matrix),sig.genes$C2$symbol)
    C2.exp <- expression.matrix[c2.genes,]
    C2 <- sig.genes$C2
    C2.sel <- C2[C2$symbol %in% c2.genes,]
    rownames(C2.sel) <- C2.sel$symbol
    C2.ord <- C2.sel[rownames(C2.exp),]

    c4.genes <- intersect(rownames(expression.matrix),sig.genes$C4$symbol)
    C4.exp <- expression.matrix[c4.genes,]
    C4 <- sig.genes$C4
    C4.sel <- C4[C4$symbol %in% c4.genes,]
    rownames(C4.sel) <- C4.sel$symbol
    C4.ord <- C4.sel[rownames(C4.exp),]

    c5.genes <- intersect(rownames(expression.matrix),sig.genes$C5$symbol)
    C5.exp <- expression.matrix[c5.genes,]
    C5 <- sig.genes$C5
    C5.sel <- C5[C5$symbol %in% c5.genes,]
    rownames(C5.sel) <- C5.sel$symbol
    C5.ord <- C5.sel[rownames(C5.exp),]

    C1.vals <- C1.exp * C1.ord$logFC
    C2.vals <- C2.exp * C2.ord$logFC
    C4.vals <- C4.exp * C4.ord$logFC
    C5.vals <- C5.exp * C5.ord$logFC

    C1.scores <- apply(C1.vals,2,sum)
    C2.scores <- apply(C2.vals,2,sum)
    C4.scores <- apply(C4.vals,2,sum)
    C5.scores <- apply(C5.vals,2,sum)

    scores <- data.frame(C1 = C1.scores, C2 = C2.scores, C4 = C4.scores, C5 = C5.scores)
    subtype.scores <- scale(scores)


    old.rownames <- rownames(subtype.scores)
    subtype.scores <- apply(subtype.scores, 2, scale)
    rownames(subtype.scores) <- old.rownames
    ind <- apply(subtype.scores, 1, which.max)
    subclasses <- colnames(subtype.scores)[ind]
    subclasses <- factor(subclasses, levels = c("C2", "C4", "C5", 
        "C1"))
    return(list(subtypes = subclasses, subtype.scores = subtype.scores))
}



