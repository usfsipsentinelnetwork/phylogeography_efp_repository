library('progress')
cluster <-
  function(data)
  {
    n <- 1
    for (i in 1 : (length(data) - 1))
    {
      if (data[i] != data[i + 1])
      {
        n[i + 1] <- n[i] + 1
      }
      else
      {
        n[i + 1] <- n[i]
      }
    }
    return(n)
  }
ancestor <-
function(tree, tip1, tip2)
{
  ance1 <- tree$edge[which(tree$edge[, 2] == tip1), 1]
  ance2 <- tree$edge[which(tree$edge[, 2] == tip2), 1]
  mi <- min(tree$edge[, 1])
  for (i in 1 : (tree$Nnode - 1))
  {
    if (ance1[i] > mi)  ance1[i + 1] <- tree$edge[which(tree$edge[,2] == ance1[i]), 1]
    else  break()
  }
  for (j in 1 : (tree$Nnode - 1))
  {
    if (ance2[j] > mi)  ance2[j + 1] <- tree$edge[which(tree$edge[,2] == ance2[j]), 1]
    else  break()
  }
  node <- max(intersect(ance1, ance2))
  return(node)
}


phylo.maker.verbose<-
function (sp.list, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL,
    output.sp.list = TRUE, output.tree = FALSE, scenarios = "S3",
    r = 1)
{
    treeX <- tree
    if (is.null(tree$node.label))
        tree$node.label <- rep("", tree$Nnode)
    dimnames(sp.list)[[2]][1:3] <- c("species", "genus", "family")
    print("applying list as factors")
    sp.list[sapply(sp.list, is.factor)] <- lapply(sp.list[sapply(sp.list,
        is.factor)], as.character)
    if (any(duplicated(sp.list$species))) {
        print("Duplicated species detected and removed.")
        print(sp.list$species[duplicated(sp.list$species)])
    }
    sp.list <- sp.list[!duplicated(sp.list$species), ]
    sp.list.original <- sp.list
    sp.list$species <- gsub(" ", "_", sp.list$species)
    sp.list$species <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$species,
        perl = TRUE)
    sp.list$genus <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$genus,
        perl = TRUE)
    sp.list$family <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$family,
        perl = TRUE)
    rnN <- data.frame(node.label = paste("N", 1:length(tree$node.label),
        sep = ""), oriN = tree$node.label, stringsAsFactors = FALSE)
    print("applying nodes as characters")
    nodes[, c("level", "family", "genus", "rn", "bn", "taxa")] <- lapply(nodes[,
        c("level", "family", "genus", "rn", "bn", "taxa")], as.character)
    tree$node.label <- paste("N", 1:length(tree$node.label),
        sep = "")
    kk <- c()
    print("Tip Labels (1/6)")
    pb1<- progress_bar$new(total = length(tree$tip.label))
    for (i in 1:length(tree$tip.label)) {
        kk <- c(kk, substring(tree$tip.label[i], 1, gregexpr("_",
            tree$tip.label[i])[[1]][1] - 1))
        pb1$tick()
    }
    m <- data.frame(num = 1:length(kk), genus = kk, species = tree$tip.label)
    m <- merge(m, nodes[, c("genus", "family")])
    mX <- m
    m <- m[, c("genus", "family")]
    m <- m[!duplicated(m$genus), ]
    dimnames(m)[[2]][2] <- "family_in_tree"
    m <- m[, c("genus", "family_in_tree")]
    m0 <- sp.list[!duplicated(sp.list$genus), c("genus", "family")]
    dimnames(m0)[[2]][2] <- "family_in_sp.list"
    mm <- merge(m0, m)
    g <- mm[which(is.na(match(paste(mm$genus, mm$family_in_sp.list,
        sep = "_"), paste(mm$genus, mm$family_in_tree, sep = "_")))),
        ]
    if (dim(g)[1] > 0) {
        print("Taxonomic classification not consistent between sp.list and tree.")
        print(g)
    }
    add.tip <- sp.list[which(is.na(match(sp.list$species, tree$tip.label))),
        ]
    status <- rep("prune", dim(sp.list)[1])
    status[which(is.na(match(sp.list$species, tree$tip.label)))] <- "bind"
    if (dim(add.tip)[1] == 0 & length(na.omit(match(sp.list$species,
        tree$tip.label))) == 0)
        stop("Incorrect format of species list.")
    if (length(setdiff(sp.list$species, treeX$tip.label)) ==
        0 & length(na.omit(match(sp.list$species, tree$tip.label))) >
        0) {
        print("All species in sp.list are present on tree.")
        splis <- sp.list.original
        treeX <- drop.tip(treeX, setdiff(treeX$tip.label, sp.list$species))
        splis$status <- "prune"
        phyloX <- list(scenario.1 = NULL, scenario.2 = NULL,
            scenario.3 = NULL, species.list = splis)
        if ("S1" %in% scenarios) {
            phyloX$scenario.1 <- treeX
        }
        if ("S2" %in% scenarios) {
            phyloX$scenario.2 <- treeX
        }
        if ("S3" %in% scenarios) {
            phyloX$scenario.3 <- treeX
        }
        phyloX[sapply(phyloX, is.null)] <- NULL
        return(phyloX)
        stop()
    }
    add.tip$sort <- ""
    add.tip$sort[which(!is.na(match(add.tip$genus, nodes[nodes$level ==
        "G", ]$genus)))] <- "G1"
    add.tip$sort[which(is.na(match(add.tip$genus, nodes[nodes$level ==
        "G", ]$genus)) & !is.na(match(add.tip$family, nodes[nodes$level ==
        "F", ]$family)))] <- "F1"
    add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort ==
        "F1", ]$genus)] <- "F2"
    a <- which(add.tip$sort == "")
    if (length(a) > 0) {
        print(paste("Note:", length(a), "taxa fail to be binded to the tree,",
            sep = " "))
        print(add.tip$species[a])
        status[match(add.tip$species[a], sp.list$species)] <- "fail to bind"
    }
    sp.list.original$status <- status
	
############ SCENARIO 1

    if ("S1" %in% scenarios) {
        t1 <- tree
        rnN1 <- rnN
        nG <- nodes[nodes$level == "G", ]
        nF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1" | add.tip$sort ==
            "F2", ]
		print("Scenario 1")	
		print("Checking F1 & F2")
		
        if (dim(data)[1] > 0) {
            print("Families (2/6)")
            pb2 <- progress_bar$new(total=dim(data)[1])
            for (i in 1:dim(data)[1]) {			
                num <- nF$bn[match(data$family[i], nF$family)]
				
				#print(paste(data$species[i], "(", data$family[i], "), num=", num))

                t1 <- at.node(t1, location.node = num, tip.label = data$species[i])
                pb2$tick()
            }
        }
        data <- add.tip[add.tip$sort == "G1", ]
        if (dim(data)[1] > 0) {
            print("Genera (3/6)")
            pb3 <- progress_bar$new(total=dim(data)[1])
            for (i in 1:dim(data)[1]) {
                num <- nG$bn[match(data$genus[i], nG$genus)]
                t1 <- at.node(t1, location.node = num, tip.label = data$species[i])
                pb3$tick()
            }
        }
        t1$edge.length <- as.numeric(t1$edge.length)
        tree1 <- t1
        tree1$node.label <- rnN1$oriN[match(tree1$node.label,
            rnN1$node.label)]
        toDrop <- setdiff(1:length(t1$tip.label), which(!is.na(match(t1$tip.label,
            sp.list$species))))
        t1 <- drop.tip(t1, tip = toDrop)
        Re <- which(!is.na(match(t1$node.label, rnN1$node.label)))
        noRe <- which(is.na(match(t1$node.label, rnN1$node.label)))
        t1$node.label[Re] <- rnN1$oriN[match(t1$node.label, rnN1$node.label)[Re]]
        t1$node.label[noRe] <- ""
    }
    else {
        t1 <- NULL
        tree1 <- NULL
    }

##############

    if ("S2" %in% scenarios) {
        t2r <- replicate(r, list())
        names(t2r) <- paste("run", 1:r, sep = ".")
        tree2r <- replicate(r, list())
        names(tree2r) <- paste("run", 1:r, sep = ".")
        for (o in 1:r) {
            t2 <- tree
            rnN2 <- rnN
            nG <- nodes[nodes$level == "G", ]
            nF <- nodes[nodes$level == "F", ]
            data <- add.tip[add.tip$sort == "F1", ]
            if (dim(data)[1] > 0) {
                for (i in 1:dim(data)[1]) {
                  n <- match(data$family[i], nF$family)
                  g <- nF$gen.n[n]
                  s <- nF$sp.n[n]
                  if (g == 1 & s == 1) {
                    num <- match(nF$taxa[n], t2$tip.label)
                    nlabel <- paste("N", t2$Nnode + 1, sep = "")
                    t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                      node.label = nlabel, position = 2/3)
                    nF$gen.n[n] <- g + 1
                    nF$sp.n[n] <- s + 1
                    x <- which(t2$node.label == nlabel)
                    xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel,
                      oriN = ""))
                    xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
                    rnN2 <- xx
                    nF$bn[n] <- nlabel
                  }
                  else {
                    num <- sample(nG$bn[which(nG$family %in%
                      data$family[i])], 1)
                    t2 <- at.node(t2, location.node = num, tip.label = data$species[i])
                    nF$gen.n[n] <- g + 1
                    nF$sp.n[n] <- s + 1
                  }
                }
            }
            data <- add.tip[add.tip$sort == "F2", ]
            if (dim(data)[1] > 0) {
                for (i in 1:dim(data)[1]) {
                  n <- grep(paste(data$genus[i], "_", sep = ""),
                    t2$tip.label)
                  nlabel <- paste("N", t2$Nnode + 1, sep = "")
                  if (length(n) == 1) {
                    num <- n
                    t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                      node.label = nlabel, position = 1/2)
                    x <- which(t2$node.label == nlabel)
                    xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel,
                      oriN = ""))
                    xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
                    rnN2 <- xx
                  }
                  if (length(n) > 1) {
                    num <- t2$edge[which(t2$edge[, 2] == n[1]),
                      1]
                    t2 <- at.node(t2, location.node = num, tip.label = data$species[i])
                  }
                }
            }
            data <- add.tip[add.tip$sort == "G1", ]
            if (dim(data)[1] > 0) {
                for (i in 1:dim(data)[1]) {
                  n0 <- match(data$genus[i], nG$genus)
                  n <- nG$sp.n[n0]
                  nlabel <- paste("N", t2$Nnode + 1, sep = "")
                  if (n == 1) {
                    num <- t2$tip.label[match(nG$taxa[n0], t2$tip.label)]
                    t2 <- ext.node(t2, location.tip = num, tip.label = data$species[i],
                      node.label = nlabel, position = 1/2)
                    x <- which(t2$node.label == nlabel)
                    xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel,
                      oriN = "", stringsAsFactors = FALSE))
                    xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
                    rnN2 <- xx
                    nG$sp.n[n0] <- nG$sp.n[n0] + 1
                  }
                  if (n > 1) {
                    num <- which(t2$node.label == nG$bn[n0]) +
                      length(t2$tip.label)
                    num1 <- which(t2$edge[, 1] %in% num)
                    part1 <- t2$edge[1:min(num1), ]
                    n1 <- max(which(part1[, 1] < num), 0) + 1
                    part2 <- t2$edge[max(num1):dim(t2$edge)[1],
                      ]
                    n2 <- min(which(part2[, 1] < num), dim(part2)[1] +
                      1) + max(num1) - 2
                    sect <- t2$edge[n1:n2, ]
                    sect <- sort(unique(c(sect[, 1], sect[, 2])))
                    sect <- sect[which(sect > length(t2$tip.label))]
                    num2 <- sect[sample(1:length(sect), 1)]
                    t2 <- at.node(t2, location.node = num2, tip.label = data$species[i])
                  }
                }
            }
            t2$edge.length <- as.numeric(t2$edge.length)
            tree2 <- t2
            tree2$node.label <- rnN2$oriN[match(tree2$node.label,
                rnN2$node.label)]
            toDrop <- setdiff(1:length(t2$tip.label), which(!is.na(match(t2$tip.label,
                sp.list$species))))
            t2 <- drop.tip(t2, tip = toDrop)
            Re <- which(!is.na(match(t2$node.label, rnN2$node.label)))
            noRe <- which(is.na(match(t2$node.label, rnN2$node.label)))
            t2$node.label[Re] <- rnN2$oriN[match(t2$node.label,
                rnN2$node.label)[Re]]
            t2$node.label[noRe] <- ""
            t2r[[o]] <- t2
            tree2r[[o]] <- tree2
        }
    }
    else {
        t2r <- NULL
        tree2r <- NULL
    }
    if ("S3" %in% scenarios) {
        t3 <- tree
        rnN3 <- rnN
        nG <- nodes[nodes$level == "G", ] # genera
        nF <- nodes[nodes$level == "F", ] # families
		
		# F1
		print("Checking F1")
		
        data <- add.tip[add.tip$sort == "F1", ]
		
		print(paste("F1", dim(data)[1], "taxa (4/6)"))
		
        if (dim(data)[1] > 0) {
            pb4 <- progress_bar$new(total=dim(data)[1])
            for (i in 1:dim(data)[1]) {
                n <- match(data$family[i], nF$family)
                g <- nF$gen.n[n]
                s <- nF$sp.n[n]
				
				#print(paste(data$species[i], "(", data$family[i], "), gen=", g, "; sp=", s))
				
                if (g == 1 & s == 1) {
                  num <- match(nF$taxa[n], t3$tip.label)
                  nlabel <- paste("N", t3$Nnode + 1, sep = "")
                  t3 <- ext.node(t3, location.tip = num, tip.label = data$species[i],
                    node.label = nlabel, position = 2/3)
                  nF$gen.n[n] <- g + 1
                  nF$sp.n[n] <- s + 1
                  x <- which(t3$node.label == nlabel)
                  xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel,
                    oriN = ""))
                  xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
                  rnN3 <- xx
                  nF$bn[n] <- nlabel
                }
                if (g == 1 & s > 1) {
                  nlabel <- paste("N", t3$Nnode + 1, sep = "")
                  if ((2/3) * nF$rn.bl[n] <= nF$bn.bl[n]) {
                    len <- (nF$rn.bl[n] - nF$bn.bl[n])/2
                  }
                  if ((2/3) * nF$rn.bl[n] > nF$bn.bl[n]) {
                    len <- nF$rn.bl[n] * 2/3 - nF$bn.bl[n]
                  }
                  port <- len/(nF$rn.bl[n] - nF$bn.bl[n])
                  t3 <- int.node(t3, location.node = nF$bn[n],
                    tip.label = data$species[i], node.label = nlabel,
                    position = port)
                  nF$gen.n[n] <- g + 1
                  nF$sp.n[n] <- s + 1
                  x <- which(t3$node.label == nlabel)
                  xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel,
                    oriN = "", stringsAsFactors = FALSE))
                  xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
                  rnN3 <- xx
                  nF$bn[n] <- nlabel
                }
                if (g > 1) {
                  t3 <- at.node(t3, 
					location.node = nF$bn[n],		#
					#location.node = as.numeric(strsplit(nF$bn[n], split="N")[[1]][2]),
                    tip.label = data$species[i])
                }
                pb4$tick()
            }
        }
        data <- add.tip[add.tip$sort == "F2", ]
		print(paste("F2", dim(data)[1], "taxa (5/6)"))

        if (dim(data)[1] > 0) {
            pb5 <- progress_bar$new(total=dim(data)[1])
            for (i in 1:dim(data)[1]) {
			
                n <- grep(paste(data$genus[i], "_", sep = ""),
                  t3$tip.label)
                nlabel <- paste("N", t3$Nnode + 1, sep = "")
				
				#print(paste(n, "; nlabel=", nlabel, "; length=", length(n)))
				#print(paste("searching ", data$genus[i], "_ in t3$tip.label", sep = ""))
				
                if (length(n) == 1) {
                  t3 <- ext.node(t3, location.tip = n, tip.label = data$species[i],
                    node.label = nlabel, position = 1/2)
                  nG$sp.n[match(data$genus[i], nG$genus)] <- length(n) +
                    1
                  x <- which(t3$node.label == nlabel)
                  xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel,
                    oriN = "", stringsAsFactors = FALSE))
                  xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
                  rnN3 <- xx
                  nG$bn[match(data$genus[i], nG$genus)] <- nlabel
                }
                if (length(n) > 1) {
                  num <- ancestor(t3, min(n), max(n))
                  t3 <- at.node(t3, location.node = num, tip.label = data$species[i])
                }
                pb5$tick()
            }
        }
        data <- add.tip[add.tip$sort == "G1", ]
		
		print(paste("G1", dim(data)[1], "taxa (6/6)"))

        if (dim(data)[1] > 0) {
            pb6 <- progress_bar$new(total=dim(data)[1])
            for (i in 1:dim(data)[1]) {
                n0 <- match(data$genus[i], nG$genus)
                s <- nG$sp.n[n0]
                nlabel <- paste("N", t3$Nnode + 1, sep = "")
				
				#print(paste(data$species[i], "(", data$family[i], "), gen=", g, "; sp=", s))
				
                if (s == 1) {
                  num <- t3$tip.label[match(nG$taxa[n0], t3$tip.label)]
                  t3 <- ext.node(t3, location.tip = num, tip.label = data$species[i],
                    node.label = nlabel, position = 1/2)
                  nG$sp.n[n0] <- nG$sp.n[n0] + 1
                  x <- which(t3$node.label == nlabel)
                  xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel,
                    oriN = "", stringsAsFactors = FALSE))
                  xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
                  rnN3 <- xx
                  nG$bn[n0] <- nlabel
                }
                if (s > 1) {
                  t3 <- at.node(t3, location.node = nG$bn[n0],
                    tip.label = data$species[i])
                }
                pb6$tick()
            }
        }
        t3$edge.length <- as.numeric(t3$edge.length)
        tree3 <- t3
        tree3$node.label <- rnN3$oriN[match(tree3$node.label,
            rnN3$node.label)]
        toDrop <- setdiff(1:length(t3$tip.label), which(!is.na(match(t3$tip.label,
            sp.list$species))))
        t3 <- drop.tip(t3, tip = toDrop)
        Re <- which(!is.na(match(t3$node.label, rnN3$node.label)))
        noRe <- which(is.na(match(t3$node.label, rnN3$node.label)))
        t3$node.label[Re] <- rnN3$oriN[match(t3$node.label, rnN3$node.label)[Re]]
        t3$node.label[noRe] <- ""
    }
    else {
        t3 <- NULL
        tree3 <- t3
    }
    if (output.sp.list == FALSE) {
        sp.list.original <- NULL
    }
    if (output.tree == FALSE) {
        tree1 <- tree2r <- tree3 <- NULL
    }
    if (r == 1) {
        t2r <- t2r$run.1
        tree2r <- tree2r$run.1
    }
    phylo <- list(scenario.1 = t1, scenario.2 = t2r, scenario.3 = t3,
        species.list = sp.list.original, tree.scenario.1 = tree1,
        tree.scenario.2 = tree2r, tree.scenario.3 = tree3)
    phylo[sapply(phylo, is.null)] <- NULL
    return(phylo)
}
ext.node <-
function (phylogeny, location.tip, tip.label, node.label = NULL,
            position = 0.5)
  {
    phylo <- reorder(phylogeny)
    if (!is.numeric(location.tip)) {
      location.tip <- which(phylo$tip.label == location.tip)
    }
    a <- location.tip
    a1 <- which(phylo$edge[, 2] == a)
    h <- phylo$edge[a1, 1]
    if (is.null(phylo$node.label)) {
      if (!is.null(node.label)) {
        nL <- rep(NA, phylo$Nnode + 1)
        n <- h - length(phylo$tip.label)
        nL[n + 1] <- node.label
      }
      if (is.null(node.label)) {
        nL <- NULL
      }
    }
    if (!is.null(phylo$node.label)) {
      if (!is.null(node.label)) {
        n <- h - length(phylo$tip.label)
        nL <- c(phylo$node.label[1:n], node.label)
        if (n < phylo$Nnode) {
          nL <- c(nL, phylo$node.label[(n + 1):phylo$Nnode])
        }
      }
      if (is.null(node.label)) {
        n <- h - length(phylo$tip.label)
        nL <- c(phylo$node.label[1:n], NA)
        if (n < phylo$Nnode) {
          nL <- c(nL, phylo$node.label[(n + 1):phylo$Nnode])
        }
      }
    }
    eG0 <- matrix(c(h + 1, h + 2, h + 2, a, h + 2, a + 1), nrow = 3,
                  byrow = T)
    eG <- matrix(phylo$edge[1:(a1 - 1), ], ncol = 2)
    eG[, 1] <- eG[, 1] + 1
    s <- which(eG[, 1] > (h + 1))
    eG[, 1][s] <- eG[, 1][s] + 1
    s <- which(eG[, 2] > a)
    eG[, 2][s] <- eG[, 2][s] + 1
    s <- which(eG[, 2] > (h + 1))
    eG[, 2][s] <- eG[, 2][s] + 1
    eG <- rbind(eG, eG0)
    tL <- c(phylo$tip.label[1:a], tip.label)
    eL <- c(phylo$edge.length[1:(a1 - 1)], phylo$edge.length[a1] *
              (1 - position), phylo$edge.length[a1] * position, phylo$edge.length[a1] *
              position)
    if (a < length(phylo$tip.label)) {
      eGn <- matrix(phylo$edge[(a1 + 1):(dim(phylo$edge)[1]),
      ], ncol = 2)
      eGn[, 1] <- eGn[, 1] + 1
      s <- which(eGn[, 1] > (h + 1))
      eGn[, 1][s] <- eGn[, 1][s] + 1
      s <- which(eGn[, 2] > a)
      eGn[, 2][s] <- eGn[, 2][s] + 1
      s <- which(eGn[, 2] > (h + 1))
      eGn[, 2][s] <- eGn[, 2][s] + 1
      eG <- rbind(eG, eGn)
      tL <- c(tL, phylo$tip.label[(a + 1):length(phylo$tip.label)])
      eL <- c(eL, phylo$edge.length[(a1 + 1):length(phylo$edge.length)])
    }
    phylo$edge <- eG
    phylo$tip.label <- tL
    phylo$edge.length <- eL
    phylo$Nnode <- phylo$Nnode + 1
    phylo$node.label <- nL
    return(phylo)
  }
 

at.node <-
function (phylogeny, location.node, tip.label)
  {
    phylo <- reorder(phylogeny)
	
	#print("at.node")
	#print(phylo)
	
    if (!is.numeric(location.node))
      location.node <- which(phylo$node.label == location.node) +
        length(phylo$tip.label)
    a <- location.node - length(phylo$tip.label); #print(paste("a =", a))
    EL <- branching.times(phylo)[a]; #print(paste("EL =", EL))
    a0 <- a + length(phylo$tip.label); #print(paste("a0 =", a0))
    a1 <- which(phylo$edge[, 1] == a0)[1]; #print(paste("a1 =", a1)) # producing NaNs
	
	if(is.na(a1)) {
		#print(paste("Did not attach", tip.label))
		return(phylo)
	}
	
    aa <- length(which(phylo$edge[1:(a1 - 1), 2] <= length(phylo$tip.label))); #print(paste("aa =", aa))
    eG0 <- matrix(c(a0 + 1, aa + 1), nrow = 1); #print(paste("eG0 =", eG0))
    eG <- matrix(phylo$edge[a1:dim(phylo$edge)[1], ], ncol = 2); #print(paste("eG =", eG))
    eG[, 1] <- eG[, 1] + 1; #print(paste("eG =", eG))
    s <- which(eG[, 2] > aa); #print(paste("s =", s))
    eG[, 2][s] <- eG[, 2][s] + 1; #print(paste("eG =", eG))
    eG <- rbind(eG0, eG); #print(paste("eG =", eG))
    eL <- c(EL, phylo$edge.length[a1:length(phylo$edge.length)]); #print(paste("eL =", eL))
    tL <- c(tip.label, phylo$tip.label[(aa + 1):length(phylo$tip.label)]); #print(paste("tL =", tL))
    if (a1 > 1) {
      eGn <- matrix(phylo$edge[1:(a1 - 1), ], ncol = 2)
      eGn[, 1] <- eGn[, 1] + 1
      s <- which(eGn[, 2] > aa)
      eGn[, 2][s] <- eGn[, 2][s] + 1
      eG <- rbind(eGn, eG)
      eL <- c(phylo$edge.length[1:(a1 - 1)], eL)
    }
    if (aa > 0) {
      tL <- c(phylo$tip.label[1:aa], tL)
    }
    phylo$edge <- eG
    phylo$tip.label <- tL
    phylo$edge.length <- eL
    return(phylo)
}

build.nodes.1.verbose<-
function (tree, tips)
{
  tips <- data.frame(as.matrix(tips), stringsAsFactors = FALSE)
  dimnames(tips)[[2]] <- tolower(dimnames(tips)[[2]])
  tips$species <- gsub("(^[[:alpha:]])", "\\U\\1", tips$species,
                       perl = TRUE)
  tips$genus <- gsub("(^[[:alpha:]])", "\\U\\1", tips$genus,
                     perl = TRUE)
  tips$family <- gsub("(^[[:alpha:]])", "\\U\\1", tips$family,
                      perl = TRUE)
  tree$node.label <- paste("N", 1:tree$Nnode, sep = "")
  node.dep <- branching.times(tree)
  tips <- tips[match(tree$tip.label, tips$species), ]
  tips$No. <- as.integer(dimnames(tips)[[1]])
  clustF <- cluster(tips$family)
  tips$numF <- clustF
  clustG <- cluster(tips$genus)
  tips$numG <- clustG
  clustF.size <- table(clustF)
  clustF.size <- data.frame(numF = as.integer(names(clustF.size)),
                            sizeF = as.integer(clustF.size))
  clustG.size <- table(clustG)
  clustG.size <- data.frame(numG = as.integer(names(clustG.size)),
                            sizeG = as.integer(clustG.size))
  tips1 <- merge(tips, clustF.size, all.x = T)
  tips1 <- merge(tips1, clustG.size, all.x = T)
  tips1 <- tips1[, c("No.", "species", "genus", "family", "numF",
                     "numG", "sizeF", "sizeG")]
  xF <- tips1[!duplicated(tips1$numF), ]
  xF <- xF[rev(order(xF$sizeF)), ]
  xF <- xF[!duplicated(xF$family), ]
  xx <- merge(tips1, data.frame(numF = unique(xF$numF)))
  xxG <- table(xx[!duplicated(xx$genus), ]$numF)
  xxG <- data.frame(numF = as.integer(names(xxG)), gen.n = as.integer(xxG))
  xxS <- table(xx$numF)
  xxS <- data.frame(numF = as.integer(names(xxS)), sp.n = as.integer(xxS))
  xxGS <- merge(xxG, xxS)
  xF1 <- xF[xF$sizeF == 1, ]
  if (dim(xF1)[1] > 0) {
    Fn1 <- data.frame(level = "F", family = xF1$family, genus = "",
                      rn = "", rn.bl = 0, bn = "", bn.bl = 0, taxa = xF1$species,
                      stringsAsFactors = FALSE)
    print("Checking Fn1 (1/4)")
    pba <- progress_bar$new(total= dim(Fn1)[1])
    for (i in 1:dim(Fn1)[1]) {
      x0 <- match(Fn1$taxa[i], tree$tip.label)
      x1 <- match(x0, tree$edge[, 2])
      x2 <- tree$edge[x1, 1] - length(tree$tip.label)
      Fn1$bn[i] <- tree$node.label[x2]
      Fn1$rn[i] <- tree$node.label[x2]
      Fn1$rn.bl[i] <- tree$edge.length[x1]
      Fn1$bn.bl[i] <- tree$edge.length[x1]
      pba$tick()
    }
  }
  if (dim(xF1)[1] == 0) {
    Fn1 <- NULL
  }
  xF2 <- xF[xF$sizeF > 1, ]
  if (dim(xF2)[1] > 0) {
    tipsF <- merge(tips1, data.frame(numF = xF2$numF))
    tF <- drop.tip(tree, setdiff(tree$tip.label, tipsF$species))
    tipsF <- tipsF[match(tF$tip.label, tipsF$species), ]
    tF$tip.label <- tipsF$family
    n1 <- which(!duplicated(tF$tip.label))
    n2 <- length(tF$tip.label) + 1 - which(!duplicated(rev(tF$tip.label)))
    nn <- sort(unique(c(n1, n2)))
    tF <- drop.tip(tF, setdiff(1:length(tF$tip.label), nn))
    Fn2 <- data.frame(level = "F", family = xF2$family, genus = "",
                      rn = "", rn.bl = 0, bn = "", bn.bl = 0, taxa = "",
                      stringsAsFactors = FALSE)
    x <- 1:length(tF$tip.label)
    print("Checking Fn2 (2/4)")
    pbb <- progress_bar$new(total= dim(Fn2)[1])
    for (i in 1:dim(Fn2)[1]) {
      h <- which(tF$tip.label == Fn2$family[i])
      tt <- drop.tip(tF, setdiff(x, h))
      Fn2$bn[i] <- tt$node.label
      Fn2$bn.bl[i] <- tt$edge.length[1]
      n <- which(tree$node.label == Fn2$bn[i]) + length(tree$tip.label)
      if (n > min(tree$edge[, 1])) {
        n1 <- tree$edge[which(tree$edge[, 2] == n), 1] -
          length(tree$tip.label)
        Fn2$rn[i] <- tree$node.label[n1]
        Fn2$rn.bl[i] <- node.dep[Fn2$rn[i]]
      }
      if (n == min(tree$edge[, 1])) {
        Fn2$rn[i] <- Fn2$bn[i]
        Fn2$rn.bl[i] <- Fn2$bn.bl[i]
      }
      pbb$tick()
    }
  }
  if (dim(xF2)[1] == 0) {
    xF2 <- NULL
  }
  Fn <- rbind(Fn1, Fn2)
  Fn <- merge(Fn, xF[, c("family", "numF")])
  Fn <- merge(Fn, xxGS)
  Fn <- Fn[, c("level", "family", "genus", "rn", "rn.bl", "bn",
               "bn.bl", "gen.n", "sp.n", "taxa")]
  xG <- tips1[!duplicated(tips1$numG), ]
  xG <- xG[rev(order(xG$sizeG)), ]
  xG <- xG[!duplicated(xG$genus), ]
  xx <- merge(tips1, data.frame(numG = unique(xG$numG)))
  xxG <- table(xx$numG)
  xxG <- data.frame(numG = as.integer(names(xxG)), gen.n = 1)
  xxS <- table(xx$numG)
  xxS <- data.frame(numG = as.integer(names(xxS)), sp.n = as.integer(xxS))
  xxGS <- merge(xxG, xxS)
  xG1 <- xG[xG$sizeG == 1, ]
  if (dim(xG1)[1] > 0) {
    Gn1 <- data.frame(level = "G", family = xG1$family, genus = xG1$genus,
                      rn = "", rn.bl = 0, bn = "", bn.bl = 0, taxa = xG1$species,
                      stringsAsFactors = FALSE)
    print("Checking Gn1 (3/4)")
    pbc <- progress_bar$new(total= dim(Gn1)[1])
    for (i in 1:dim(Gn1)[1]) {
      x0 <- match(Gn1$taxa[i], tree$tip.label)
      x1 <- match(x0, tree$edge[, 2])
      x2 <- tree$edge[x1, 1] - length(tree$tip.label)
      Gn1$bn[i] <- tree$node.label[x2]
      Gn1$rn[i] <- tree$node.label[x2]
      Gn1$rn.bl[i] <- tree$edge.length[x1]
      Gn1$bn.bl[i] <- tree$edge.length[x1]
      pbc$tick()
    }
  }
  if (dim(xG1)[1] == 0) {
    Gn1 <- NULL
  }
  xG2 <- xG[xG$sizeG > 1, ]
  if (dim(xG2)[1] > 0) {
    tipsG <- merge(tips1, data.frame(numG = xG2$numG))
    tG <- drop.tip(tree, setdiff(tree$tip.label, tipsG$species))
    tipsG <- tipsG[match(tG$tip.label, tipsG$species), ]
    tG$tip.label <- tipsG$genus
    n1 <- which(!duplicated(tG$tip.label))
    n2 <- length(tG$tip.label) + 1 - which(!duplicated(rev(tG$tip.label)))
    nn <- sort(unique(c(n1, n2)))
    tG <- drop.tip(tG, setdiff(1:length(tG$tip.label), nn))
    Gn2 <- data.frame(level = "G", family = xG2$family, genus = xG2$genus,
                      rn = "", rn.bl = 0, bn = "", bn.bl = 0, taxa = "",
                      stringsAsFactors = FALSE)
    x <- 1:length(tG$tip.label)
    print("Checking Gn2 (4/4)")
    pbd <- progress_bar$new(total= dim(Gn2)[1])
    for (i in 1:dim(Gn2)[1]) {
      h <- which(tG$tip.label == Gn2$genus[i])
      tt <- drop.tip(tG, setdiff(x, h))
      Gn2$bn[i] <- tt$node.label
      Gn2$bn.bl[i] <- tt$edge.length[1]
      n <- which(tree$node.label == Gn2$bn[i]) + length(tree$tip.label)
      if (n > min(tree$edge[, 1])) {
        n1 <- tree$edge[which(tree$edge[, 2] == n), 1] -
          length(tree$tip.label)
        Gn2$rn[i] <- tree$node.label[n1]
        Gn2$rn.bl[i] <- node.dep[Gn2$rn[i]]
      }
      if (n == min(tree$edge[, 1])) {
        Gn2$rn[i] <- Gn2$bn[i]
        Gn2$rn.bl[i] <- Gn2$bn.bl[i]
      }
      pbd$tick()
    }
  }
  if (dim(xG2)[1] == 0) {
    Gn2 <- NULL
  }
  Gn <- rbind(Gn1, Gn2)
  Gn <- merge(Gn, xG[, c("genus", "numG")])
  Gn <- merge(Gn, xxGS)
  Gn <- Gn[, c("level", "family", "genus", "rn", "rn.bl", "bn",
               "bn.bl", "gen.n", "sp.n", "taxa")]
  nodes <- rbind(Fn, Gn)
  return(nodes)
}
