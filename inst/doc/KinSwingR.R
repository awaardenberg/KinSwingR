## ----eval=FALSE----------------------------------------------------------
#  library(KinSwingR)

## ----eval=FALSE----------------------------------------------------------
#  data(example_phosphoproteome)
#  data(phosphositeplus_human)

## ----eval=FALSE----------------------------------------------------------
#  head(example_phosphoproteome)
#  
#  ##                              annotation peptide          fc        pval
#  ## 1      A0A096MJ61|NA|89|PRRVRNLSAVLAART      NA -0.08377538 0.218815889
#  ## 2 A0A096MJB0|Adcy9|1296|LDKASLGSDDGAQTK      NA  0.03707147 0.751069301
#  ## 3  A0A096MJB0|Adcy9|610|PRGQGTASPGSVSDL      NA -0.06885408 0.594494965
#  ## 4  A0A096MJB0|Adcy9|613|QGTASPGSVSDLAQT      NA -0.29418446 0.002806832
#  ## 5   A0A096MJN4|Sept4|49|ILEPRPQSPDLCDDD      NA  0.09097982 0.078667811
#  ## 6   A0A096MJN4|Sept4|81|FCPPAPLSPSSRPRS      NA -0.12246661 0.078619010

## ----eval=FALSE----------------------------------------------------------
#  head(phosphositeplus_human)
#  ##      kinase    substrate
#  ## [1,] "EIF2AK1" "MILLSELSRRRIRSI"
#  ## [2,] "EIF2AK1" "RILLSELSR______"
#  ## [3,] "EIF2AK1" "IEGMILLSELSRRRI"
#  ## [4,] "PRKCD"   "MKKKDEGSYDLGKKP"
#  ## [5,] "PRKCD"   "FPLRKTASEPNLKVR"
#  ## [6,] "PRKCD"   "PLLARSPSTNRKYPP"
#  

## ----eval=FALSE----------------------------------------------------------
#  annotated.data <- clean.annotation(input.data = example_phosphoproteome,
#                                     annotation.delimiter = "|",
#                                     multi.protein.delimiter = ":",
#                                     multi.site.delimiter = ";",
#                                     seq.number = 4,
#                                     replace = TRUE,
#                                     replace.search = "X",
#                                     replace.with = "_")
#  
#  head(annotated.data)
#  ##                              annotation         peptide          fc        pval
#  ## 1      A0A096MJ61|NA|89|PRRVRNLSAVLAART PRRVRNLSAVLAART -0.08377538 0.218815889
#  ## 2 A0A096MJB0|Adcy9|1296|LDKASLGSDDGAQTK LDKASLGSDDGAQTK  0.03707147 0.751069301
#  ## 3  A0A096MJB0|Adcy9|610|PRGQGTASPGSVSDL PRGQGTASPGSVSDL -0.06885408 0.594494965
#  ## 4  A0A096MJB0|Adcy9|613|QGTASPGSVSDLAQT QGTASPGSVSDLAQT -0.29418446 0.002806832
#  ## 5   A0A096MJN4|Sept4|49|ILEPRPQSPDLCDDD ILEPRPQSPDLCDDD  0.09097982 0.078667811
#  ## 6   A0A096MJN4|Sept4|81|FCPPAPLSPSSRPRS FCPPAPLSPSSRPRS -0.12246661 0.078619010

## ----eval=FALSE----------------------------------------------------------
#  pwms <- build.pwm(phosphositeplus_human)

## ----eval=FALSE----------------------------------------------------------
#  head(pwms$kinase)
#  ##    kinase   n
#  ## 1 EIF2AK1   3
#  ## 2   PRKCD 138
#  ## 3    PIM2  10
#  ## 4  CAMK2A 189
#  ## 5 CSNK2A1 590
#  ## 6    VRK2   9

## ----eval=FALSE----------------------------------------------------------
#  scores <- score.sequences(input.data = annotated.data,
#                            pwm.in = pwms,
#                            threads = 4)

## ----eval=FALSE----------------------------------------------------------
#  swing.out <- swing(input.data = annotated.data,
#                    pwm.in = pwms,
#                    pwm.scores = scores,
#                    threads = 4)
#  
#  # This will produce two tables, one is a network for use with e.g. Cytoscape and the other is the scores. To access the scores:
#  
#  head(swing$scores)
#  ##      kinase pos neg all        pk        nk swing.raw   n    swing   p.greater   p.less
#  ## 78  CSNK2A1  12   6  18 0.6666667 0.3333333 12.590843 590 2.661780 0.009990010 0.989011
#  ## 33   CAMK2A  42  24  66 0.6363636 0.3636364 12.078387 189 2.577006 0.006993007 0.988012
#  ## 276   PRKG2  34  15  49 0.6938776 0.3061224  8.466674  15 1.979529 0.000999001 0.999001
#  ## 130   IKBKB  26  15  41 0.6341463 0.3658537  7.991465  52 1.900917 0.015984016 0.984016
#  ## 133    INSR  14   9  23 0.6086957 0.3913043  5.500727  52 1.488880 0.045954046 0.953047
#  ## 111     FGR  21  10  31 0.6774194 0.3225806  5.434642   8 1.477948 0.002997003 0.995005

## ----eval=FALSE----------------------------------------------------------
#  one.call <- swing.master(kinase.table = phosphositeplus_human,
#                           input.data = annotated.data,
#                           threads = 4)

