#' View motif
#'
#' @description View information content for each position of the PWM.
#' Information content is modelled using Shannon's Entropy Model. The maximum
#' information content is therefore log2(n), where n is the number of amino
#' acids. Colors of Amino Acids are in accordance with the Lesk scheme.
#' @param pwm_in View a PWM provided using the buildPWM. Default = NULL
#' @param which_pwm If pwms are input (outputs of buildPWM), a kinase name must
#' match a name in pwms$kinase$kinase list of names. Default = NULL
#' @param fontsize Font size to use on x and y axis. Default = 10
#' @param view_pwm View the PWM. Default = FALSE
#' @param color_scheme Which color scheme to use for Amino Acid Groups. Options
#' are "lesk" or "shapely". Default = "shapely"
#' @param convert_PWM pwm_in is a matrix of counts at position. TRUE will 
#' convert this matrix to a PWM. Default = FALSE
#' @param pseudo Small amount added to the PWM model, where zero's exist, to 
#' avoid log zero. Default = 0.01
#' @param correction_factor Number of sequences used to infer the PWM. 
#' This can be used where a small number of sequences were used to build the 
#' model and included as E_n in the Shannon's Entropy Model. Default = NULL
#'
#' @examples
#' ## Build PWM models from phosphositeplus data with default of minimum
#' ## of 10 substrate sequences for building a PWM model.
#' data(phosphositeplus_human)
#' ##randomly sample 1000 substrates for demonstration.
#' set.seed(1)
#' sample_pwm <- phosphositeplus_human[sample(nrow(phosphositeplus_human), 
#' 1000),]
#' pwms <- buildPWM(sample_pwm)
#'
#' ## Data frame of models built and number of sequences used to build each
#' ## PWM model:
#' head(pwms$kinase)
#' ## Will not visualise the motif
#' CAMK2A_motif <- viewPWM(pwm_in = pwms, 
#'                         which_pwm = "CAMK2A",
#'                         view_pwm = FALSE)
#' # Use view_pwm = TRUE to view the motif
#' @return Visualisation of a motif, scaled on bits and two tables. 1) pwm: 
#' corresponding to the PWM from pwm and 2) pwm_bits: corresponding to the
#' conversion to bits.
#'
#' @export viewPWM
#' 
#' @importFrom grid grid.newpage grid.polygon gpar pushViewport plotViewport
#' @importFrom grid dataViewport grid.xaxis gpar grid.yaxis
#' @importFrom grDevices rgb

#this function is for delivering a PWM to a frame:
viewPWM <- function(pwm_in = NULL, 
                    which_pwm = NULL,
                    fontsize = 10,
                    view_pwm = FALSE,
                    pseudo = 0.01,
                    convert_PWM = FALSE,
                    color_scheme = "shapely",
                    correction_factor = NULL){
  
  #----------------------------------------------
  if (is.null(pwm_in) && convert_PWM == FALSE)
    stop(
      "pwm_in not provided; you must provide an input table containing
      computed position weight matrices using buildPWM()"
    )
  if (!is.null(pwm_in) && is.null(which_pwm))
    stop(
      "which_pwm not provided; you must provide a PWM to search for"
    )
  if (!is.null(pwm_in) && !is.list(pwm_in))
    stop(
      "pwm_in is not a list format; something has gone wrong. Make sure
      you compute the position weight matrices using buildPWM()"
    )
  if (!(color_scheme %in% c("shapely", "lesk")))
    stop("color_scheme must be either shapely or lesk")
  
  #----------------------------------------------
  
  # PWM already converted
  if(convert_PWM == FALSE){
    # if a list is input and to be extracted:
    rownames(pwm_in$kinase) <- seq_len(nrow(pwm_in$kinase))
    which_pwm <- pwm_in$kinase[pwm_in$kinase$kinase==which_pwm,]
    which_pwm <- pwm_in$pwm[[as.numeric(rownames(which_pwm))]]
  }
  # data is input as count table
  if(convert_PWM == TRUE){
      which_pwm <- pwm_in
      #which_pwm <- SPRK2_motif_1
      wild_card <- "_"
      uniq_AA <- c(wild_card, "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
      missing_AA <- setdiff(uniq_AA, rownames(which_pwm))
      missing_AA <- data.frame(matrix(data = 0, ncol = ncol(which_pwm), 
                           nrow = length(missing_AA)), row.names = missing_AA)
      colnames(missing_AA) <- colnames(which_pwm)
      which_pwm <- rbind(missing_AA, which_pwm)
      # remove wildcard from motif scoring
      which_pwm <- which_pwm[!rownames(which_pwm) %in% c(wild_card),]
      which_pwm <- which_pwm[order(rownames(which_pwm)),]
      # col sums
      col_counts <- apply(which_pwm, 2, sum, na.rm = TRUE)
      # conver to PPM
      ppm <- sapply(1:ncol(which_pwm), function(i) 
        which_pwm[,i] / col_counts[i])
      rownames(ppm) <- rownames(which_pwm)
      # Generate Position Weight Matrix
      # addition of pseudo count
      which_pwm <- log2((ppm / (1 / nrow(ppm)))+pseudo)
      
      colnames(which_pwm) <- paste("p", seq_len(ncol(which_pwm)), sep = "")
    }
    
  # select window width = size of motif
  window_size <- ncol(which_pwm)
  half_window <- (window_size-1)/2
  # inverse log
  pwm <- 2^(which_pwm)/20
  # max bits, uncertainty with 20 AAs
  max_bits <- log2((1/0.05)+pseudo) 
  # Calculate Shannon Entropy, uncertainty at position (l)
  H_l <- pwm * log2(pwm)
  H_l <- -apply(H_l, 2, sum)
  if(is.null(correction_factor)){
    R_l <- max_bits - (H_l)# + e_n)
  }
  if(!is.null(correction_factor)){
    e_n <- correction_factor
    R_l <- max_bits - ((H_l) + e_n)
  }
  # multiplication step for vector, each col, to find Height (information)
  Height <- sapply(1:ncol(pwm), function(i) 
                    pwm[,i] * R_l[i])
  colnames(Height) <- colnames(pwm)
  # transform pwm to the Height
  pwm <- Height
  y_max <- 1
  # convert to 0-1 scale
  pwm <- pwm/max_bits
  
  if(view_pwm == TRUE){
    # start new grid for viewing and set up margins
    grid.newpage()
    pushViewport(plotViewport(margins = c(4,3,2,2)))
    pushViewport(dataViewport(xscale = c(0,1), yscale = c(0,1), name="vp1"))
    # label axis
    grid.xaxis(at = (seq(0,(window_size-1))/window_size)+((1/(window_size-1))/2),
               label = c(-half_window:half_window), 
               gp = gpar(fontsize = fontsize))
    grid.yaxis(at = c(0, 0.5, 1),
               label = c(0, round(max_bits/2, 2), round(max_bits, 2)), 
               gp = gpar(fontsize = fontsize))
    # add letters to the grid
    sapply(seq_len(ncol(pwm)), function(i)
      sapply(seq_len(nrow(pwm)), function(j)
        add.letter(letter = rownames(pwm[order(pwm[,i]),])[j], 
                   x.position = i, 
                   y.position = sum(pwm[order(pwm[,i]),][1:j-1,i]),
                   prob = pwm[order(pwm[,i]),][j,i],
                   window.size = window_size,
                   color_scheme = color_scheme
                   )
      )
    )
  }
return(list("pwm" = which_pwm, 
            "pwm_bits" = pwm))
}


# helper function to select the letters
add.letter <- function(letter, 
                       x.position, 
                       y.position, 
                       prob, 
                       window.size,
                       color_scheme){

  if(letter == "A"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 255,
                      green = 150,
                      blue = 0, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 200,
                      green = 200,
                      blue = 200, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)/10, 
                y = c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)/10, 
                id = rep(1, 13), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
                )
    }
  if(letter == "C"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 230,
                      green = 230,
                      blue = 0, 
                      maxColorValue = 255)
    }
    ang <- seq((pi/2), (2*pi)+(pi/2), length = 360)
    x.out <- ((sin(ang))+1)/2; x.out <- x.out[40:(length(x.out)-40)]
    y.out <- ((cos(ang))+1)/2; y.out <- y.out[40:(length(y.out)-40)]
    x.in <- ((sin(ang)*0.7)+1)/2;x.in <- x.in[40:(length(x.in)-40)]
    y.in <- ((cos(ang)*0.7)+1)/2;y.in <- y.in[40:(length(y.in)-40)]
    add.polygon(x = c(x.out, rev(x.in)), 
                y = c(y.out, rev(y.in)), 
                id = rep(1,length(c(x.out, rev(x.in)))), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
                )
    }
  if(letter == "D"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 230,
                      green = 10,
                      blue = 10, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 230,
                      green = 10,
                      blue = 10, 
                      maxColorValue = 255)
    }
    ang <- seq((pi), (2*pi)+(pi), length = 360)
    x.out <- (((sin(ang)))/2)+0.5;x.out <- x.out[180:(length(x.out))]
    y.out <- ((cos(ang))+1)/2;y.out <- y.out[180:(length(y.out))]
    x.out <- c(x.out, 0, 0);y.out <- c(y.out, 0, 1)
    x.in <- (((((sin(ang)))/2)+0.5)*0.7)+0.15;x.in <- x.in[180:(length(x.in))]
    y.in <- ((((cos(ang))*0.7)+1)/2);y.in <- y.in[180:(length(y.in))]
    x.in <- c(x.in, 1-max(x.in), 1-max(x.in));y.in <- c(y.in, min(y.in), max(y.in))
    add.polygon(x = c(x.out, x.in), 
                y = c(y.out, y.in), 
                id = c(rep(1,length(x.out)), 
                       rep(2,length(x.in))), 
                fill = c(fill_rgb, 
                         rgb(red = 255,
                             green = 255,
                             blue = 255, 
                             maxColorValue = 255)),
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "E"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 230,
                      green = 10,
                      blue = 10, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 230,
                      green = 10,
                      blue = 10, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,0,10,10,2,2,6,6,2,2,10,10)/10, 
                y = c(0,10,10,8,8,6,6,4,4,2,2,0)/10, 
                id = rep(1, 12), 
                fill = fill_rgb, 
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "F"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 50,
                      green = 50,
                      blue = 170, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,0,10,10,2,2,10,10,2,2)/10, 
                y = c(0,10,10,8,8,6,6,4,4,0)/10, 
                id = rep(1, 10), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "G"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 255,
                      green = 150,
                      blue = 0, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 235,
                      green = 235,
                      blue = 235, 
                      maxColorValue = 255)
    }
    ang <- seq((pi/2), (2*pi)+(pi/2), length = 360)
    x.out <- ((sin(ang))+1)/2;x.out <- x.out[40:(length(x.out)-40)]
    y.out <- ((cos(ang))+1)/2;y.out <- y.out[40:(length(y.out)-40)]
    x.in <- ((sin(ang)*0.7)+1)/2;x.in <- x.in[40:(length(x.in)-40)]
    y.in <- ((cos(ang)*0.7)+1)/2;y.in <- y.in[40:(length(y.in)-40)]
    x <- c(x.in[1], x.in[1], 0.5,0.5, max(x.out),  x.out, rev(x.in))
    y <- c( y.in[1],0.4,0.4,0.6,0.6, y.out, rev(y.in))
    add.polygon(x = x, 
                y = y, 
                id = rep(1,length(x)), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "H"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 255,
                      green = 255,
                      blue = 255, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 130,
                      green = 130,
                      blue = 210, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,0,2,2,8,8,10,10,8,8,2,2,0)/10, 
                y = c(0,10,10,6,6,10,10,0,0,4,4,0,0)/10, 
                id = rep(1, 13), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "I"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,0,4,4,0,0,10,10,6,6,10,10,0)/10, 
                y = c(0,2,2,8,8,10,10,8,8,2,2,0,0)/10, 
                id = rep(1, 13), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "K"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 20,
                      green = 90,
                      blue = 255, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 20,
                      green = 90,
                      blue = 255, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,0,2,2,7,10,4,10,7,2,2,0)/10, 
                y = c(0,10,10,6,10,10,5,0,0,4,0,0)/10, 
                id = rep(1, 12), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "L"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,0,2,2,10,10,0)/10, 
                y = c(0,10,10,2,2,0,0)/10, 
                id = rep(1, 7), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "M"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 230,
                      green = 230,
                      blue = 0, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,0,2,5,8,10,10,8,8,5,2,2,0)/10, 
                y = c(0,10,10,6,10,10,0,0,5,2,5,0,0)/10, 
                id = rep(1, 13), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "N"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 255,
                      green = 0,
                      blue = 255, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 0,
                      green = 220,
                      blue = 220, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,0,3,8,8,10,10,7,2,2,0)/10, 
                y = c(0,10,10,3,10,10,0,0,7,0,0)/10, 
                id = rep(1, 11), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "P"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 250,
                      green = 150,
                      blue = 130, 
                      maxColorValue = 255)
    }
    ang <- seq((pi),(2*pi)+(pi), length = 360)
    x.out <- (((sin(ang)))/2)+0.5
    x.out <- c(x.out[180:(length(x.out))],0.2,0.2,0,0) 
    y.out <- (((cos(ang))+1)/4)+0.5
    y.out <- c(y.out[180:(length(y.out))],0.5,0,0,1)
    x.in <- (((((sin(ang)))/2)+0.5)*0.5)+0.30
    x.in <- c(x.in[180:(length(x.in))],0.2,0.2)
    y.in <- (((((cos(ang))+1)/4)+0.5)*0.5)+0.375
    y.in <- c(y.in[180:(length(y.in))],min(y.in), max(y.in))
    add.polygon(x = c(x.out, x.in), 
                y = c(y.out, y.in), 
                id = c(rep(1,length(x.out)), 
                       rep(2,length(x.in))), 
                fill = c(fill_rgb,
                         rgb(red = 255,
                             green = 255,
                             blue = 255, 
                             maxColorValue = 255)),
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "Q"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 255,
                      green = 0,
                      blue = 255, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 0,
                      green = 220,
                      blue = 220, 
                      maxColorValue = 255)
    }
    ang <- seq((pi/2), (2*pi)+(pi/2), length = 360)
    x.out <- ((sin(ang))+1)/2
    y.out <- ((cos(ang))+1)/2
    x.in <- ((sin(ang)*0.7)+1)/2
    y.in <- ((cos(ang)*0.7)+1)/2
    x <- c(x.out, rev(x.in))*0.9
    y <- c(y.out, rev(y.in))
    x.q <- c(0.9, 0.55, 0.65, 1)
    y.q <- c(0, 0.35, 0.45, 0.1)
    add.polygon(x = c(x, x.q), 
                y = c(y, y.q), 
                id = c(rep(1,length(x)), 
                       rep(2,length(x.q))), 
                fill = c(fill_rgb,
                         fill_rgb),
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "R"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 20,
                      green = 90,
                      blue = 255, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 20,
                      green = 90,
                      blue = 255, 
                      maxColorValue = 255)
    }
    ang <- seq((pi), (2*pi)+(pi), length = 360)
    x.out <- (((sin(ang)))/2)+0.5
    x.out <- c(x.out[180:(length(x.out))], 0.2, 0.2, 0, 0)
    y.out <- (((cos(ang))+1)/4)+0.5
    y.out <- c(y.out[180:(length(y.out))], 0.5, 0, 0, 1)
    x.in <- (((((sin(ang)))/2)+0.5)*0.5)+0.30
    x.in <- c(x.in[180:(length(x.in))],0.2, 0.2)
    y.in <- (((((cos(ang))+1)/4)+0.5)*0.5)+0.375
    y.in <- y.in[180:(length(y.in))];
    y.in <- c(y.in, min(y.in), max(y.in))
    x.r <- c(0.8, 0.4, 0.6, 1);y.r <- c(0, 0.55, 0.55, 0)
    add.polygon(x = c(x.out, x.in, x.r), 
                y = c(y.out, y.in, y.r), 
                id = c(rep(1,length(x.out)), 
                       rep(2,length(x.in)), 
                       rep(3,length(x.r))), 
                fill = c(fill_rgb,
                         rgb(red = 255,
                             green = 255,
                             blue = 255, 
                             maxColorValue = 255),
                         fill_rgb),
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "S"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 250,
                      green = 150,
                      blue = 0, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 250,
                      green = 150,
                      blue = 0, 
                      maxColorValue = 255)
    }
    angle1 <- seq((pi), (2*pi)+(pi), length = 360)
    #outer and inner
    x.out <- (((sin(angle1)))/2)+0.5
    x.out <- (x.out[100:(length(x.out))]*-1)+1
    y.out <- (((cos(angle1))+1)*0.3)
    y.out <- y.out+(1-max(y.out))
    y.out <- y.out[100:(length(y.out))]
    x.out.b <- (x.out*-1)+1
    y.out.b <- (y.out*-1)+1
    x.in <- (x.out*0.5)+0.25
    y.in <- (y.out*0.4)
    y.in <- y.in + (max(y.out.b)-min(y.in))
    x.in.b <- (x.in*-1)+1
    y.in.b <- (y.in*-1)+1
    x <- c(x.out, rev(x.in.b), x.out.b, rev(x.in))
    y <- c(y.out, rev(y.in.b), y.out.b, rev(y.in))
    add.polygon(x = x, 
                y = y, 
                id = rep(1,length(x)), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "T"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 250,
                      green = 150,
                      blue = 0,
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 250,
                      green = 150,
                      blue = 0, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0,10,10,6,6,4,4,0)/10, 
                y = c(10,10,8,8,0,0,8,8)/10, 
                id = rep(1, 8), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "V"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15,
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(0.4,0,0.2,0.5,0.8,1,0.6), 
                y = c(0,1,1,0.2,1,1,0), 
                id = rep(1, 7), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "W"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15, 
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 180,
                      green = 90,
                      blue = 180,
                      maxColorValue = 255)
    }
    add.polygon(x = c(2,0,2,3.5,5,6.5,8,10,8,6,5,4)/10, 
                y = c(0,10,10,3,6,3,10,10,0,0,3,0)/10, 
                id = rep(1, 12), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "Y"){
    if(color_scheme == "lesk"){
      fill_rgb <- rgb(red = 15,
                      green = 130,
                      blue = 15,
                      maxColorValue = 255)
    }
    if(color_scheme == "shapely"){
      fill_rgb <- rgb(red = 50,
                      green = 50,
                      blue = 170, 
                      maxColorValue = 255)
    }
    add.polygon(x = c(4,4,0,2,5,8,10,6,6)/10, 
                y = c(0,5,10,10,6.5,10,10,5,0)/10, 
                id = rep(1, 9), 
                fill = fill_rgb,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "_"){
    add.polygon(x = , 
                y = , 
                id = , 
                fill = ,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
  if(letter == "X"){
    add.polygon(x = c(0,0,10,10)/10, 
                y = c(0,10,10,0)/10, 
                id = rep(1, 4), 
                fill = 1,
                x.factor = 1/window.size, 
                x.pos = x.position, 
                y.pos = y.position, 
                y.scale = prob
    )
  }
}

add.polygon <- function(x, y, id, fill, x.factor, x.pos, y.pos, y.scale){
  x.pos = x.pos-1
  grid.polygon(x=(x*x.factor)+(x.pos*x.factor),
               y=(y*y.scale)+y.pos,
               id = id,
               gp = gpar(fill = fill, col="transparent"))
}