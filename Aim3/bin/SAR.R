################################################################################
#                                                                              #
# Functions for generating species-area relationships                          #
# Written by: Ken Locey                                                        #
#                                                                              #
# Notes: This code provides functions to construct species-area relationship   #
#                                                                              #
# To-Do List:                                                                  #
#         1.                                                                   #
#         2. Add warnings                                                      #
#                                                                              #
################################################################################

iterations = 10

# A function to generate the species-area relationship by
# accumulating area according to distance

SAR.accum.dist <- function(com, geo.dist){
  Alist <- c()
  Slist <- c()
  num.ponds <- c(1, 2, 4, 6, 8, 12, 16, 24, 32, 42, 48)
  
  for (i in num.ponds) {   
    areas <- c() # hold iterated area values 
    Ss <- c() # hold iterated S values
    
    for(j in 2:iterations){
      pondID <- sample(51, size = 1)
      while(pondID == 29 | pondID == 31){ 
        pondID <- sample(51, size = 1) 
      }
      
      area <- as.numeric(pond.areas[pondID]) # aggregating area
      cum.abs <- com[pondID, ]
      used <- c()
      
      for (k in 2:i) { # Loop through ponds
        sdata <- subset(coord.dist.ls, FALSE == is.element(NBX, used) & FALSE == is.element(NBY, used))
        sdata <- subset(sdata, NBX == pondID | NBY == pondID)
        sdata <- subset(sdata, NBX != 29 | NBY != 31)
        sdata <- subset(sdata, geo.dist == min(sdata[, 3]))
        
        if (dim(sdata)[1] > 1) {
          x <- sample(dim(sdata)[1], size=1)
          sdata <- sdata[x,]
          }
        
        sdata <- t(as.matrix(as.numeric(as.matrix(sdata))))
        used <- c(used, as.integer(pondID))
        area <- area + as.numeric(pond.areas[pondID]) # aggregating area
        cum.abs <- cum.abs + com[pondID, ]
        
        # update pondID to find next nearest neighbor
        if (sdata[1] == pondID) {
          pondID <- sdata[2]
        } else {
          pondID <- sdata[1]
        }
        
      }
      Ss <- c(Ss, length(cum.abs[cum.abs > 0]))
      areas <- c(areas, area)
    }
    # End random pond samples loop  
    Alist <- rbind(Alist, mean(areas))
    Slist <- rbind(Slist, mean(Ss))
    #print(c(mean(areas), mean(Ss)))
  }
  
  return(cbind(log10(Alist), log10(Slist)))
  #return(cbind(Alist, Slist))
}


# A function to generate the species-area relationship by
# Random Accumulating Sites

SAR.rand.accum <- function(com){
  Alist <- c()
  Slist <- c()
  
  num.ponds <- c(1, 2, 4, 6, 8, 12, 16, 24, 32, 42, 48)
  for (i in num.ponds) {   
    areas <- c() # hold iterated area values 
    Ss <- c() # hold iterated S values
    
    for(j in 2:iterations){
      pond.sample <- sample(51, replace = FALSE, size = i)
      pond.sample <- pond.sample[pond.sample != 29]
      pond.sample <- pond.sample[pond.sample != 31]
      
      area <- 0
      cum.abs <- vector(length = length(com[1, ]))
      
      for (k in pond.sample) { # Loop through each randomly drawn pond
        area <- area + pond.areas[k] # aggregating area
        cum.abs <- cum.abs + com[k, ]
      } # End random pond samples loop
      
      Ss <- c(Ss, length(cum.abs[cum.abs > 0]))
      areas <- c(areas, area)
    }
    
    Alist <- rbind(Alist, mean(areas))
    Slist <- rbind(Slist, mean(Ss))
    #print(c(mean(areas), mean(Ss)))
  }
  
  return(cbind(log10(Alist), log10(Slist)))
  #return(cbind(Alist, Slist))
}