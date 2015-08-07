################################################################################
#                                                                              #
# Functions for calculating metrics of diversity, evenness, rarity, etc.       #
# These are not included in other diversity packages, e.g., Vegan              #
#                                                                              #
################################################################################
#                                                                              #
# Written by: Ken Locey                                                        #
#                                                                              #
################################################################################
#                                                                              #
# Recent Changes:                                                              #
#                                                                              #
# Future Changes (To-Do List):                                                 #
#         1.                                                                   #
#         2. Add warnings                                                      #
#                                                                              #
################################################################################

#### A function to generate observed richness
S.obs <- function(x = ""){ rowSums(x > 0) * 1}

