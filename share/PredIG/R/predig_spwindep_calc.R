#### PREDIG SPW INDEP (SPW 25)

### LIBRARIES
library(dplyr, quietly = T)
library(xgboost, quietly = T)
library(argparser, quietly = T)

### ARGUMENTS
# Create a parser
p <- arg_parser("PREDIG INDEP SPW25")

# Add command line arguments
p <- add_argument(p, "--input", help= "Path to Predig Feature dataset (.csv format)", type="character", default = "../datasets/predig_spwindep_model_features.csv")
p <- add_argument(p, "--seed", help= "Seed number", type = "numeric",  default= 123)
p <- add_argument(p, "--model", help= "Path to PredIG pretrained XGBoost model. File extension '.model'", type = "character",  default= ".")

# Parse the command line arguments
argv <- parse_args(p)

### LOAD MODEL
predig_model <- xgb.load(argv$model)

## PREDICT based on TT - SPW Indep
predig_input <- read.csv(argv$input)
predig_input <- as.matrix(predig_input)
predig_score = predict(predig_model, predig_input)
predig_df <- cbind(predig_input, predig_score)
predig_df <- as.data.frame(predig_df)

### Export predig test
write.csv(predig_df, file = paste(gsub(".csv", "", argv$input),"_predig.csv", sep = ""), quote = F, row.names = F)
