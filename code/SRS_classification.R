library(SepstratifieR)

test_data <- read.csv(file = './pysudo_Davenport_logcpm.csv',row.names = 1)# data is pysudo bulk, logcpm

predictions <- stratifyPatients(test_data)

write.csv(predictions@SRS, "Davenport_predictions.csv", row.names = TRUE)

