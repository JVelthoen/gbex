### script to generate a txt file with all parameter values that need to be used
setwd("/Users/jjvelthoen/Documents/GitHub/gbex/documentation/Examples/step_function.png/cluster run/normal_noise/")
file_name <- "parameters.txt"

file.create(file_name)

n = c(500,1000,2000)
tau_thresh = c(0.8,0.9)
reps <- 20
lines_to_write = ""
for(ncnt in n){
  for(tcnt in tau_thresh){
    for(i in 1:reps){
      lines_to_write = paste0(lines_to_write,paste0("n = ",ncnt,"\n","tau_thresh = ",tcnt,"\n"))
    }
  }
}


write(lines_to_write,file=file_name)

nr_simulations <- reps*length(tau_thresh)*length(n)
print(paste("Total number of simulations",nr_simulations))
