
#pdf("Sims_estimating_mu_over_s_DIFF_s_values.pdf")
pdf("TRY.pdf")

for (Ne in c(10000)){
#currently Ne cannot be changed in the sims 
#DataOverview<-data.frame("N"=Ne,"mu"=0,"cost"=0,"num_runs"=0,"datapointsperrun"=0,"equiPi"=0,"VarPi"=0,"expectedPi"=0,"t_half"=0)
	n=0
	NUMRUNS=1; numoutputs=10000
	system("./Code_and_shellscript/make_HIV1site")	#compile the code 
for (mu in c(0.000002,0.00002,0.0002)){
#	for (mu in c(0.00002)){
#for (cost in c(0.001,0.005,0.01,0.05,0.1,0.2)){
		for (cost in c(0.01,0.05,0.1)){
			print("")
			print(paste("cost",cost))
			print("")
			print(paste("mu",mu))
			print("")
			avepivalues<-vector()	
			for (seed in 1:NUMRUNS){
#make script 
				x<-"#!/bin/bash"
				x<-c(x,paste("mu=",mu,sep=""))
				x<-c(x,paste("cost=",cost,sep=""))
				outputfrequency=min(c(2*Ne,ceiling(5/cost)))
				x<-c(x,paste("output_every_Xgen=",outputfrequency,sep=""))
				x<-c(x,paste("numgen_inN=",(numoutputs+2)*outputfrequency/Ne,sep=""))
				x<-c(x,paste("start_output=",2*outputfrequency/Ne,sep=""))
				x<-c(x,paste("for seed in",seed))
				x<-c(x,"do",
					 "echo \"", "$seed", "$mu", "$cost",
					 "$output_every_Xgen", "$numgen_inN", "$start_output",
					 paste("\" | ./Code_and_shellscript/HIVevolution_HIV1site >../Data/Link",cost,"_",seed,".txt",sep=""), 
					 "done")
				write(x,file="./Code_and_shellscript/tempscript.sh")
				system("chmod 775 ./Code_and_shellscript/tempscript.sh")
#Run tempscript.sh
				system("./Code_and_shellscript/tempscript.sh")
				
#READ FREQS FILE
				read.csv(paste("../Data/Link",cost,"_",seed,".txt",sep=""),sep="\t",header=TRUE)->simdata
#in stead of the real frequencies, lets assume we have a sample from each patient of, say, 100, seqs. 
				plot(c(0,0),col=0,xlim=c(-0.5,3.5),ylim=c(-1*mu/cost,8*mu/cost),xlab="Num Patients",xaxt="n",ylab="95% of observed average freq of allele",main=paste("Ne",Ne,", mu",mu,", Theta",2*Ne*mu ,", cost",cost))
				axis(1, at=log10(samplesizes), labels=samplesizes)
				abline(h=mu/cost,lty=2)
				
				diffnumseqs <- c(50,100,200)
				for (num_seqs_per_patient in diffnumseqs){
				for (i in 1:length(simdata$freq)){
					simdata$est_freq[i]<-rbinom(1,num_seqs_per_patient,simdata$freq[i])/num_seqs_per_patient}
				
#system(paste("rm ","../Data/Link",cost,"_",seed,".txt",sep=""))
				samplesizes<-c(1,3,10,30,100,300,1000)
				for (num_patients in samplesizes){
					list_averages<-vector()
					for (i in 1:1000){ 
						list_averages<-c(list_averages,mean(sample(simdata$est_freq,num_patients)))}
					co=which(diffnumseqs==num_seqs_per_patient)
					X=(which(diffnumseqs==num_seqs_per_patient)-2)*0.1
					print(paste(num_seqs_per_patient,num_patients,log10(num_patients)-0.02+X))
					rect(log10(num_patients)-0.02+X,sort(list_averages)[25],log10(num_patients)+0.02 +X,sort(list_averages)[975],col=co)
#	text(log10(num_patients)+X,sort(list_averages)[975]+0.15*mu/cost+(X+0.2)/200,paste(round(sort(list_averages)[25]/(mu/cost),2),"-",round(sort(list_averages)[975]/(mu/cost),2)),cex=0.6)
				
					text(log10(num_patients)+X,(-1+co/8)*mu/cost,paste(round(sort(list_averages)[25]/(mu/cost),2),"-",round(sort(list_averages)[975]/(mu/cost),2)),cex=0.5)	
				}
				text(log10(num_patients),6.5*mu/cost,"bars and",cex=0.8)
				text(log10(num_patients),6*mu/cost,"numbers indicate range",cex=0.8)
				text(log10(num_patients),5.5*mu/cost,"of 95% of estimates",cex=0.8)
				text(log10(num_patients),5.*mu/cost,"black:50, red: 100, ",cex=0.8)
				text(log10(num_patients),4.5*mu/cost,"green: 200 sequences/pat",cex=0.8)
				text(-0.3,1.1*mu/cost,"mu/s",cex=0.8)
				
				}}}}}

dev.off()

#rbinom(n, size, prob)
