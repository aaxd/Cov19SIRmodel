library(shiny)
library(deSolve)

# Define UI ----
ui <- fluidPage(
  
  plotOutput(outputId = "distPlot"),
  #titlePanel("title panel"),
  
fluidRow(
    column(4,
      
      # Input: Slider for the number of bins ----
           sliderInput(inputId = "droptime",
              label = "Days for trans to drop 0.2 to lower rate:",
              min = 10,
              max = 90,
              value = 15),
      
      # Slider for the lower rate of transmission
      sliderInput(inputId = "sectrans",
                  label = "Lower rate of transmission:",
                  min = 0.01,
                  max = 0.05,
                  value = 0.01,
                  step=0.01),
      
     
    ),
       # Slider for the initial Pop size infected
    column(4,  
     sliderInput(inputId = "InfectS",
                    label = "Initial Infected Pop size:",
                    min = 10,
                    max = 10000,
                    value = 1000,
                    step=10),
     sliderInput(inputId = "population",
                    label = "Population size:",
                    min = 1e+6,
                    max = 350e+6,
                    value = 5e+6,
                    step=1e+6),
    ),
    column(4,
           # Input: Selector for choosing dataset ----
           selectInput(inputId = "dataset",
                       label = "Choose a dataset:",
                       choices = c("Wisconsin", "Portugal")),
      )
    )
)

# Define server logic ----
server <- function(input, output) {
  require(deSolve)
  
  
  #Differential equation model. The code below requires this model
  
  # rate functions 
  	Model<-function(times,y,parms){
  	dS.dt<-(-(parms[1]*y[1]*y[2]))  		 #Susceptible: dS/dt = -(trans*Suscept*Infec)   
  	dI.dt<-((parms[1]*y[1]*y[2])-(parms[2]*y[2])) 	#Infected: dI/dt =  (trans*Suscept*Infec) - Recov*Infect
  	dR.dt<-parms[2]*y[2] 				 #Recovered/dead: dR/dt = Recov*Infect
   list(c(dS.dt,dI.dt,dR.dt))
   }
 
  	
  	 # reading real data
  	#if(input$dataset=="Wisconsin"){ dataset<-"coronaWI.txt" }
  	#if(input$dataset=="Portugal"){ dataset<-"CoronaPortugal.txt" }
  	
  	# Return the requested dataset ----
  	datasetInput <- reactive({
  	  switch(input$dataset,
  	         
  	         "Portugal" = "CoronaPortugal.txt",
  	         "Wisconsin" = "coronaWI.txt")
  	         
  	})
  	
  
  

  output$distPlot <- renderPlot({
    #reading observed data
    regionFile<-datasetInput()
    scanD<-scan(paste0("./Data/",regionFile),what=character())
    scanD2<-gsub(pattern=",",replacement="",scanD[seq(1,length(scanD),2)])
    numCols<-ifelse(regionFile=="coronaWI.txt",8,11)
    scanD3<-matrix(scanD2,byrow=T,ncol=numCols)
    if(regionFile=="CoronaPortugal.txt"){ scanD3<-scanD3[,-8]}
    
    readDF<-data.frame(date=as.Date(scanD3[,1],"%m/%d/%y"),
                       cases=as.numeric(scanD3[,4]),
                       casesToday=as.numeric(scanD3[,5]),
                       deaths=as.numeric(scanD3[,6]),
                       deathsToday=as.numeric(scanD3[,7]),
                       active=as.numeric(scanD3[,8])
    )
    readDF<-readDF[!duplicated(readDF$date),]
    readDF$casesToday<-c(NA,diff(readDF$cases))
    readDF$deathsToday<-c(NA,diff(readDF$deaths))
    readDF$TransRate<-c(NA,readDF$casesToday[2:nrow(readDF)]/readDF$active[1:(nrow(readDF)-1)])
  
    
    sdate<-ifelse(regionFile=="coronaWI.txt",as.Date("2020/3/29"),as.Date("2020/3/21"))
    sdateI<-which(readDF$date==sdate)
   
     # vector of Trans
    Total.Time<-300   			# Total days modeled
    drop.Time<-input$droptime 				# how fast does the transmission rate falls from 20% to 1% SOCIAL ISOLATION EFFICIENCY!!!
    remaining.Time<-Total.Time-drop.Time	# how much time left after the drop in transmission rate
    secTransR<-input$sectrans
    Tvec<-c(seq(0.2,secTransR,length.out=drop.Time),seq(secTransR,0,length.out=remaining.Time))
    t<-c(1,2)
    
    
    #Starting values
    recov<-.03       # recovery rate, this one is hard to estimate, I used some data from Wiscosin and China, they seem similar, not sure I am doing this part right
    Pop<-input$population
    Istart<-input$InfectS/Pop #1e+3/WiPop #InfectS/5e+6  # the infected population size in Wisconsin on March 30
    Rstart<-0          # The recovered or dead pop in Wisconsin, not exactly zero but close
    #t<-seq(0,300,1)
    Sstart<-1-Istart
    
    out.DF<-NULL
    
    for(time in 1:Total.Time){         # This loop makes the parameters of the model change with Time. You can play with the variable drop.
                                       #Time above to set how long in days the virus transmission takes to come down from 0.2 to 0.01 transmission
      
      trans<-Tvec[time]
      #ode(y,times,func,parms,method)
      
      out<-ode(y=c(S=Sstart,I=Istart,R=Rstart),
               times=t,
               func= Model,
               parms=c(trans,recov)
      )
      
      Istart<-out[2,3]
      Rstart<-out[2,4]
      Sstart<-out[2,2]
      trans<-Tvec[time]
      res<-c(out[2,],trans)
      out.DF<-rbind(out.DF,res)
      
    }
    
    colnames(out.DF)<-c("time","S","I","R","Trate")
    out.DF<-data.frame(out.DF)     # This table contains the results of the simulation
    if(regionFile=="coronaWI.txt"){
        out.DF[,1]<-seq(as.Date("2020/3/29"),by="day",length.out=Total.Time)}else{
        out.DF[,1]<-seq(as.Date("2020/3/21"),by="day",length.out=Total.Time)}
    
    MaxYval<-c(max(readDF$active),max(out.DF[,4]*Pop),max(out.DF[,3]*Pop)) 
    ylimVal<-MaxYval[which.max(MaxYval)]
    par(mar=c(5.1, 6.1, 4.1, 4.1))
    par(cex=1.2)
    plot(out.DF$time,out.DF[,4]*Pop,lwd=3,type="l",main="Covid 19 Model",col="darkblue",yaxt="n",bty="n",xaxt="n",ylab="",xlab="date",ylim=c(0,ylimVal))
    axis(1, at=seq(from=as.Date("2020/3/21"),to=as.Date("2021-01-22"),by="month"),labels=seq(from=as.Date("2020/3/21"),to=as.Date("2021-01-22"),by="month"))
    #axis(2,pos=as.Date("2020/3/15"),col="darkblue",col.axis="darkblue")
    
    #par(new=T)
   # plot(out.DF$time,out.DF[,3]*Pop,lwd=3,bty="n",type="l",axes=F,xlab="",ylab="",ylim=c(0,max(out.DF$time,out.DF[,4]*Pop)))
    lines(out.DF$time,out.DF[,3]*Pop,lwd=3)
    axis(2,pos=sdate)
    grid()
    #if() {# switch to print real data or not 
    #real data
    if(regionFile=="coronaWI.txt"){
     lines(x=seq(as.Date("2020/3/29"),readDF$date[nrow(readDF)],by="day"),y=readDF$active[sdateI:nrow(readDF)],lwd=3,col="brown1")}else{
     lines(x=seq(as.Date("2020/3/21"),readDF$date[nrow(readDF)],by="day"),y=readDF$active[sdateI:nrow(readDF)],lwd=3,col="brown1")
        }
    
    # Plot of susceptible in model
    #par(new=T)
   # plot(seq(as.Date("2020/3/29"),by="day",length.out=Total.Time),out.DF[,2]*Pop,lwd=3,type="l",col="grey75",bty="n",ylab="",xlab="",yaxt="n",xaxt="n")
    #axis(4,pos=as.Date("2021/2/1"),col="grey75",col.axis="grey50")
    
    
    par(new=T)#seq(as.Date("2020/3/29"),by="day",length.out=Total.Time)
    plot(out.DF$time,out.DF[,5],lwd=3,lty=3,col="coral4",axes=F,xlab="",ylab="",type="l")
    axis(4,col.axis="coral4",pos=as.Date("2021/1/23"))
    if(regionFile=="coronaWI.txt"){
         lines(seq(as.Date("2020/3/29"),readDF$date[nrow(readDF)],by="day"),readDF$TransRate[sdateI:nrow(readDF)],col="brown1",lty=3,lwd=3)}else{
         lines(seq(as.Date("2020/3/21"),readDF$date[nrow(readDF)],by="day"),readDF$TransRate[sdateI:nrow(readDF)],col="brown1",lty=3,lwd=3)
        }
    
    legend(x=18570,y=0.165,legend=c("Infected Model","Infected Real","Recov/Dead","Rate of Trans model","Rate of Trans real"),
           lty=c(1,1,1,3,3),col=c("black","brown1","darkblue","coral4","brown1"),lwd=3)
    
    
    
     # x    <- faithful$waiting
   # bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    #hist(x, breaks = bins, col = "#75AADB", border = "orange",
    #     xlab = "Waiting time to next eruption (in mins)",
     #    main = "Histogram of waiting times")
    
  })
  
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
