
shiny_ui <- 
  fluidPage(tags$head(tags$style(HTML("body { max-width: 1250px !important; }"))),
            titlePanel("NMCSDE v0.0"),
            sidebarLayout(
              sidebarPanel(width=3, tags$head(tags$style(type="text/css", ".well { max-width: 300px; }")),
                           sliderInput("deg", HTML("Degree of B-spline Basis:"), min = 0, max = 5, value = 3, step = 1, width="250px"),
                           uiOutput("df", width="250px"), 
                           tags$hr(style="border-color: red;", width="150px"),
                           column(12,sliderInput("dimn", HTML("# of Common Basis <br/>(# of Clusters for Clustering):"), min = 1, max = 10, value = 3, step = 1, width="210px")),
                           column(8,radioButtons('lambda', "Tuning Parameter:", c("S" = "sing", "M" = "mult", "F"="fix"), inline = TRUE),
                                  checkboxInput('cbs', 'Comm.Bas.Spec.', value=FALSE)),
                           column(4,sliderInput('init.lam', "2^i, i:", value = -55, min = -55, max = 45, step = 1, width="150px")),
                           column(12,checkboxGroupInput('model.est', NULL, choices=c("MNCSDE"="mncsde","Separate Est."="sep"), selected=c("mncsde","sep"), inline = TRUE, width="250px")), 
                           tags$hr(style="border-color: red;", width="150px"),
                           checkboxInput('log.psd', 'Log Scale "Y" : (SD-read)', value=FALSE, width="400px"),
                           checkboxInput('clust.Asd', 'Clustering: Scores vs SDs', value=TRUE, width="400px"), tags$hr(style="border-color: red;", width="150px"),
                           uiOutput("sh.tr", width="400px"), checkboxInput('sh.all', 'Show All Spec. (Var)', value=FALSE, width="400px")),
              mainPanel(width=9, tags$style(type="text/css", ".shiny-output-error { visibility: hidden; }",".shiny-output-error:before { visibility: hidden; },"),#".nav-tabs {font-size: 10px}"),
                        tabsetPanel(id="Panel", type = "tabs",
                                    tabPanel("Data",
                                             column(12,uiOutput("ts.selected", align = "center"), style="color:red;"),
                                             fluidRow(column(4,radioButtons('f.choice', 'Choose a file:', c("Server" = "server", "Upload" = "upload", "Simulation" = "sim"), selected= "sim", inline = TRUE, width="250px")),
                                                      column(2,uiOutput("which.mts")), column(3,uiOutput("s.choice"),uiOutput("file")), column(3, uiOutput("freq"))),
                                             conditionalPanel(condition="output.flag_sim",
                                                              column(4,checkboxGroupInput('model', 'Models:', choices=c("Model I"="1","Model II"="2","Model III"="3"), selected=c("1","2","3"), inline = TRUE, width="250px")),
                                                              column(3,sliderInput("t.len","Length of TS", min = 100, max = 400, value = 200, width="250px")),
                                                              column(3,sliderInput("rep","Repetition", min = 2, max = 50, value = 4, width="250px")),
                                                              column(2,sliderInput("rho","Correlation", min = -0.99, max = 0.99, value = 0.9, step = 0.01, width="250px"))),
                                             column(12,tableOutput('data'))),
                                    tabPanel("Data Description",
                                             column(3,uiOutput("desc", width="250px")), column(3,uiOutput("which.var")), column(3,uiOutput("ams.choice", width="400px")), column(3,uiOutput("sts.choice")),
                                             plotOutput("data.desc", height = 600, width = 600)),
                                    tabPanel("Basis Functions",
                                             column(8,plotOutput("basis.desc", height = 600, width = 600)), column(4,uiOutput("basis.n", width="300px"))),
                                    tabPanel("NMCSDE",
                                             column(3,uiOutput("n.iter", width="200px")), column(3,uiOutput("which.spec", width="200px")), column(2,uiOutput("runit", width="200px")), column(4, fluidRow(column(6,uiOutput("diffP")), column(6,uiOutput("m"))),tags$hr(style="border-color: red;", width="200px", align="center")),
                                             column(8,plotOutput("MNCSDE.out", height = 600, width = 600)),
                                             column(4,uiOutput("results"),uiOutput("s.iter"),uiOutput("pair.dend.dens"),uiOutput("ams.MNCSDE"),uiOutput("sts.MNCSDE"),uiOutput("download.b"))),
                                    tabPanel("Manual", includeMarkdown(system.file("shiny/rmd", "report.Rmd", package = "NMSDE")))))
            )
  )

shiny_server <- function(input, output, session) {
  
  data("servshiny"); desc.choice <- list("Elbow Method"="elbow", Data=c("Time Series" = "ts"));
  options(shiny.maxRequestSize=100*1024^2); iX <- reactiveVal(list()); imodelPar <- reactiveVal(NULL); Ts <- reactiveVal(NULL); track.data <- reactiveVal(list());
  output$flag_sim <- reactive("sim"%in%input$f.choice); outputOptions(output, "flag_sim", suspendWhenHidden = FALSE)
  
  simul <- function() {
    library(dse)
    sigma<-matrix(c(1,input$rho,input$rho,1), 2, 2); set.seed(5)
    Phi1<-list(diag(2), -diag(c(0.5,-0.3)), -diag(c(0,-0.5)))
    Phi2<-list(diag(2), -matrix(c(0,-0.6,0.5,-0.2), 2, 2))
    Phi3<-list(diag(2), -matrix(c(0.7,0.2,0.7,0.2), 2, 2))
    Phis<-list(Phi1, Phi2, Phi3); ARcoeffs <- models <- X <- list()
    for (k in 1:length(Phis)) {
      ARcoeffs[[k]] <- array(0,dim=c((length(Phis[[k]])), 2, 2))
      for(i in 1:2) for(j in 1:(length(Phis[[k]]))) {ARcoeffs[[k]][j,,i]<-Phis[[k]][[j]][,i]}
      models[[k]] <- ARMA(A=ARcoeffs[[k]], B=diag(2), TREND=rep(0,2))
    }
    idx <- as.numeric(sample(rep(input$model,input$rep)))
    for (i in 1:length(idx)) {
      varsim <- simulate(models[[idx[i]]], sampleT=input$t.len, Cov=sigma)
      X[[i]] <- matrix(varsim$output, nr=2, nc=input$t.len, byrow=T)
    }
    return(list(X=X, modelPar=list(paras=Phis, sigma=sigma, rnd=idx)))
  }
  
  output$df <- renderUI({
    sliderInput("df", HTML("Deg. of freedom of B-spline Basis:"), min = input$deg+1, max = min(40,input$freq[2]-input$freq[1]), value = 20, step = 1, width="250px")
  })
  
  output$sh.tr <- renderUI({
    if (is.null(imodelPar())) return(NULL)
    checkboxInput('sh.tr', 'Show True Spec.', value=FALSE, width="400px")
  })
  
  observeEvent(input$model, {  updateSliderInput(session, "dimn", value=length(input$model))  })
  
  output$which.mts <- renderUI({
    if (input$f.choice=="upload" && is.null(input$file)) return(NULL)
    m.choices <- 1:length(iX()); names(m.choices) <- names(iX())
    selectInput('which.mts', "Which MTS", choices=m.choices, width="125px")
  })
  
  output$s.choice <- renderUI({
    if (input$f.choice!="server") return();
    s.choices <- 1:length(dat); names(s.choices) <- names(dat)
    selectInput("s.choice","Select a file from server: ", choices = s.choices, width="250px");
  })
  
  output$file <- renderUI({
    if (input$f.choice!="upload") return();
    fileInput('file', 'Choose .Rdata File', accept=c('Rda/Rdata', 'Rda/Rdata', '.Rda'))
  })
  
  output$freq = renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || !length(iX())) return(NULL);
    freq.l<-floor((ncol(iX()[[1]])-1)/2)
    sliderInput("freq", "Frequency (Order):", min = 1, max = freq.l, value = c(1,freq.l), dragRange=TRUE, step = 1, width="400px")
  })
  
  output$ts.selected = renderText({
    if (input$f.choice=="server" && is.null(input$s.choice)) return()
    if (input$f.choice=="upload" && is.null(input$file)) return("<b>Select a '.Rda' file that contain the multivariate time series as a list in rows</b>")
    if (input$f.choice=="upload") {
      load(input$file$datapath);  iX(X); imodelPar(NULL)
    } else if (input$f.choice=="sim") {
      sim <- simul(); iX(sim$X); imodelPar(sim$modelPar)
    } else {
      i <- as.numeric(input$s.choice); iX(dat[[i]]); tr.par <- datPar[[i]]
      if (is.null(unlist(tr.par))) imodelPar(NULL) else imodelPar(tr.par)
    }
    updateSliderInput(session, "dimn", max=min(10,length(iX())))
    text <- paste("<b>",length(iX()),"Multivariate Time series of length",ncol(iX()[[1]]),"</b>"); return(text)
  })
  
  output$data <- renderTable({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$which.mts)) return(NULL)
    return(head(t(iX()[[as.numeric((input$which.mts))]]),15))
  })
  
  output$desc <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file))) return();
    selectInput("desc","Description: ", choices = desc.choice, selected = "ts", width="250px")
  })
  
  output$ams.choice <- renderUI({
    if (is.null(input$desc)) return(NULL)
    if (input$f.choice=="upload" && is.null(input$file) || input$desc%in%c("scree","elbow")) return(NULL)
    radioButtons('ams.choice', 'Plot Choices:', c("All" = "all", "Multiple" = "multiple", "Single" = "single"), selected= "all", inline = TRUE, width="400px")
  })
  
  output$sts.choice <- renderUI({
    if (is.null(input$desc) || is.null(input$ams.choice)) return(NULL)
    if (!(input$f.choice=="upload" && is.null(input$file)) && (!input$desc%in%c("scree","elbow") && input$ams.choice!="all"))
      if (input$ams.choice=="single") {sliderInput("sts.choice", "Choose TS:", min = 1, max = length(iX()), value = 1, step = 1, width="400px")} else {
        sliderInput("sts.choice", "Choose a range( '< 10' ):", min = 1, max = length(iX()), value = c(1,min(10,length(iX()))), dragRange=TRUE, step = 1, width="400px")
      }
  })
  
  output$which.var <- renderUI({
    if (input$f.choice=="upload" && is.null(input$file) || is.null(input$desc)) return(NULL)
    v.choices <- 1:nrow(iX()[[1]]); names(v.choices) <- rownames(iX()[[1]])
    selectInput('which.var', "Which Var", choices=v.choices, width="125px")
  })
  observeEvent(input$which.var, {
    ts <- rapply(iX(), classes = 'matrix', how = 'list', f = function(x) x[as.numeric(input$which.var),, drop = FALSE])
    Ts(t(do.call(rbind,ts)))
  })
  
  output$data.desc = renderPlot({
    if (input$f.choice=="upload" && is.null(input$file) || is.null(input$which.var)) return(NULL)
    if (input$f.choice=="upload") {fname <- input$file$name} else {fname <- names(iX())};
    indx <- as.numeric(input$sts.choice); par(mgp=c(1.5,.5,0),mar=c(3,3,1.5,1.75),cex=1.5);
    withProgress(message = 'Please wait......', value = 0, {
      if (input$desc == "ts") {
        clust <- hclust(dist(t(Ts())),method="ward.D"); clcol <- cutree(clust,k=input$dimn);
        if (input$ams.choice=="multiple") {
          if (length(input$sts.choice)==1) return();
          plot.ts(Ts()[,seq(indx[1],indx[2])],ann=F,ylim=range(Ts()));
          title(main=paste("Time Series -",fname,"- from TS:",indx[1],"to TS:",indx[2]), xlab="Time")
        } else {
          if (input$ams.choice=="all" || is.null(input$sts.choice)) {indx <- 1:ncol(Ts()); m.lab <- fname} else {m.lab <- paste(fname,"- TS:",indx)}
          ts.plot(Ts()[,indx], col=clcol[indx], main=paste("Time Series -",m.lab), ylab="", ylim=range(Ts()), gpars=list(xaxt="n"));
          axis(1,trunc(summary(1:nrow(Ts()))[-4]))
        }
      } else {
        isolate(init <- initial(X=iX(), nbasis=input$df, norder=input$deg+1, R=input$dimn, pen="diff", a=2, theta_P2=!input$cbs, tru=imodelPar(), f.r=input$freq))
        if (!input$desc%in%c("scree","elbow")) {
          if (input$desc == "pgram") {resp <- init$I}
          if (input$desc == "s.pgram") {resp <- exp(init$B%*%init$Theta.tA)}
          if (input$desc == "tsvd.pgram") {resp <- svd(log(init$I)); resp <- exp(resp$u[,1:input$dimn]%*%diag(resp$d[1:input$dimn],input$dimn)%*%t(resp$v[,1:input$dimn]))}
          if (input$desc == "s.tsvd.pgram") {resp <- init$svd.Theta.tA; resp <- exp(init$B%*%resp$u[,1:input$dimn]%*%diag(resp$d[1:input$dimn],input$dimn)%*%t(resp$v[,1:input$dimn]))}
          clust <- hclust(dist(resp),method="ward.D"); clcol <- cutree(clust,k=input$dimn);
          if (input$ams.choice=="multiple") {
            if (length(input$sts.choice)==1) return();
            plot.ts(resp[,seq(indx[1],indx[2])],ann=F,axes=FALSE,frame=TRUE,ylim=range(resp), log=ifelse(input$log.psd,"y",""));
            title(main=paste("Periodogram -",fname,"- from P:",indx[1],"to P:",indx[2]), xlab="Frequency")
          } else {
            if (input$ams.choice=="all" || is.null(input$sts.choice)) {indx <- 1:ncol(Ts()); m.lab <- fname} else {m.lab <- paste(fname,"- P:",indx)}
            ts.plot(resp[,indx], col=clcol[indx], main=paste("Periodogram -",m.lab), xlab="Frequency", ylab="", ylim=range(resp), log=ifelse(input$log.psd,"y",""), gpars=list(xaxt="n"));
            axis(1,summary(1:nrow(resp))[-4],summary(init$freq)[-4])#c(0,.125,.25,.375,.5))
          }
        } else {
          if (input$desc == "elbow") plot(init$elbow, size=5, cex=4)
          else{
            percentage <- paste0("\n The variation explained in first ",input$dimn," adaptive basis of the initial est. is ", sum(init$scree[1:input$dimn]),"%")
            plot(init$scree, xlab="Component", ylab="Percentage", main=paste("Scree Plot -",fname,percentage), pch=20, frame.plot=F, ylim=c(0,100), cex.main=.75)
            points(init$scree[1:input$dimn],col=2)
          }
        }
      }})
  })
  
  output$basis.desc = renderPlot({
    if (is.null(input$basis.n)) return()
    B <- fda::bsplineS(seq(0,0.5,length.out=1000),breaks=seq(0,.5,length.out=input$df-input$deg+1),norder=input$deg+1)
    ts.plot(B, col=8, main="Bspline Basis", xlab="Frequency", gpars=list(xaxt="n"))
    points(B[,input$basis.n], type="l", lwd=2, col=2)
    axis(1,summary(1:nrow(B))[-4],c(0,.125,.25,.375,.5))
  })
  
  output$basis.n <- renderUI({
    sliderInput("basis.n", "Basis #:", min = 1, max = input$df, value = 1, step=1, width="400px")
  })
  
  output$n.iter <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file))) return();
    numericInput('n.iter',"Maximum # of iterations:", value = 100, min = 1, max = 100, width="200px")
  })
  
  output$which.spec <- renderUI({
    if (input$f.choice=="upload" && is.null(input$file)) return(NULL)
    s.choices <- 1:nrow(iX()[[1]])^2;# names(s.choices) <- rownames(iX()[[1]])
    selectInput('which.spec', "Which Spec", choices=s.choices, width="125px")
  })
  
  output$runit <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file))) return();
    if (input$f.choice!="upload") {fname <- names(iX())} else {fname <- input$file$name};
    actionButton('runit', paste0('RUN (',fname,')'))
  })
  
  output$diffP <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file))) return();
    checkboxInput('diffP', 'Diff. Pen.', TRUE, width="200px")
  })
  
  output$m <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$diffP)) return(); if (!input$diffP) return()
    sliderInput('m','Pen. Diff. Ord.', value = 2, min = 0, max=min(4,input$df-1), step=1, width="200px")
  })
  
  output$results <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$runit)) return(); if (!input$runit) return()
    choice <- NULL; if("mncsde"%in%input$model.est) choice <- "Collective(MNCSDE)"
    if("sep"%in%input$model.est) choice <- c(choice,"Separate(MNSDE)","tSVD.MNSDE");
    selectInput('results', 'Select Results: ', choices = choice, width="200px")
  })
  
  output$s.iter <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$results)) return(); if (!input$runit) return()
    if (input$results%in%c("Separate(MNSDE)","tSVD.MNSDE")) return();
    if (!exists("res")) res <- list(ans=list(r.iter=input$n.iter))
    return(sliderInput("s.iter", "MNCSDE Iteration #:", min = 1, max = res$ans$r.iter, value = res$ans$r.iter-1, step = 1, width="200px"))
  })
  
  output$pair.dend.dens <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$runit)) return(); if (!input$runit) return()
    epd.dens <- list("Estimated Densities"="Densities", Clustering = c("Scores"="Scores", "Dendrogram"="Dendrogram"));
    selectInput('pair.dend.dens', 'Select Plot: ', choices = epd.dens, selected="Scores", width="200px")
  })
  
  runit <- eventReactive(input$runit, {
    if (input$f.choice!="upload") {fname <- names(iX())} else {fname <- input$file$name};
    withProgress(message = 'MNCSDE: Newton-Raphson', value = 0, {
      lam <- 2^input$init.lam-2^(-55)
      init <- initial(X=iX(), nbasis=input$df, norder=input$deg+1, R=input$dimn, pen=ifelse(input$diffP,"diff","deriv"), a=ifelse(input$diffP,input$m,2), theta_P2=!input$cbs, tru=imodelPar(), f.r=input$freq)
      if (input$cbs && input$lambda=="mult") lam <- rep(lam,input$dimn)
      else if (!input$cbs && input$lambda=="mult") lam <- rep(lam,nrow(iX()[[1]])^2) #matrix(rep(lam,input$dimn*nrow(iX()[[1]])^2),nr=input$dimn)
      # else if (!input$cbs && input$lambda=="sing") lam <- rep(lam,nrow(iX()[[1]])^2)
      ans <- ansi <- anst <- list();
      if ("mncsde" %in% input$model.est) {
        try(ans <- iterate(init=init, lambda=lam, n.iter=input$n.iter, update.lam=(input$lambda!="fix"))); #m=ifelse(input$diffP,input$m,2), eps=0.1
        ans$dens <- aperm(array(apply(tA2G(ans$thetas[[ans$r.iter]],ans$As[[ans$r.iter]],init),3,inChol2spec,init=init),dim=c(init$P^2,init$K,init$m)),c(2,3,1))
        updateSliderInput(session, "s.iter", max = ans$r.iter, value = ans$r.iter-1);
      }
      if ("sep" %in% input$model.est) {
        ansi$dens <- anst$dens <- array(0,dim=c(input$K,length(iX()),init$P^2));
        ansi.Theta.tA <- array(0,dim=c(input$df,length(iX()),init$P^2));
        ansi.stA <- array(0,dim=c(init$P^2,init$K,length(iX())));
        anst$As <- array(0,dim=c(length(iX()),input$dimn,init$P^2));
        #if (input$cbs && input$lambda=="mult") lam <- lam[1]
        for (i in 1:length(iX())) {
          try(ansi[[i]] <- iterate(init=init, lambda=lam, n.iter=input$n.iter, update.lam=(input$lambda!="fix"), which.spec=i))   ##### Separate Estimations
          for (j in 1:init$P^2) {
            if (input$cbs) {
              ansi.Theta.tA[,i,j] <- ansi[[i]]$thetas[[ansi[[i]]$r.iter]]%*%t(ansi[[i]]$As[[ansi[[i]]$r.iter]][,,j])
            } else {
              ansi.Theta.tA[,i,j] <- ansi[[i]]$thetas[[ansi[[i]]$r.iter]][,,j]%*%t(ansi[[i]]$As[[ansi[[i]]$r.iter]][,,j])
            }
          }
          if (input$cbs) {
            svd.t <- svd(ansi.Theta.tA[,i,])
            ansi.Theta.tA[,i,] <- svd.t$u[,1]%*%diag(svd.t$d[1],1)%*%t(svd.t$v[,1])
          }
          ansi.stA[,,i] <- t(init$B%*%ansi.Theta.tA[,i,])
        }
        ansi$dens <- aperm(array(apply(ansi.stA,3,inChol2spec,init=init) ,dim=c(init$P^2,init$K,init$m)),c(2,3,1))
        if (input$cbs) {
          merg.theta.tA <- matrix(ansi.Theta.tA,nr=input$df)
          svd.merg.theta.tA<-svd(merg.theta.tA)
          thetas<-as.matrix(svd.merg.theta.tA$u[,1:input$dimn])
          mergA<-as.matrix(svd.merg.theta.tA$v%*%diag(svd.merg.theta.tA$d,ncol(svd.merg.theta.tA$v))[,1:input$dimn])
          for(j in 1:init$P^2) anst$As[,,j]<-mergA[((j-1)*length(iX())+1):(j*length(iX())),]
        } else {
          svd.P2<- list(); thetas <- array(0,dim=c(input$df,input$dimn,init$P^2))
          for (j in 1:init$P^2) {
            svd.P2[[j]] <- svd(ansi.Theta.tA[,,j]); thetas[,,j] <- svd.P2[[j]]$u[,1:input$dimn]
            anst$As[,,j] <- as.matrix(svd.P2[[j]]$v%*%diag(svd.P2[[j]]$d,length(svd.P2[[j]]$d))[,1:input$dimn])
          }
        }
        anst$dens <- aperm(array(apply(tA2G(thetas,anst$As,init),3,inChol2spec,init=init),dim=c(init$P^2,init$K,init$m)),c(2,3,1))
      }
      return(list(ans=ans,ansi=ansi,anst=anst,init=init))
    })
  })
  
  output$ams.MNCSDE <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$pair.dend.dens)) return(); if (!input$runit) return()
    if (input$pair.dend.dens%in%c("EEG", "Scores", "Dendrogram")) return();
    radioButtons('ams.MNCSDE', 'Plot Choices:', c("All" = "all", "Multiple" = "multiple", "Single" = "single", "Clusters" = "clust"), selected= "all", inline = TRUE, width="200px")
  })
  
  output$sts.MNCSDE <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$ams.MNCSDE)) return(); if (!input$runit) return()
    if ((input$pair.dend.dens%in%c("EEG", "Scores", "Dendrogram")) || input$ams.MNCSDE=="all") {return()}
    if (input$ams.MNCSDE=="multiple") {sliderInput("sts.MNCSDE", "Choose a range( '< 10' ):", min = 1, max = length(iX()), value = c(1,min(10,length(iX()))), dragRange=TRUE, step = 1, width="200px")}
    else if (input$ams.MNCSDE=="single") {return(sliderInput("sts.MNCSDE", "Choose TS:", min = 1, max = length(iX()), value = 1, step = 1, width="200px"))}
    else return(sliderInput("sts.MNCSDE", "Choose Cluster:", min = 0, max = input$dimn, value = 0, step = 1, width="200px"))
  })
  
  output$download.b <- renderUI({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$runit)) return(); if (!input$runit) return()
    if (input$f.choice!="upload") {fname <- names(iX())} else {fname <- input$file$name};
    downloadButton('saveit', paste0('Download results (',fname,')'));
  })
  
  output$saveit <- downloadHandler('saveit', filename = function() {paste0(ifelse(input$f.choice!="upload",input$s.choice,input$file$name),".RData")}, content = function(file) {
    init = res$init; ans = res$ans; ansi = res$ansi;
    clust <- clusters <- NULL; clust$A <- hclust(dist(scale(ans$A[[input$s.iter]])), method="ward.D");
    clust$dens <- hclust(dist(t(init$B%*%ans$Theta[[input$s.iter]]%*%t(ans$A[[input$s.iter]]))), method="ward.D");
    clusters$A <- cutree(clust$A,k=input$dimn); clusters$dens <- cutree(clust$dens,k=input$dimn);
    save(init, ans, ansi, clust, clusters, file=file)
  })
  
  output$MNCSDE.out = renderPlot({
    if ((input$f.choice=="upload" && is.null(input$file)) || is.null(input$model.est)) return();
    if (input$f.choice!="upload") {fname <- names(iX())} else {fname <- input$file$name};
    res <<- runit()
    flag <- TRUE; if (!is.null(input$ams.MNCSDE)) if (input$ams.MNCSDE=="multiple") flag <- FALSE;
    if (input$sh.all && flag && (input$pair.dend.dens!="Scores" || input$dimn==1)) {
      js <- 1:res$init$P^2; par(mfrow=c(res$init$P,res$init$P),mgp=c(1.5,.5,0),mar=c(3,3,1.5,1.75))
    } else { js <- as.numeric(input$which.spec); par(mgp=c(1.5,.5,0),mar=c(3,3,1.5,1.75),cex=1.5) }
    if (input$results=="Collective(MNCSDE)") dens.P2 <- aperm(array(apply(tA2G(res$ans$thetas[[input$s.iter]],res$ans$As[[input$s.iter]],res$init),3,inChol2spec,init=res$init),dim=c(res$init$P^2,res$init$K,res$init$m)),c(2,3,1))
    for (j in js) {
      if (input$results=="Collective(MNCSDE)") {
        A <- as.matrix(res$ans$As[[input$s.iter]][,,j])
        dens <- dens.P2[,,j]
      } else if (input$results=="Separate(MNSDE)") {
        A <- cmdscale(dist(t(res$ansi$dens[,,j])),input$dimn);
        dens <- res$ansi$dens[,,j]
      } else {
        A <- as.matrix(res$anst$As[,,j]);
        dens <- res$anst$dens[,,j]
      }
      rownames(A) <- colnames(dens) <- colnames(Ts());
      if (input$clust.Asd) clus.inp <- scale(A) else clus.inp <- log(t(dens))
      clust <- hclust(dist(clus.inp), method="ward.D"); clcol <- cutree(clust,k=input$dimn);
      if (input$pair.dend.dens=="Scores") {
        if (ncol(A)==1) {plot(A, col=clcol, main=paste("Score Plot -",j,"-", fname), pch=20, xaxt="n", ylab="Score", xlab=""); axis(1,1:nrow(A),clust$labels,las=2)}
        else {pairs(A, col=clcol, main=paste("Matrix Plot -",j,"-", fname), cex.lab=2, label=paste("Score",1:ncol(A)))}
      } else if (input$pair.dend.dens=="Dendrogram") {
        plot(color_labels(as.dendrogram(clust), k = input$dimn, col=unique(clcol)), main=paste("Dendrogram -",j,"-",fname), cex.axis=0.75)
      } else {
        if (is.null(input$ams.MNCSDE)) return(); indx <- as.numeric(input$sts.MNCSDE);
        if (input$ams.MNCSDE=="multiple") {
          if (is.null(input$sts.MNCSDE)) return(); if (length(input$sts.MNCSDE)==1) return();
          plot.ts(dens[,seq(indx[1],indx[2])],ann=F,axes=FALSE,frame=TRUE,ylim=range(dens),log=ifelse(input$log.psd,"y",""));
          title(main=paste(input$results,"-",j,"-",fname,"- from TS:",indx[1],"to","TS:",indx[2]), xlab="Frequency")
        } else {
          range.dens <- c(ifelse(input$log.psd,min(abs(dens)),min(dens)),max(dens)); ind.clust <- FALSE;
          if (input$ams.MNCSDE=="all") {indx <- 1:length(iX()); m.lab <- fname}
          else if (input$ams.MNCSDE=="single") {m.lab <- paste("- TS:",indx)}
          else {
            m.lab <- ifelse(input$sts.MNCSDE, paste0("Cluster : ",input$sts.MNCSDE, " (", sum(clcol==input$sts.MNCSDE), " SDs)"), "Averaged Clusters")
            if (is.null(input$sts.MNCSDE)) return();
            if (!input$sts.MNCSDE[1]) {
              den <- NULL; indx <- unique(clcol)
              for (i in indx) den <- cbind(den, apply(as.matrix(dens[,which(clcol==i)]), 1, mean));
              dens <- den; clcol <- indx; ind.clust <- TRUE
            } else indx <- which(clcol==input$sts.MNCSDE)
          }
          ts.plot(dens[,indx], col=clcol[indx], main=paste(input$results,"-",j,"-",m.lab), xlab="Frequency", ylab="", ylim=range.dens, gpars=list(xaxt="n"), log=ifelse(input$log.psd,"y",""));
          axis(1,summary(1:nrow(dens))[-4],summary(res$init$omegas)[-4])#c(0,.125,.25,.375,.5))
          if (!is.null(imodelPar())) if (input$sh.tr) {
            if (ind.clust) ind <- indx else ind <- unique(imodelPar()$rnd[indx])
            for (i in ind) lines(res$init$truspec[j,,i],type="l",lty=2,col="grey")
          }
          axis(1,summary(1:nrow(dens)),summary(res$init$omegas))#c(0,.125,.25,.375,.5))
        }
      }
    }
  })
}



#' laucnh a shiny app
#'
#' @param ... Any possible variables
#'
#' @return open a shiny app
#' @export
launch_app <- function(...) {
  shiny_env <- new.env()
  environment(shiny_ui) <- shiny_env
  environment(shiny_server) <- shiny_env
  shiny::shinyApp(
    ui = shiny_ui,
    server = shiny_server
  )
}
