##########################################
# Bayesian random effects consensus calculations 
# for inter-laboratory comparisons 
# Author: Juris Meija, NRC Canada
# Version 6, Nov 2022
##########################################

# Brief description
# This calculator fits random effects model to inter-laboratory data

require(shiny)
require(shinyjs) # Javascript for shiny
require(shinyFeedback)
require(rhandsontable) # Hot tables in shiny
require(rjags)
require(logspline) # smooth density kernel
require(DT)
require(dplyr)
require(nortest)  # Normality test
require(symmetry) # Symmetry test
require(runjags)
require(knitr) # Generate reports
require(rmarkdown)
# student t distribution
dstudent_t = function (x, df, mu = 0, sigma = 1) {
  if (isTRUE(any(sigma <= 0))) { stop2("sigma must be greater than 0")}
    dt((x - mu)/sigma, df = df)/sigma
}

# Initial dataset (CCQM K148a, purity of bisphenol A)
dd='BAM	TRUE 995.4	0.46
BIPM	FALSE 993.3	1.7
CENAM	FALSE 977.02	0.26
EXHM	TRUE 994.23	0.64
GLHK	TRUE 996.3	2.5
HSA	TRUE 995.2	1.5
INMETRO	TRUE 995.7	0.6
KRISS	TRUE 995.87	0.82
LGC	TRUE 995.8	1.2
NIM	TRUE 996.41	1.08
NIMT	FALSE 987.8	2.76
NMIA TRUE	997	0.9
NMIJ	TRUE 996.1	0.5
NMISA	FALSE 989.6	4
NRC	TRUE 993.7	2.4
UME	TRUE 996.64	3.03
VNIIM	FALSE 997.75	0.146'
dd=matrix(strsplit(dd,'[\n\t ]')[[1]],4)

## INITIAL DATASET
init.df = data.frame(include=as.logical(dd[2,]), lab=dd[1,], result=as.double(dd[3,]), uncertainty=as.double(dd[4,]) )

### server file
server <- function(input, output, session) {
  
  v <- reactiveValues(ready=FALSE, mu=NULL, tau=NULL, x=NULL, DoE=NULL, 
                      summary=NULL, failedconvergence=FALSE, data = NULL, digits = 2,
                      chi2 = 1, chi2upper = 1, ad = 1, sym = 1)
  
  output$hot <- renderRHandsontable({
    
    if(is.null(input$hot)) { DF = init.df } else { DF = hot_to_r(input$hot) }  
    
    names(DF) <- c('Include?','Laboratory','Result','Uncertainty')
    myindex = which(DF[,1]=='FALSE')-1
    rhandsontable(DF, readOnly = FALSE, stretchH = "all", selectCallback = TRUE, myindex = myindex) %>%
      hot_context_menu(allowColEdit = FALSE ) %>%
      hot_validate_numeric(cols=3:4) %>%
      hot_col(2, format = '', halign = 'htCenter', valign = 'htTop') %>%
      hot_col(1, halign = 'htCenter') %>%
      hot_col(c(3,4), format = "0.000000", halign = 'htCenter') %>%
      #hot_col(c(3,4), format = paste0("0.", paste0(rep(0, 5), collapse='')), halign = 'htCenter') %>%
      hot_col(c(2,3,4), renderer = "function(instance, td, row, col, prop, value, cellProperties) {
            Handsontable.renderers.TextRenderer.apply(this, arguments);
            if (instance.params) {
            hrows = instance.params.myindex
            hrows = hrows instanceof Array ? hrows : [hrows] 
          }
          if (instance.params && hrows.includes(row)) {td.style.background = 'lightpink';}
          else {td.style.background = 'darkseagreen';  }
        }
        ") %>%
    hot_col(1, type = 'checkbox')
      
  })
  
  observeEvent(input$hot, {
    z=hot_to_r(input$hot)
    incl = as.logical(z[,1])
    updateNumericInput(session, "mu_prior", value = round(median(as.double(z[incl,3]), na.rm=T), 1+v$digits+input$digits) )
    updateNumericInput(session, "u_prior", value = round(abs(median(as.double(z[incl,3]), na.rm=T)), 1+v$digits+input$digits) )
    updateNumericInput(session, "tau_prior", value = round(mad(as.double(z[incl,3]), na.rm=T), 1+v$digits+input$digits) )
  }, ignoreInit = TRUE)

  observeEvent(input$Niter,  { if( input$Niter >= 50000 ) { hideFeedback("Niter") } else { showFeedbackDanger("Niter", "Invalid number of iterations") }  })
  observeEvent(input$mu_prior,  { if( is.numeric(input$mu_prior) ) { hideFeedback("mu_prior") } else { showFeedbackDanger("mu_prior", "Invalid value") }  })
  observeEvent(input$u_prior,  { if( is.numeric(input$u_prior)&input$u_prior>0 ) { hideFeedback("u_prior") } else { showFeedbackDanger("u_prior", "Invalid value") }  })
  observeEvent(input$tau_prior,  { if( is.numeric(input$tau_prior)&input$tau_prior>0 ) { hideFeedback("tau_prior") } else { showFeedbackDanger("tau_prior", "Invalid value") }  })
  
  observeEvent(input$fitmodel, {
      output$adtest <- renderUI({  NULL })
      if (is.null(input$hot)) { z = init.df } else { z = hot_to_r(input$hot) }
      
      v$data <- list(lab=z[,2], include=as.logical(z[,1]), result=as.numeric(z[,3]), uncertainty=as.numeric(z[,4]),
           muprior = input$mu_prior,
           uprior = input$u_prior,
           tauprior = input$tau_prior
      )
    
    # validate input data
    if(sum(unlist(do.call(rbind,input$hot$params$data)[,1]))<2) {
      output$comment <- renderUI(HTML(paste("<p style='color:red'><b>WARNING: At least two measured values are required!</b></p>")))
      show("comment")
      hide("f")
      hide("table")
      hide("consensus")
      hide("modeltext")
      hide("tabletitle")
      v$ready <- FALSE
    }
    
    d=do.call(rbind, input$hot$params$data)
    if( any(sapply(c(d[,3:4]), is.null))  ) {
        output$comment <- renderUI(HTML(paste("<p style='color:red'><b>WARNING: All input results must be numeric!</b></p>")))
        show("comment")
        hide("f")
        hide("table")
        hide("consensus")
        hide("modeltext")
        hide("tabletitle")
        v$ready <- FALSE
      }  
    
    validate(
      need(input$Niter >= 50000, 'NITER'),
      need(is.numeric(input$mu_prior), 'MU_PRIOR'),
      need(is.numeric(input$u_prior)&input$u_prior>0, 'U_PRIOR'),
      need(is.numeric(input$tau_prior)&input$tau_prior>0, 'TAU_PRIOR'),
      need(!any(sapply(c(d[,3:4]), is.null)), 'RESULTS'),
      need(length(unlist(do.call(rbind,input$hot$params$data)[,3]))==length(unlist(do.call(rbind,input$hot$params$data)[,4])), 'RESULTS'),
      need(sum(unlist(do.call(rbind,input$hot$params$data)[,1]))>1, 'RESULTS')
    )
    hide("comment")
    show("f")
    show("table")
    show("consensus")
    show("tabletitle")
    show("modeltext")
    
    # TEST 1 # Overdispersion test of data (chi2)
    v$chi2 <- sum((v$data$result[v$data$include] - median(v$data$result[v$data$include]))^2/v$data$uncertainty[v$data$include]^2)
    v$chi2upper <- qchisq(0.05, df=length(v$data$result[v$data$include]) - 1, lower.tail = FALSE)
    
    # TEST 2 # Normality test of data (Anderson-Darling test on standardized data)
    if(length(v$data$result[v$data$include]) > 7) v$ad <- ad.test((v$data$result[v$data$include] - median(v$data$result[v$data$include]))/v$data$uncertainty[v$data$include])$p.value
    
    # TEST 3 # Symmetry test of data around their median (Miao, Gel, Gaswirth 2006)
    v$sym <- symmetry_test(x = v$data$result[v$data$include], stat = 'MGG')$p.value
    
    # setup the statistical model
    # random effects (inference)
    m1.n  = paste0("for(i in incl) { theta[i] ~ dnorm(mu, 1.0/ut^2) }")
    m1.t4 = paste0("for(i in incl) { theta[i] ~ dt(mu, 1.0/ut^2, 4) }")
    m1.sn = paste0("for(i in incl) { z[i] ~ dlnorm(0, 1.0/ut^2) }",
                   "for(i in incl) { theta[i] <- mu + z[i] - zmean}", 
                   "zmean <- exp(0.5*ut^2)", collapse="\n")
    m1.l  = paste0("for(i in incl) { theta[i] ~ ddexp(mu, 1.0/(ut/sqrt(2))) }")
    
    # random effects (predictive)
    m2.n  = paste0("for(i in 1:N) { e_r[i] ~ dnorm(0, 1.0/ut^2) }")
    m2.t4 = paste0("for(i in 1:N) { e_r[i] ~ dt(0, 1.0/ut^2, 4) }")
    m2.sn = paste0("for(i in 1:N) { z_r[i] ~ dlnorm(0, 1.0/ut^2) }",
                   "for(i in 1:N) { e_r[i] <- z_r[i] - zmean }", collapse="\n")
    m2.l  = paste0("for(i in 1:N) { e_r[i] ~ ddexp(0, 1.0/(ut/sqrt(2))) }")
    
    modelstring=paste("model {",
      c(m1.n, m1.t4, m1.sn, m1.l)[as.double(input$reffect)],
      paste0("for(i in incl) { x[i] ~ dnorm(theta[i], 1.0/u[i]^2) }"),
      paste0("mu ~ dnorm(", input$mu_prior, ",1.0/(",input$u_prior,")^2 )"),
      paste0("ut ~ dt(0, 1.0/(", input$tau_prior, ")^2, 1)T(0,)"),
      paste0("for(i in 1:N) { e_f[i] ~ dnorm(0, 1.0/u[i]^2) }"),
      c(m2.n, m2.t4, m2.sn, m2.l)[as.double(input$reffect)],
      paste0("for(i in 1:N) { DoE[i] <- x[i] - mu + e_f[i] + e_r[i] }"),
      paste0("}"),
      collapse="\n")
    
    jags.inits = list(mu = median(v$data$result[which(v$data$include)]), ut = input$tau_prior)
    out = list()
    model=jags.model(textConnection(modelstring), 
                     data=list(x=v$data$result, u=v$data$uncertainty, N=length(v$data$result), incl = which(v$data$include)
                               ), 
                     inits = jags.inits, quiet=TRUE, n.chains = 4 )
    withProgress(message="Fitting the model...", {
      for(i in 1:5) {
        update(model, n.iter=input$Niter%/%5, progress.bar='none')
        incProgress(1/10, detail=paste0(10*i,'%'))
        }
      for(i in 1:5) {
        out[[i]] <- combine.mcmc(coda.samples(model=model, variable.names=c("mu","ut","DoE","theta"), n.iter=input$Niter%/%5, thin=20, progress.bar='none'))
        incProgress(1/10, detail=paste0(50+10*i,'%'))}
    })
    
    dic <- dic.samples(model, n.iter=1000)
    
    output <- combine.mcmc(out)
    out <- as.mcmc(output)
       
    # MCMC draws from the posterior pdf of measurement uncertainty
    v$mu <- c(out[,'mu'])
    v$digits <- 2+abs(log10(sd(v$mu)))%/%1
    
    # MCMC draws from the posterior pdf of the dark uncertainty
    v$tau <- c(out[,'ut'])
    # MCMC draws from the posterior pdf of degrees-of-equivalence
    v$DoE <- out[,grepl('DoE',colnames(out))]
    v$x <- v$data$result
    # summary output
    v$summary <- summary(output)
    # MCMC convergence test
    v$failedconvergence <- min(p.adjust(2*pnorm(-abs(geweke.diag(out)[[1]])), method="BH")) < 0.05
    v$ready <- TRUE
    
  })
  
  main_plot <- reactive( { 
    
    if( v$ready & length(v$mu)>1 & length(v$tau)>1){  
      par(mar=c(10,2,3,3))
      layout(matrix(c(1,2,2,3,3), nrow = 1, ncol = 5, byrow = TRUE))  
      
      qq <- quantile(v$mu, c(0.01, 0.025, 0.975, 0.99, 0.50))
      
      # smoothed kernel density using splines to approximate the log-density of mu
      m <- oldlogspline(v$mu)
      x0.logspline = seq(qq[1],qq[4],length.out = 200)
      x1.logspline = seq(qq[2],qq[3],length.out = 200)
      y0.logspline = sapply(x0.logspline, function(x) doldlogspline(x, m))
      y1.logspline = sapply(x1.logspline, function(x) doldlogspline(x, m))
      
      plot(m, main='Consensus value', pch='', ylim=c(0, 1.05*max(y0.logspline)),
           col='steelblue3', yaxt='n', ylab='', xlab='', xlim=c(qq[1], qq[4]), cex.axis=1.5, cex.main=2)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey98")
      # 5-95 % credible interval for the posterior of mu
      polygon(c(x1.logspline, rev(x1.logspline)), c(y1.logspline, rep(0,length(x1.logspline))), col='#4F94CD7F', border='gray98')
      # draw prior distribution for the consensus value
      if(input$plotprior){
        xx <- seq(qq[1], qq[4], length.out = 200)
        lines(x=xx, y=dnorm(xx, input$mu_prior, input$u_prior), col='tomato', lwd=4)  
      }
      
      # posterior distribution for the consensus value
      lines(x=x0.logspline,y=y0.logspline, lwd=4, col='steelblue3')    
      segments(qq[2],0,qq[3],0)
      points(x=qq[5],y=0, pch=19, cex=2)
      if(input$plotprior){
        legend("topright", legend = c('prior','posterior'), col = c('tomato','steelblue3'), lwd = c(4, 4), text.col = c('tomato','steelblue3'), text.font = c(2, 2), border = 'gray40')
      }
      
      graphics::box()
      
      # Individual laboratory plot
      incl = v$data$include
      u = v$data$uncertainty
      N=length(v$data$result)
      
      cc=rep('tomato',N)
      cc[incl] <- 'gray5'
      
      o = order(v$data$result)
      ys = v$data$result[o]
      # get robust median of the u_tau
      Tys = sqrt(u[o]^2 + mean(v$tau^2))
      Uys = u[o]
      Umu = abs(diff(quantile(v$mu, c(0.025, 0.975))))/2
      
      df=data.frame(lab=v$data$lab[o], value=ys,
                    Ulow=ys-2*Uys,Uhigh=ys+2*Uys,Tlow=ys-2*Tys,Thigh=ys+2*Tys,
                    order=1:N, color=cc[o])
      
      plot(x=1:N,y=ys, main='Results', pch=19, col='steelblue3', 
           ylab='', xaxt='n',xlab='', ylim=range(df$Tlow,df$Thigh),
           cex.axis=1.5, cex.main=2)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey98")
      segments(1:N,par("usr")[3],1:N,mean(v$mu)-1.01*Umu,lty=3,lwd=0.7,col='gray50')
      segments(1:N,mean(v$mu)+1.01*Umu,1:N,par("usr")[4],lty=3,lwd=0.7,col='gray50')
      rect(0, mean(v$mu)-Umu, 2*N, mean(v$mu)+Umu, col = "#4F94CD7F", border='#4F94CD7F')
      abline(h=mean(v$mu), col='white')
      
      segments(1:N,df$Tlow,1:N,df$Thigh,lwd=2,col=df$color)
      segments(1:N,df$Ulow,1:N,df$Uhigh,lwd=5, col=df$color)
      points(x=1:N,y=df$value, pch=21,bg='white',cex=2)
      
      d=subset(df,color=='gray5')
      axis(side=1,at=d$o,labels=d$lab, font=2, las=2, col.axis='gray5', cex.axis=1.5)
      d=subset(df,color=='tomato')
      axis(side=1,at=d$o,labels=d$lab, font=2, las=2, col.axis='tomato', cex.axis=1.5)
      graphics::box()
      
      # Degrees of equivalence
      tt = t(apply(v$DoE,2, function(x) c(mean(x), abs(quantile(x, 0.025) - quantile(x, 0.975))/2 )))
      ys=tt[o,1]
      Uys=tt[o,2]
      plot(x=1:N,y=ys, main='Degrees of equivalence', pch=19, col='steelblue3', 
           ylab='', xaxt='n',xlab='', ylim=range(ys-Uys,ys+Uys),
           cex.axis=1.5, cex.main=2)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey98")
      segments(1:N,par("usr")[3],1:N,par("usr")[4],lty=3,lwd=0.7,col='gray50')
      abline(h=0, col='tomato',lwd=2)
      segments(1:N,ys-Uys,1:N,ys+Uys,lwd=2)
      points(x=1:N,y=ys, pch=21,bg='white',cex=2)
      axis(side=1,at=1:N,labels=v$data$lab[o],las=2, font=2, col.lab='black', cex.axis=1.5)
      graphics::box()
    }
    
    } )
  
  output$f <- renderPlot({ main_plot() })
  
  output$table <- renderDT( if(v$ready & !is.null(v$mu) & !is.null(v$tau)){  
      
      tt = data.frame( lab = v$data$lab,
                       DoE = v$x - mean(v$mu),
                         u = apply(v$DoE,2,sd),
                       U95 = apply(v$DoE,2,function(x) abs(diff(quantile(x,c(0.025,0.975))))/2  ),
                       row.names = NULL
                       )
      colnames(tt) <- c('Laboratory','DoE','Uncertainty', 'U95')
      datatable(tt[order(tt$DoE),], rownames= FALSE) %>% formatRound(2:4, digits=v$digits+input$digits) 
      
          } else ( return(NULL) ) )
  
   modeltext <- eventReactive(input$fitmodel, {
     updateNumericInput(session, "mu_prior", value = round(v$data$muprior, v$digits+input$digits) )
     updateNumericInput(session, "u_prior", value = round(v$data$uprior, v$digits+input$digits) )
     updateNumericInput(session, "tau_prior", value = round(v$data$tauprior, v$digits+input$digits) )
     
    if(v$ready) HTML(paste('<b>Statistical model</b>',
               ifelse(input$reffect==1,'Hierarchical random effects model (normal-normal):',
               ifelse(input$reffect==2,'Hierarchical random effects model (t4-normal):',
               ifelse(input$reffect==3,'Hierarchical random effects model (skew normal-normal):',
                                       'Hierarchical random effects model (laplace-normal):'))),
               paste0('&nbsp;&nbsp;&nbsp; result[i] ~ normal(mean = <i>m</i> + <i>t</i>, sd=unc[i]) for each laboratory part of the consensus building <i>i</i> (<i>i</i> = 1...', length(which(v$data$include)),")"),
               ifelse(input$reffect==1, paste0('&nbsp;&nbsp;&nbsp; <i>t</i> ~ normal(mean = 0, sd=<i>u</i><sub>t</sub>)'),
               ifelse(input$reffect==2, paste0('&nbsp;&nbsp;&nbsp; <i>t</i> ~ t4(mean = 0, scale=<i>u</i><sub>t</sub>)'),
               ifelse(input$reffect==3, paste0(''),
                             paste0('&nbsp;&nbsp;&nbsp; <i>t</i> ~ laplace(mean = 0, sd=<i>u</i><sub>t</sub>)')
                             ))),
               'The following priors were used for the model parameters:',
               paste0("&nbsp;&nbsp;&nbsp; <i>m</i> ~ normal(mean = ", input$mu_prior, ", sd=",input$u_prior,")"),
               paste0("&nbsp;&nbsp;&nbsp; <i>u</i><sub>t</sub> ~ halfCauchy(median =",input$tau_prior,")"),
               'Degree-of-equivalence (DoE) is the difference between the laboratory result and the consensus value:',
               paste0('&nbsp;&nbsp;&nbsp; DoE[j] = result[j] - <i>m</i> + <i>e</i><sub>1</sub>[j] + <i>e</i><sub>2</sub> for all laboratories <i>j</i> (<i>j</i> = 1...', length(v$data$result),")"),
               paste0('&nbsp;&nbsp;&nbsp; e<sub>1</sub>[j] ~ normal(mean = 0, sd=unc[j])'),
               ifelse(input$reffect==1, paste0('&nbsp;&nbsp;&nbsp; e<sub>2</sub> ~ normal(mean = 0, sd=<i>u</i><sub>t</sub>)'),
               ifelse(input$reffect==2, paste0('&nbsp;&nbsp;&nbsp; e<sub>2</sub> ~ t4(mean = 0, scale=<i>u</i><sub>t</sub>)'),
               ifelse(input$reffect==3, paste0(''),
                              paste0('&nbsp;&nbsp;&nbsp; e<sub>2</sub> ~ laplace(mean = 0, sd=<i>u</i><sub>t</sub>)')
                              ))),
               paste0(''),
               paste0("The model was fit to data using Markov-chain Monte Carlo sampling in R using rjags"),
               paste0("All uncertainty bars in the above plots correspond to 95 % confidence."),
               paste0(""),
               paste0("NRC Interlaboratory Consensus Calculator (2020-2022) Juris Meija, v6 November 2022"),
               sep="<br/>"))
    })
   tabletitle <- eventReactive(input$fitmodel, { if(v$ready) HTML(paste('<b>Unilateral Degrees-of-equivalence</b>')) })
   consensustext <- eventReactive(input$fitmodel, { if(v$ready) HTML(paste0("<p style='color:#36648b'><b>Consensus value: ", formatC(mean(v$mu), digits=v$digits+input$digits, format='f'),' (<i>u</i> = ',formatC(sd(v$mu),digits=v$digits+input$digits,format='f'),') with expanded uncertainty <i>U</i><sub>95 &#37;</sub> = ', formatC(abs(quantile(v$mu, 0.025) - quantile(v$mu, 0.975))/2, digits=v$digits+input$digits, format='f'),' from ',length(v$data$result[v$data$include]),' laboratory results</b></p>')) })
   darkunctext <- eventReactive(input$fitmodel, { if(v$ready) HTML(paste0("<p style='color:#36648b'><b>Dark uncertainty: <i>u</i><sub>t</sub> = ", formatC(mean(v$tau), digits=v$digits+input$digits, format='f'),'</b> (<i>u</i><sup>2</sup><sub>t</sub> = ', formatC(mean(v$tau^2), digits=v$digits+input$digits, format='f')  ,')</p>')) })
   
   # MCMC convergence test
   convtext <- eventReactive(input$fitmodel, { 
     if (v$ready & v$failedconvergence) HTML( paste("<p style='color:red'>WARNING: MCMC may not have reached equilibrium. Results are not reliable. Try increasing the number of iterations.</p>" ) )
   })
 
   # The results of the three tests on data (overdispersion, normality, symmetry)
   test1.text <- eventReactive(input$fitmodel, { 
     chi2res = c(v$chi2<=(length(v$data$include)-1), v$chi2>(length(v$data$include)-1)& v$chi2<=v$chi2upper, v$chi2>v$chi2upper)
     if (v$ready) HTML(paste('Overdispersion of the random effects: ', c('No overdispersion present', 'Borderline consistency', 'Clear evidence of overdispersion')[chi2res], collapse=''))
   })
   
   test2.text <- eventReactive(input$fitmodel, { 
     if (v$ready) HTML(paste('Normality test of the random effects: ', paste(ifelse(v$ad<0.05, 'NOT NORMAL, try Laplace model', 'NORMAL'), '( p-value:', formatC(v$ad, digits=2),')') ))
   })
   
   test3.text <- eventReactive(input$fitmodel, { 
     if (v$ready) HTML(paste('Symmetry test of the random effects: ', paste(ifelse(v$sym<0.05, 'NOT SYMMETRIC, try skew normal model', 'SYMMETRIC'), '( p-value:', formatC(v$sym, digits=2),')') ))
  })
   
   output$consensus <- renderUI({  consensustext() })
   output$darkuncertainty <- renderUI({  darkunctext() })
   output$modeltext <- renderUI({  modeltext() })
   output$tabletitle <- renderUI({  tabletitle() })
   
   output$test1 <- renderUI({ test1.text()  })
   output$test2 <- renderUI({ test2.text()  })
   output$test3 <- renderUI({ test3.text()  })
   
   output$convtext <- renderUI({  convtext() })
   
   # not yet implemented
   output$downloadReport <- downloadHandler(filename = "report.html",
                                            content = function(file){
                                              # generate PDF
                                              knit2html("reports/report.Rmd")
                                              
                                              # copy pdf to 'file'
                                              file.copy("report.html", file)
                                              
                                              # delete generated files
                                              file.remove("report.html")
                                              
                                              # delete folder with plots
                                              unlink("figure", recursive = TRUE)
                                            },
                                            contentType = "text/html"
   )
     
   
   shinyjs::onclick("toggleextra", shinyjs::toggle(id = "filterextra", anim = TRUE))
  
  }

### user interface
ui <- fluidPage(
  useShinyFeedback(),
  titlePanel( title="Interlaboratory Consensus Calculator (v6)" ),
  shinyjs::useShinyjs(),
  sidebarLayout(
    
    sidebarPanel(
      h5(tags$b("Enter (paste) the observed results")),
      rHandsontableOutput("hot"),
      helpText("right-click to add or delete rows"),
      h5(tags$div(HTML(' <i class="fas fa-exclamation" style = "color:red"></i> Excluded laboratories are not used to obtain the consensus value'))),
      textOutput("error"),
      h5("MODEL PARAMETERS ", a(id = "toggleextra", "show/hide")),
      shinyjs::hidden(div(id = "filterextra",
                          fluidRow(
                            column(6, numericInput("Niter", label = "MCMC draws", value = 100000, min = 50000, step = 10000, width =  '85%')),
                            column(6, sliderInput("digits", label = "Additional decimal digits", value = 0, max = 4, min = 0, step = 1, ticks = FALSE, width='85%')),
                            column(12, h5('Parameters for prior distributions of the fixed effects')),
                            column(6, numericInput("mu_prior", label = "location prior (N)", value = 1, width = '85%')),
                            column(6, numericInput("u_prior", label = "scale prior (hC)", value = 1, min = 0, width = '85%')),
                            column(6, checkboxInput("plotprior","Plot prior distribution of m?",value = FALSE), offset=0.5),
                            column(12, h5('Parameters for prior distributions of the random effects')),
                            column(12, selectInput("reffect","Distribution",choices = c('Normal'=1,'Student t4'=2,'Laplace'=4), selected=1, width='40%')),
                            column(6, numericInput("tau_prior", label = "scale prior (hC)", value = 1, min = 0, width = '85%')),
                            conditionalPanel(condition = "input.reffect == 3", 
                                             column(6, numericInput("skew_prior", label = "skew prior (N)", value = 0, width = '85%'))
                                             )
                            ))),
      br(),
      conditionalPanel(condition = "!$('html').hasClass('shiny-busy')", actionButton(inputId = "fitmodel", label = "Fit model", icon = icon('chart-line'))),
      conditionalPanel(condition = "$('html').hasClass('shiny-busy')",  actionButton(inputId = "fitmodel", label = "busy...", icon = icon('hourglass-half'))), 
      br()
    ),
    
    mainPanel(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      
                 htmlOutput("comment"),br(),
                 plotOutput("f", width="95%"),br(),
                 htmlOutput("consensus"),htmlOutput("darkuncertainty"),
                 # currently in development
				 # downloadButton("downloadReport", "Generate PDF report"),
                 htmlOutput("convtext"),htmlOutput("test1"),htmlOutput("test2"),htmlOutput("test3"),br(),
                 htmlOutput("tabletitle"),br(),
                 DTOutput("table", width="95%"),br(),
                 htmlOutput("modeltext"),br()
                 
      # fluidRow(column(helpText("This calculator fits hierarchical random effects model to the data using Bayesian method. The random effects can be specified to follow normal distribution or Laplace distribution."), width=11)),
      
    )
  )
)

shinyApp(ui = ui, server = server)