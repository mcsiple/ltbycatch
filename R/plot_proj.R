#' Plot projections
#' @description plots outputs from several projections that result from a Projections() call.
#'
#' @param high A list containing output from Projections() (including a matrix of simulation trajectories) - corresponding to a high bycatch level (this is at the high end of the range defined by the user in the Shiny app)
#' @param med A list containing output from Projections() (including a matrix of simulation trajectories) - corresponding to a medium bycatch level (this is at the median of the low and high end of the range defined by the user)
#' @param low A list containing output from Projections() (including a matrix of simulation trajectories) - corresponding to a high bycatch level (this is at the high end of the range defined by the user)
#' @param years.plot - how many years to plot
#' @param spaghetti - either FALSE or a number, where the number is how many simulations to show from the same scenario
#' @param years.to.plot - how many years to plot on the x axis
#' @param ylims
#' @param K1plus
#' @param InitDepl
#' @param color.palette
#' @param multiplot
#'
#' @return A plot of 50 percent and 90 percent confidence intervals of population projections (if \code{spaghetti ==FALSE}) or a spaghetti plot (if \code{spaghetti ==TRUE}),  from \code{Projections()}.
#'
#' @export
plot_proj <- function(high,
                      med,
                      low,
                      years.plot = 50,
                      ylims,
                      spaghetti = FALSE,
                      years.to.plot=50,
                      K1plus=9000,
                      InitDepl=0.8,
                      color.palette = c("forestgreen","darkorange","red"),
                      multiplot=FALSE){
  high.col <-color.palette[3]
  med.col <-color.palette[2]
  low.col <- color.palette[1]

  #print(K1plus)
  #print(InitDepl)

  if(spaghetti){
    ts.length <- years.plot
    spag.df <- data.frame(year = as.numeric(rep(1:(ts.length+1),times=spaghetti)),
                          sim = as.factor(rep(1:spaghetti,each=ts.length+1)),
                          high = as.vector(t(high$trajectories[1:spaghetti,1:(ts.length+1)]))/K1plus,
                          med = as.vector(t(med$trajectories[1:spaghetti,1:(ts.length+1)]))/K1plus,
                          low = as.vector(t(low$trajectories[1:spaghetti,1:(ts.length+1)]))/K1plus
    )

    sdf <- melt(spag.df,id.vars=c('year','sim'))
    sdf$sim <- as.factor(sdf$sim)
    dlab1 <- paste("N[0] == ",round(InitDepl,digits=2),"*K")

    sdf <- sdf %>% mutate(variable = recode(variable,
                                            high = 'High end of \nbycatch range',
                                            med = 'Midpoint of \nbycatch range',
                                            low = 'Low end of \nbycatch range'))

    p <- sdf %>%
      ggplot(aes(x=year,y=value,group=sim,colour=variable)) +
      xlim(0,ts.length)+
      geom_path(lwd=1,alpha=0.5) +
      scale_colour_manual("Bycatch level",values = rev(color.palette)) +
      theme_classic(base_size=18) +
      xlab('Year') + ylab('Abundance (N/K)') +ylim(0,1)+
      annotate("label", x=ts.length*.8, y=.95, label = dlab1,parse=TRUE,size=6) +
      labs(title = 'Model projections',
           subtitle = 'Median population size and quantiles',
           caption = '95%: lightly shaded; 50%: heavily shaded') +
      theme(legend.position = 'bottom',
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18),
            axis.title.y = element_text(size=20))


    # if(multiplot){
    #   p <- p + theme(legend.position = "none")
    # }
    p

  }else{
    probs <- c(0.05,0.25,0.5,0.75,0.95)

    # low med and high refer to low med and high bycatch
    summary.high <- apply(high$trajectories[,1:years.to.plot],MARGIN = 2,FUN = quantile,probs=probs,na.rm=T)
    summary.med <- apply(med$trajectories[,1:years.to.plot],MARGIN = 2,FUN = quantile,probs=probs,na.rm=T)
    summary.low <- apply(low$trajectories[,1:years.to.plot],MARGIN = 2,FUN = quantile,probs=probs,na.rm=T)
    ts.length <- ncol(summary.high)


    t.high <- data.frame(t(summary.high),blvl="high",year=1:years.to.plot)
    t.med <- data.frame(t(summary.med),blvl="med",year=1:years.to.plot)
    t.low <- data.frame(t(summary.low),blvl="low",year=1:years.to.plot)
    all <- rbind(t.high,t.med,t.low)

    colnames(all) <- c("lo90","lo50","median","hi50","hi90","blvl","year")

    all <- all %>%
      mutate(blvl = factor(blvl)) %>%
      mutate(blvl = recode_factor(blvl,
                                  high = 'High end of \nbycatch range',
                                  med = 'Midpoint of \nbycatch range',
                                  low = 'Low end of \nbycatch range'))
    #print(levels(all$blvl))
    dlab <- paste("N[0] == ",round(InitDepl,digits=2),"*K")

    p <- all %>%
      ggplot(aes(x=year)) +
      geom_ribbon(aes(ymin=lo90/K1plus,ymax=hi90/K1plus,fill=blvl),alpha=0.5) +
      geom_ribbon(aes(ymin=lo50/K1plus,ymax=hi50/K1plus,fill=blvl),alpha=0.5) +
      geom_line(aes(y=median/K1plus,colour=blvl),lwd=1.1) +
      ylim(0,1) +
      scale_fill_manual("Bycatch level",values = rev(color.palette)) +
      scale_colour_manual("Bycatch level",values = rev(color.palette)) +
      xlab("Year") +
      ylab(paste("Population size relative to K")) +
      labs(title = 'Model projections',
           subtitle = 'Median population size and quantiles',
           caption = '95%: lightly shaded; 50%: heavily shaded') +
      theme_classic(base_size=18) +
      annotate("label", x=years.to.plot*.8, y=.95, label = dlab,parse=TRUE,size=6) +
      theme(legend.position = 'bottom',
            axis.text.x = element_text(size=18),
            axis.text.y = element_text(size=18),
            axis.title.y = element_text(size=20),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size=16))

    # if(multiplot){
    #   p <- p + theme(legend.position = "none")
    # }
    p

  }
}
