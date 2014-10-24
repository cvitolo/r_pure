loglogplot <- function(d.f){

  ##  Open a new default device.

  #get( getOption( "device" ) )()

  ##  Plot the data, hiding the points for now to prevent the calls to
  ##  abline() from drawing over the points.

  plot(
    y ~ x,
    data = d.f,
    type = "n",
    log  = "xy",
    main = "Log-log Plot",
    xlab = "Return period",
    ylab = "[mm/d]"#,
    #xlim = c( 1, 1000 ),
    #ylim = c( 1, 1000 )
  )

  ##  Put grid lines on the plot, using a light blue color ("lightsteelblue2").

  abline(
    h   = c( seq( 1, 9, 1 ), seq( 10, 90, 10 ), seq( 100, 1000, 100 ) ),
    lty = 3,
    col = colors()[ 440 ] )

  abline(
    v   = c( seq( 1, 9, 1 ), seq( 10, 90, 10 ), seq( 100, 1000, 100 ) ),
    lty = 3,
    col = colors()[ 440 ] )

  ##  Draw the points over the grid lines.

  points( y ~ x, data = d.f )
  lines( y ~ x, data = d.f )

  ##  Redraw the plot box over the grid lines.

  box()

}
