    # Gnuplot script file for plotting data in file "data(i).dat"
    # This file is called   plotall.p    
    
    set term pdfcairo
    set output 'compareNew.pdf'

    # set terminal postscript eps enhanced color 
    # set output '|ps2pdf - outputfile.pdf'

    set   autoscale                        # scale axes automatically
    unset logscale                          # remove any log-scaling
    set logscale y                         
    unset label                            # remove any previous labels
    set xtic auto                          # set xtics automatically
    set ytic auto                          # set ytics automatically
    set title "LU Decomposition Duration Comparison for different n"
    set title offset 0,-0.5
    set xlabel "m"
    set xlabel offset 0,1
    set ylabel "time / s"
    set ylabel offset 2.1
    #set key at 0.01,100
    set key font ",10"
    set key bottom right
    set label "n = 2000" at 20,6
    set label "n = 3000" at 30,18
    set label "n = 4000" at 40,37
    set label "n = 5000" at 50,70
    #set arrow from 0.0028,250 to 0.003,280
    set xr [1:700]
    set yr [2.7:80]
    set ytics 1,2,80
    plot "data5.dat" using 1:2 title 'Block, n = 5000' w l, \
    "data5.dat" using 1:3 title 'BLAS, n = 5000' w l,\
    "data4.dat" using 1:2 title 'Block, n = 4000' w l , \
    "data4.dat" using 1:3 title 'BLAS, n = 4000' w l,\
    "data3.dat" using 1:2 title 'Block, n = 3000' w l , \
    "data3.dat" using 1:3 title 'BLAS, n = 3000' w l,\
    "data2.dat" using 1:2 title 'Block, n = 2000' w l , \
    "data2.dat" using 1:3 title 'BLAS, n = 2000' w l
    
    set out                                 # gives output as declared at the top (e.g. pdf)
#command in gnuplot prompt: load "plotall.p"