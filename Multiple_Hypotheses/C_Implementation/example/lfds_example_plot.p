# Gnuplot script for plotting the lfds of the example problem
# Run 'lfds_example' to generate the file 'lfds_example.dat' first

# number od densities
N = 3

set xlabel "w"
set ylabel "PDF"

set terminal pdf enhanced color linewidth 2
set output "lfds_example.pdf"

plot for [col=2:N+1] 'lfds_example.dat' using 1:col with lines title 'q_'.(col-2)
