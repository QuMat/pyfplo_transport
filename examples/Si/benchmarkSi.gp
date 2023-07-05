set term epslatex color size 8.5 cm, 12 cm
epsfile='benchmarkSi.eps'
set out epsfile

set multiplot
set size 1,0.5
set origin 0,0.5

set xlabel '$\varepsilon$ (eV)' offset 0,0.8
set ylabel '$\sigma$ (BoltzWann units)' offset 1.5,0

#here 3.67 comes from translating our units to Boltzann units, which involves the number np.pi * 7.7e-5/0.658

p [-5:5] "wien2k/gga_so_sp/wan2/wan2_tdf.dat"  w l lc rgb 'orange' ti 'Wien2k + Wannier90', "k_80_80_80/TDF_0_0.dat" u ($1+0.5):($2/3.67e-2) w l lc rgb 'red' ti 'FPLO'

set ylabel '$\sigma$ ($\Omega$m)$^{-1}$' offset 1.5,0

set origin 0,0
set xlabel '$\mu$ (eV)' offset 0,0.8
p [-1:1] "wien2k/gga_so_sp/wan2/wan2_elcond.dat" u 1:($2==300 ? $3 : 1/0)  w p lc rgb 'orange' ti 'Wien2k + Wannier90', "k_80_80_80/sigma_vs_mu_0_0.dat" u ($2+0.5):($3*1e8) w l lc rgb 'red' ti 'FPLO'
