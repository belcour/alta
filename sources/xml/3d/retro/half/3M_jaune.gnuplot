set term epslatex standalone color font 8

set xlabel "incidence elevation (in degrees)"
set ylabel "BRDF"
set key on inside left

# output ABC fits
set output "yellow_retro_abc.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_abc_retro.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "retro ABC fit" w points pointtype 12, "./results/3d/retro/half/3M_jaune_abc_back.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "back ABC fit"

# output Beckmann fits
set output "yellow_retro_beck.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_beck_retro.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "retro Beckmann fit" w points pointtype 12, "./results/3d/retro/half/3M_jaune_beck_back.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "back Beckmann fit"

# output Blinn fits
set output "yellow_retro_blinn.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_blinn_retro.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "retro Blinn fit" w points pointtype 12, "./results/3d/retro/half/3M_jaune_blinn_back.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "back Blinn fit"

# output Lafotune fit
set output "yellow_retro_laf.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_laf.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Lafortune fit"


# output rational fits
set term wxt

set ylabel "BRDF x cosine"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR_dense__nbsgrid_162.alta" using (180/pi*$2):($2 > 0.0 && $3 > -0.01 && $3 < 0.01 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_rat_TK2TL.dat" using (180/pi*$2):($2 > 0.0 && $3 > -0.01 && $3 < 0.01) ? $4 : 1/0 title "rational interpolation"

set term epslatex standalone color font 8
set output "yellow_retro_rat.tex"
set yrange [0:GPVAL_Y_MAX]
replot
