#!/bin/zsh

for f in dump_confinement_f_0.6.txt dump_confinement_f_1.0.txt dump_confinement_Df_0.6.txt dump_confinement_Df_1.0.txt dump_elastic_f_0.6.txt dump_elastic_f_1.0.txt dump_elastic_Df_0.6.txt dump_elastic_Df_1.0.txt dump_pressure_f_0.6.txt dump_pressure_f_1.0.txt dump_pressure_Df_0.6.txt dump_pressure_Df_1.0.txt dump_finite_differences_matrices_0.6.txt dump_finite_differences_matrices_1.0.txt dump_inner_coords_0.6.txt dump_inner_coords_1.0.txt; do;
    echo $f
    cat $f | sed -e 's/.*Sparse.*//' -e 's/.*Float64.*//' > $f.filtered
done;

#vimdiff dump_confinement_f_0.6.txt.filtered dump_confinement_f_1.0.txt.filtered
#vimdiff dump_confinement_Df_0.6.txt.filtered dump_confinement_Df_1.0.txt.filtered
#vimdiff dump_elastic_f_0.6.txt.filtered dump_elastic_f_1.0.txt.filtered
#vimdiff dump_elastic_Df_0.6.txt.filtered dump_elastic_Df_1.0.txt.filtered
#vimdiff dump_pressure_f_0.6.txt.filtered dump_pressure_f_1.0.txt.filtered
#vimdiff dump_pressure_Df_0.6.txt.filtered dump_pressure_Df_1.0.txt.filtered
vimdiff dump_inner_coords_0.6.txt.filtered dump_inner_coords_1.0.txt.filtered
#vimdiff dump_finite_differences_matrices_0.6.txt.filtered dump_finite_differences_matrices_1.0.txt.filtered
