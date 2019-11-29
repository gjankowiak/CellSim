import DelimitedFiles
import PyPlot

"""
plot_csv(filename, axis_labels, column; args...)

plots the data from CSV file filename
axis_labels: use the CSV headers as axis labels
column: use the corresponding column as plotted values (y-axis)
"""
function plot_csv(file, l=false, column=2; args...)
    println(file)

    d, h = DelimitedFiles.readdlm(file, ','; header=true)

    if size(d,1) == 1 && d[1] == "-"
        PyPlot.axhline(d[column]; args...)
    else
        sorted_idx = sortperm(d[:,1])
        println(d)
        PyPlot.plot(d[sorted_idx,1], d[sorted_idx,column]; args...)
    end

    if l
        PyPlot.xlabel(h[1])
        PyPlot.ylabel(h[column])
    end
end


PyPlot.figure()

# kb
PyPlot.subplot(241)
plot_csv("plot_data/before_recompute_nucleus/v_kb_nuc_new.csv", true, label="nuc new μn = 20")
plot_csv("plot_data/before_recompute_nucleus/v_kb_nuc_new_mun60.csv", true, label="nuc new µn=60")
plot_csv("plot_data/before_recompute_nucleus/v_kb_nuc_old.csv", true, label="nuc old μn = 20")
plot_csv("plot_data/before_recompute_nucleus/v_kb_nuc_old_mun60.csv", true, label="nuc old μn=60")
plot_csv("plot_data/before_recompute_nucleus/v_kb_nuc_old_mun100.csv", true, label="nuc old µn=100")

plot_csv("plot_data/before_recompute_nucleus/v_kb_nonuc_mun100.csv", false, label="nonuc μn=100", ls="dotted")
plot_csv("plot_data/before_recompute_nucleus/v_kb_nonuc_mun20.csv", false, label="nonuc μn=20", ls="dotted")
PyPlot.axhline(0.0, lw=0.5, color="black")
PyPlot.legend()
PyPlot.title("Velocity vs kb")

# mu
PyPlot.subplot(242)
plot_csv("plot_data/before_recompute_nucleus/v_A_mu_nuc_old.csv", true, label="nuc old")
plot_csv("plot_data/before_recompute_nucleus/v_A_mu_nonuc.csv", true, label="nonuc", ls="dotted")
PyPlot.axhline(0.0, lw=0.5, color="black")
PyPlot.legend()
PyPlot.title("Velocity vs µn")

PyPlot.subplot(243)
plot_csv("plot_data/before_recompute_nucleus/v_A_mu_nuc_old.csv", true, 3, label="nuc old")
PyPlot.axhline(0.0, lw=0.5, color="black")
PyPlot.legend()
PyPlot.title("Area vs µn")

# kb
PyPlot.subplot(244)
plot_csv("plot_data/before_recompute_nucleus/v_A_KMT_nuc_old.csv", true, label="nuc old")
PyPlot.axhline(0.0, lw=0.5, color="black")
PyPlot.legend()
PyPlot.title("Velocity vs kMT")

# channel depth
PyPlot.subplot(245)
plot_csv("plot_data/before_recompute_nucleus/v_beta_nuc_old.csv", true, label="nuc old")
plot_csv("plot_data/before_recompute_nucleus/v_beta_nuc_new.csv", true, label="nuc new")
plot_csv("plot_data/before_recompute_nucleus/v_beta_nonuc.csv", false, label="nonuc", ls="dotted")
PyPlot.axhline(0.0, lw=0.5, color="black")
PyPlot.legend()
PyPlot.title("Velocity vs f_β")

# channel width
PyPlot.subplot(246)
plot_csv("plot_data/before_recompute_nucleus/v_width_nuc_old.csv", true, label="nuc old")
plot_csv("plot_data/before_recompute_nucleus/v_width_nuc_new.csv", true, label="nuc new")
plot_csv("plot_data/before_recompute_nucleus/v_width_nonuc.csv", false, label="nonuc", ls="dotted")
PyPlot.axhline(0.0, lw=0.5, color="black")
PyPlot.legend()
PyPlot.title("Velocity vs f_width")

# channel pulsation
PyPlot.subplot(247)
plot_csv("plot_data/before_recompute_nucleus/v_omega_nuc_old.csv", true, label="nuc old")
plot_csv("plot_data/v_omega_nuc.csv", true, label="nuc")
PyPlot.axhline(0.0, lw=0.5, color="black")
PyPlot.legend()
PyPlot.title("Velocity vs f_ω")

PyPlot.show()
