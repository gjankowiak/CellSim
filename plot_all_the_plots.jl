import DelimitedFiles
import PyPlot

function plot_csv(file, l=false; args...)
    println(file)

    d, h = DelimitedFiles.readdlm(file, ','; header=true)

    sorted_idx = sortperm(d[:,1])

    println(d)

    PyPlot.plot(d[sorted_idx,1], d[sorted_idx,2]; args...)

    if l
        PyPlot.xlabel(h[1])
        PyPlot.ylabel(h[2])
    end
end


PyPlot.figure()


# kb
PyPlot.subplot(231)
plot_csv("plot_data/v_kb_nuc_new.csv", true, label="nuc new μn = 20")
plot_csv("plot_data/v_kb_nuc_new_mun60.csv", true, label="nuc new µn=60")
plot_csv("plot_data/v_kb_nuc_old.csv", true, label="nuc old μn = 20")
plot_csv("plot_data/v_kb_nuc_old_mun60.csv", true, label="nuc old μn=60")
plot_csv("plot_data/v_kb_nuc_old_mun100.csv", true, label="nuc old µn=100")
PyPlot.axhline(0.0, lw=0.5)
PyPlot.title("Velocity vs kb")

# channel depth
PyPlot.subplot(234)
plot_csv("plot_data/v_beta_nuc_old.csv", true, label="nuc old")
plot_csv("plot_data/v_beta_nonuc.csv", false, label="nonuc")
PyPlot.axhline(0.0, lw=0.5)
PyPlot.legend()
PyPlot.title("Velocity vs f_β")

# channel width
PyPlot.subplot(235)
plot_csv("plot_data/v_width_nuc_old.csv", true, label="nuc old")
plot_csv("plot_data/v_width_nuc_new.csv", true, label="nuc new")
plot_csv("plot_data/v_width_nonuc.csv", false, label="nonuc")
PyPlot.axhline(0.0, lw=0.5)
PyPlot.legend()
PyPlot.title("Velocity vs f_width")

# channel pulsation
PyPlot.subplot(236)
plot_csv("plot_data/v_omega_nuc_old.csv", true, label="nuc old")
plot_csv("plot_data/v_omega_nuc_new.csv", true, label="nuc new")
PyPlot.axhline(0.0, lw=0.5)
PyPlot.legend()
PyPlot.title("Velocity vs f_ω")

PyPlot.show()
