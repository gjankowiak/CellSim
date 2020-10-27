push!(LOAD_PATH, "src")

import CellSim
import Dates

# The configuration file from which parameters are taken
# base_config = "configs_metrics/nucleus/varia_mu.yaml"
base_config = "configs_metrics/nucleus/varia_mu.yaml"

# The output prefix, this overrides the one set in the configuration file
output_prefix = "nuc_mu_kb_fine_1_test"

# Load configuration
P, F, config = CellSim.read_config(base_config)

# Override the output prefix and create it
config["output_prefix"] = output_prefix
mkpath(output_prefix)

# Set batch mode (only closes Figures at the end to avoid clutter)
config["batch"] = true

# Log file, used to list parameter value so we don't have to
# compute them by hand. We cannot use the YAML files as Julia can
# only read YAML files, not write to them.

# Log file name
parameter_log_file = string(output_prefix, "/parameter_log.txt")

# Open the log file and write the date on top
param_log = open(parameter_log_file, "a")
date = string(Dates.now())
write(param_log, string("\n", date, "\n\n"))

# The range in which the varying parameter is taken
# mus = 10 .^range(-4, log10.(400); length=20)
# mus = collect(range(1, stop=50, length=15))
# kbs = 10 .^range(log10.(5e-4), log10.(0.5); length=15)
#
mus = collect(range(1, stop=50, length=15))[11]
kbs = 10 .^range(log10.(5e-4), log10.(0.5); length=15)[10]

# fine 1
mus = collect(range(23.6, stop=34, length=8))[1:1]
kbs = 10 .^range(log10.(0.0204), log10.(0.035); length=4)[1:1]

# fine 2
# mus = collect(range(34.3, stop=50, length=6))
# kbs = 10 .^range(log10.(0.035), log10.(0.0575); length=4)

param_range = Base.Iterators.product(kbs, mus)
# param_range = Base.Iterators.drop(param_range, 72)

catch_errors = false

# Loop
for (kb, mu) in param_range
    println("kb: ", kb, ", mu: ", mu)

    # Set the parameter value
    P.N_mu = mu
    P.N_kb = kb

    # Write the current value to the log file
    date = string(Dates.now())
    write(param_log, string(date, "\t", mu, "\t", kb, "\n"))

    # Run the simulation
    if !catch_errors
        CellSim.launch(P, F, config; force_date_string=date)
    else
        try
            CellSim.launch(P, F, config; force_date_string=date)
        catch
            println("error")
        end
    end
end

close(param_log)
