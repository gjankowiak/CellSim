push!(LOAD_PATH, "src")

import CellSim
import Dates

# The configuration file from which parameters are taken
base_config = "configs_metrics/no_nucleus/varia_beta_param.yaml"

# The output prefix, this overrides the one set in the configuration file
output_prefix = "nonuc_beta"

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
betas = collect(range(0.0, 0.3; length=21))[2:end]

# Loop
for beta in betas
    # Set the parameter value
    P.f_β = beta

    # Write the current value to the log file
    write(param_log, string(beta, "\n"))

    # Run the simulation
    CellSim.launch(P, F, config)
end

close(param_log)
