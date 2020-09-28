push!(LOAD_PATH, "src")

import CellSim
import Dates

# The configuration file from which parameters are taken
base_config = "configs_metrics/nucleus/varia_width_param.yaml"

# The output prefix, this overrides the one set in the configuration file
output_prefix = "nuc_width"

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
# If it already exists, we do not delete it, but write at the end
param_log = open(parameter_log_file, "a")
date = string(Dates.now())
write(param_log, string("\n", date, "\n\n"))

# The range in which the varying parameter is taken
ws = range(0.3, 0.8; length=20)[10:end]

# Loop
for w in ws
    # Set the parameter value
    P.f_width = w

    # Write the current value to the log file
    write(param_log, string(w, "\n"))

    # Run the simulation
    CellSim.launch(P, F, config)
end

close(param_log)
