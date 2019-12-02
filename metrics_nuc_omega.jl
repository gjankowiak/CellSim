push!(LOAD_PATH, "src")

import CellSim
import Dates

# The configuration file from which parameters are taken
base_config = "configs_metrics/nucleus/varia_omega_param.yaml"

# The output prefix, this overrides the one set in the configuration file
output_prefix = "nuc_omega"

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
omegas = collect(range(2.0, 12.0; length=11))

# Loop
for omega in omegas
    # Set the parameter value
    P.f_ω0 = omega

    # Write the current value to the log file
    write(param_log, string(omega, "\n"))

    # Run the simulation
    CellSim.launch(P, F, config)
end

close(param_log)
