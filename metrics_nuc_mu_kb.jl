using Distributed

@everywhere push!(LOAD_PATH, "src")

@everywhere import Pkg
@everywhere Pkg.activate(".")

@everywhere import CellSim
@everywhere import CellSimCommon

@everywhere import Dates

function trigrid(xmin::Real, xmax::Real, nx::Int64,
                 ymin::Real, ymax::Real, ny::Int64, as_tuples=true)

    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    x_aligned = collect(range(xmin, step=dx, length=nx+1))
    x_non_aligned = collect(range(xmin+0.5*dx, step=dx, length=nx))

    d, r = divrem(ny, 2)

    n = (2nx+1)*(d+1) - (1-r)*nx

    coords = zeros(n, 2)

    k = 1

    # remaining rows
    for j in 0:ny

        if j % 2 == 0
            coords[k:k+nx, 1]  = x_aligned
            coords[k:k+nx, 2] .= ymin + j*dy
            k += nx + 1
        else
            coords[k:k+nx-1, 1] = x_non_aligned
            coords[k:k+nx-1, 2] .= ymin + j*dy
            k += nx
        end
    end

    if as_tuples
        t = [(coords[i,1], coords[i,2]) for i in 1:size(coords, 1)]
        return t
    end
    return coords
end

# The configuration file from which parameters are taken
# base_config = "configs_metrics/nucleus/varia_mu.yaml"
base_config = "configs_metrics/nucleus/varia_mu.yaml"

# The output prefix, this overrides the one set in the configuration file
output_prefix = "/scratch/scratch/jankowiak/cellsim_results/nuc_mu_kb_fine_1-1/"

# Load configuration
P_global, F_global, config_global = CellSim.read_config(base_config)

# Override the output prefix and create it
config_global["output_prefix"] = output_prefix
mkpath(output_prefix)

# Set batch mode (only closes Figures at the end to avoid clutter)
config_global["batch"] = true

# Log file, used to list parameter value so we don't have to
# compute them by hand. We cannot use the YAML files as Julia can
# only read YAML files, not write to them.

# Log file name
# parameter_log_file = string(output_prefix, "/parameter_log.txt")

# Open the log file and write the date on top
# param_log = open(parameter_log_file, "a")
# date = string(Dates.now())
# write(param_log, string("\n", date, "\n\n"))

mus = collect(range(1, stop=50, length=30))
kbs = 10 .^range(log10.(5e-4), log10.(0.5); length=30)

mu_min = 1
mu_max = 50
mu_n = 30

kb_min = 5e-4
kb_max = 0.5
kb_n = 30

param_range = trigrid(mu_min, mu_max, mu_n,
            log10.(kb_min), log10.(kb_max), kb_n)
param_range = [(p[1], 10 .^p[2]) for p in param_range]

# param_range = collect(Base.Iterators.product(kbs, mus))
# param_range = Base.Iterators.drop(param_range, 72)

catch_errors = true

# Loop
@everywhere function mukb_to_velocity(t::Tuple{Real,Real})
    (mu, kb) = t
    println("kb: ", kb, ", mu: ", mu)

    P = CellSimCommon.copy($P_global)
    config = deepcopy($config_global)

    # Set the parameter value
    P.N_mu = mu
    P.N_kb = kb

    # Write the current value to the log file
    date = string(Dates.now(), "-", myid())
    # write(param_log, string(date, "\t", mu, "\t", kb, "\n"))

    # Run the simulation
    got_error = false
    try
        metrics = CellSim.launch(P, $F_global, config; force_date_string=date)
        # println(P)
        # println($F_global)
        # println(config)
        # println(date)
        # metrics = Dict([("average_velocity", 0.5), ("status","finished")])
        return (mu=mu, kb=kb, metrics=metrics, date=date)
    catch e
        println("error")
        got_error = true
        fn = "$(config["output_prefix"])Run_$(date)_error.txt"
        f = open(fn, "w")
        write(f, string(e))
        write(f, "\n")
        close(f)
        return (mu=mu, kb=kb, metrics=nothing, date=date)
    end
end

results = pmap(mukb_to_velocity, param_range)

parameter_log_file = string(output_prefix, "/parameter_log.txt")
f = open(parameter_log_file, "w")
for l in results
    if !isnothing(l.metrics)
        av = get(l.metrics, "average_velocity", -1)
        status = l.metrics["status"]
    else
        av = -2
        status = "error"
    end
    write(f, string(l.date, ", ", l.mu, ", ", l.kb, ", ", av, ", ", status, "\n"))
end
close(f)

println(parameter_log_file)


# close(param_log)
