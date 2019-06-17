push!(LOAD_PATH, "src")

import Plotting
import Metrics

filename = ARGS[1]
path = join(split(ARGS[1], "/")[1:end-1], "/")
tail = split(ARGS[1], "/")[end]
prefix = join(split(tail, "_")[1:2], "_")
full_prefix = string(path, "/", prefix)

println("prefix: ", full_prefix)

m = Metrics.read_metrics(full_prefix)
Plotting.plot_metrics(m)
