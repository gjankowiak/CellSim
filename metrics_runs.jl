push!(LOAD_PATH, "src")

import CellSim

base_config = "configs_metrics/without_nucleus/varia_width_param.yaml"

output_prefix = "nonuc_width"

P, F, config = CellSim.read_config(base_config)

config["output_prefix"] = output_prefix

ws = range(4, 14; length=20)

for w in ws
    P.f_ω0 = w
    CellSim.launch(P, F, config)
end
