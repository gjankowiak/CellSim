import JankoUtils

d_depth = 0.1
d_w0 = 8
d_width = 0.4

width = d_width
depth = d_depth
prefix="w0"

w0_range = collect(range(4, stop=12; step=0.25))

for (i,w0) in enumerate(w0_range)
    outputfile = "autoconfigs/metrics_$(width)_$(depth)_$(w0).yaml"
    run(pipeline(`sed -e s/§§prefix§§/$prefix/ -e s/§§width§§/$width/ -e s/§§w0§§/$w0/ -e s/§§depth§§/$depth/ metrics.yaml.in`, stdout=outputfile))
    go(outputfile)
end
