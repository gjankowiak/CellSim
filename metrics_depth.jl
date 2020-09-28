import JankoUtils

d_depth = 0.1
d_w0 = 8
d_width = 0.4

width = d_width
depth = d_depth
w0 = d_w0
prefix="depth"

for dd in range(0, stop=0.7; step=0.0125)
    width = 0.3 + dd
    depth = dd
    outputfile = "autoconfigs/metrics_$(width)_$(depth)_$(w0).yaml"
    run(pipeline(`sed -e s/§§prefix§§/$prefix/ -e s/§§width§§/$width/ -e s/§§w0§§/$w0/ -e s/§§depth§§/$depth/ metrics.yaml.in`, stdout=outputfile))
    go(outputfile)
end

