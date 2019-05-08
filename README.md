# CellSim

## Installation and first run

- Install Julia 1.0 or newer
- Clone this repository:
```
$ git clone https://gitea.oknaj.eu/gjankowiak/CellSim
```

- Switch directories and branch:

```
    $ cd CellSim
    $ git checkout nucleus
```

- Download necessary packages: on the Julia prompt, hit `]`, you should get a prompt:

```
    pkg>
```
- You can now do

```
    pkg> activate .
```
and

```
    (CellSim) pkg> instantiate
```
- At this point you should be able to run the code with

```
    $ julia --project run.jl configs/straight.yaml
```
It can take a really long time (30 seconds or more) to start because of the initial loading required.

## Configuration files

The models parameters are defined in YAML files located in the `configs` directory. You can of course copy and edit those. There are a large number of parameters, which are explained in short comments.

## Troubleshooting

If you get an error like the following:

```
[NSApplication _setup:]
: unrecognized selector sent to instance 0x7f9183705170
2019-01-30 17:23:30.788 julia[2002:154898] *** Terminating app due to uncaught exception 'NSInvalidArgumentException', reason: '-[NSApplication _setup:]: unrecognized selector sent to instance 0x7f9183705170'

in expression starting at /Users/chiaragiverso/CellSim/run.jl:5
__pthread_kill at /usr/lib/system/libsystem_kernel.dylib (unknown line)
Allocations: 104904598 (Pool: 104879133; Big: 25465); GC: 227
Abort trap: 6
```

you can try installing PyPlot via julia directly:
```
    julia> ]
    pkg> add Conda
    julia> 
```
