# CellSim

## Installation and first run

- Install Julia 1.0 or newer
- Clone this repository:

    `$ git clone https://gitea.oknaj.eu/gjankowiak/CellSim`

- Switch directories and branch:

    `$ cd CellSim`
    `$ git checkout nucleus`

- Download necessary packages: on the Julia prompt, hit `]`, you should get a prompt:

    `pkg>`

You can now do

    `pkg> activate .`

and

    `(CellSim) pkg> instantiate`

- At this point you should be able to run the code with

    `$ julia --project run.jl configs/straight.yaml`

It can take a really long time (30 seconds or more) to start because of the initial loading required.

## Configuration files

The models parameters are defined in YAML files located in the `configs` directory. You can of course copy and edit those. There are a large number of parameters, which are explained in short comments.
