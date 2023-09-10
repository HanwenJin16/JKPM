include(raw"../src/JKPM.jl")
using .JKPM

using StatsBase

# Your partition as an example
partition = [2, 2, 3, 1, 1, 3, 3, 2, 2]

# Count the atoms in each group
group_counts = StatsBase.countmap(partition)

println(group_counts)  # This will output something like: Dict(2 => 4, 3 => 3, 1 => 2)