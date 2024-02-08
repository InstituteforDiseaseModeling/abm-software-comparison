using BenchmarkTools
include("ModelA.jl")
import .ModelA
include("ModelB.jl")
import .ModelB

suite = BenchmarkGroup()

suite["A"] = BenchmarkGroup()

suite["A"][4] = @benchmarkable ModelA.run_model(10^4)
suite["A"][6] = @benchmarkable ModelA.run_model(10^6)

suite["B"] = BenchmarkGroup()
suite["B"][4] = @benchmarkable ModelB.run_model(10^4)
suite["B"][6] = @benchmarkable ModelB.run_model(10^6)
tune!(suite)
results = run(suite, verbose = true)

BenchmarkTools.save("output.json", median(results))