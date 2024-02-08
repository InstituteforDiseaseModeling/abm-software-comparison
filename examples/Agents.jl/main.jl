using BenchmarkTools

import .ModelA

suite = BenchmarkGroup()

suite["A"] = BenchmarkGroup()

suite["A"][4] = @benchmarkable ModelA.run_model(10^4)

tune!(suite)
results = run(suite, verbose = true)

BenchmarkTools.save("output.json", median(results))