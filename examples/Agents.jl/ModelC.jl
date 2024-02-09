"""
Single node simulation
- 1e6 agents
- 35 years 
- sigma = 1/8 (exposure recovery)
- gamma = 1/8 (infectious recovery)
- R0 = 10
- mu = 0.03 kpop (crude birth rate)
- nu = 1/80 / pop(mortality rate)
"""

module ModelC
using Agents

struct DiseaseParams
    sigma::AbstractFloat
    gamma::AbstractFloat
    mu::AbstractFloat
    nu::AbstractFloat
    R0::AbstractFloat
end
DiseaseParams() = DiseaseParams(1/8,1/8,0.03,1/80,10)

# Agent structure
@agent PoorSoul AbstractAgent begin
    # id, days_infected, status
    status::Symbol  # 1:S, 2:E, 3:I, 4:R
end

end

