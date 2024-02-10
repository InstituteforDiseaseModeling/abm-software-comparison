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
using Random

# Agent structure
mutable struct PoorSoul <: AbstractAgent    # id, days_infected, status, timer
    id::Integer
    status::Symbol  # 1:S, 2:E, 3:I, 4:R
    timer::Integer
end

abstract type AbstractDiseaseParams end
struct ExposureDiseaseParams <: AbstractDiseaseParams
    mean::Real
    sigma::Real
end
struct InfectiousDiseaseParams <: AbstractDiseaseParams
    mean::Real
    sigma::Real
end

function initialize(E, I, R0::Real;
    num_agents::Integer = 10^6, 
    initial_infections::Integer = 100, 
    seed::Integer = 42)

    rng = Xoshiro(seed)

    # Fix properties of the model using a named tuple
    properties = Dict(
        "E" => E,
        "I" => I,
        "R0" => R0,
        "mu" => 0.03,
        "nu" => 0.0125,
        "foi" => 0
    )
    properties["mu_daily"] = log(1 + properties["mu"]) / 365
    properties["nu_daily"] = log(1 - properties["nu"]) / 365

    # Initialize the model
    model = ABM(PoorSoul; properties, rng)

    # Add initial susceptible individuals
    for p in 1:num_agents
        add_agent!(model, :S, 0) # Susceptible
    end

    # initialize recovered fraction
    num_recovered = Integer(num_agents - num_agents / model.properties["R0"])
    # S = np.uint32(np.round(params.pop_size / params.r_naught))  # * 1.0625
    for n in 1:num_recovered
        agent = model[rand(1:num_agents)]
        agent.status = :R # Recovered
    end

    # add initial outbreak
    for n in 1:initial_infections
        agent = model[rand(1:num_agents)]
        expose!(agent, model)
    end
    # initialize force of infection
    foi!(model)

    return model
end
initialize() = initialize(ExposureDiseaseParams(8, 1), InfectiousDiseaseParams(8, 1), 10)
initialize(; num_agents) = initialize(ExposureDiseaseParams(8, 1), InfectiousDiseaseParams(8, 1), 10, num_agents=num_agents)

function foi!(model)
    num_infections = 0
    for (id, agent) in model.agents
        if agent.status == :I
            num_infections += 1
        end
    end
    model.properties["foi"] = max(1, num_infections) * model.properties["R0"] / model.properties["I"].mean / length(model.agents)
end

function agent_step!(agent, model)
    update!(agent, model)
end

function model_step!(model)
    births!(model)
    deaths!(model)
    foi!(model)
end

expose!(agent, model) = start_timer!(agent, model.properties["E"])
function start_timer!(agent, params::ExposureDiseaseParams)
    agent.status = :E
    agent.timer = round(typeof(agent.timer), max(0, params.sigma*randn() + params.mean))
end

infectious!(agent, model) = start_timer!(agent, model.properties["I"])
function start_timer!(agent, params::InfectiousDiseaseParams)
    agent.status = :I
    agent.timer = round(typeof(agent.timer), max(0, params.sigma*randn() + params.mean))
end

recover!(agent, model) = agent.status = :R

function update!(agent, model)
    if agent.status == :S
        if rand() < model.properties["foi"]
            expose!(agent, model)
        end
    elseif agent.status == :E
        if agent.timer == 0
            infectious!(agent, model)
        else
            agent.timer -= 1
        end
    elseif agent.status == :I
        if agent.timer == 0
            recover!(agent, model)
        else
            agent.timer -= 1
        end
    end
end

function births!(model)
    num_births = round(Int, model.properties["mu_daily"] * length(model.agents))
    for i in 1:num_births
        add_agent!(model, :S, 0)
    end
end

function deaths!(model)
    num_deaths = round(Int, model.properties["nu_daily"] * length(model.agents))
    for i in 1:num_deaths
        agent = rand(1:length(model.agents))
        remove_agent!(agent, model)
    end
end
#########################

function run_model(;num_days::Integer=365, num_agents::Integer = 10^4)

    # Initialize
    model = initialize(num_agents=num_agents)

    # function to count states
    susceptible(x) = count(i == :S for i in x)
    exposed(x) = count(i == :E for i in x)
    infectious(x) = count(i == :I for i in x)
    recovered(x) = count(i == :R for i in x)

    # Data collection
    to_collect = [(:status, f) for f in (susceptible, exposed, infectious, recovered)]

    # Run
    data, _ = run!(model, agent_step!, model_step!, num_days; adata=to_collect)

    return data
end

end

import .ModelC
using CairoMakie

@time data = ModelC.run_model(num_days=365, num_agents=10^6)
# blue, gold, green, purple
fig1 = Figure()
ax = Axis(fig1[1,1])
lines!(ax, data.step, data.susceptible_status, label="susceptble")
lines!(ax, data.step, data.exposed_status, label="exposure")
lines!(ax, data.step, data.infectious_status, label="infectious")
# lines!(ax, data.step, data.recovered_status, label="recovered")
fig1
# save("fig1.png", fig1)