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

function initialize(E::NamedTuple, I::NamedTuple, R0::Real;
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
        "nu" => 1/80,
        "foi" => 0
    )
    properties["mu_daily"] = 0 # TODO
    properties["nu_daily"] = 0 # TODO

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
        infect!(agent, model)
    end
    # initialize force of infection
    foi!(model)

    return model
end
initialize() = initialize((mean=8, sigma=3), (mean=8, sigma=3), 10)

function foi!(model)
    num_infections = 0
    for (id, agent) in model.agents
        if agent.status == :I
            num_infections += 1
        end
    end
    model.properties["foi"] = num_infections * model.properties["R0"] * model.properties["I"].mean / length(model.agents)
end

function agent_step!(agent, model)
    update!(agent, model)
end

function model_step!(model)
    births!(model)
    deaths!(model)
    foi!(model)
end

function infect!(agent, model)
    agent.status = :E
    agent.timer = round(typeof(agent.timer), max(0, model.properties["E"].sigma*randn() + model.properties["E"].mean))
end

function update!(agent, model)
    if agent.status == :S
        if rand < model.properties["foi"]
            infect!(agent, model)
        end
    elseif agent.status == :E
    elseif agent.status == :I
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

end

import .ModelC

model = ModelC.initialize()
