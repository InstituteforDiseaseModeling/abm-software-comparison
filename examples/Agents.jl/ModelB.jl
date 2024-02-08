"""
"Create an agent based SEIR disease model with 1 million agents distributed unevenly across 25 nodes with 
an initial outbreak of 100 agents. β and mean incubation period are inputs to the model. Run for 1 simulated 
year and plot the aggregated SEIR curves for the total population and the SEIR curves for each node in a 5 by 5 grid of plots."

This model aggregates the force of infection per node and calculates resulting infections.
"""
module ModelB

using Agents
using Agents.Graphs
using Random
using StatsBase
using CairoMakie


# Agent structure
@agent PoorSoul GraphAgent begin
    # id, pos, days_infected, status
    status::Symbol  # 1:S, 2:E, 3:I, 4:R
end

"""
Initialize the model by:
1. Setting up the spatial structure. 

"""
function initialize(β::AbstractFloat, σ::AbstractFloat,  γ::AbstractFloat;
    num_agents::Integer = 10^6, 
    initial_infections::Integer = 100, 
    seed::Integer = 42)

    rng = Xoshiro(seed)

    # Spatial setup (NxN grid)
    grid_dims = (5, 5)
    space = GraphSpace(grid(grid_dims))

    # Fix properties of the model using a named tuple
    properties = Dict(
        "β" => β,
        "γ" => γ,
        "σ" => σ,
        "migration_rate" => 0.001
        "foi" => zeros(nv(space.graph))
    )

    # Initialize the model
    model = UnremovableABM(PoorSoul, space; properties, rng)

    # Add initial individuals
    W = rand(nv(space.graph))
    W ./= sum(W)
    for p in 1:num_agents
        node = sample(1:length(W), Weights(W))
        add_agent!(node, model, :S) # Susceptible
    end

    # add initial outbreak
    # outbreak_node = rand(1:nv(space.graph))
    outbreak_node = 1
    inds = ids_in_position(outbreak_node, model)
    for n in 1:initial_infections
        agent = model[rand(inds)]
        agent.status = :I # Infected
    end

    return model
end

function agent_step!(agent, model)
    migrate!(agent, model)
    transmit!(agent, model)
    update!(agent, model)
end

function model_step!(model)
    foi!(model)
end

# function to count states
infectious(x) = count(i == :I for i in x)
susceptible(x) = count(i == :S for i in x)
exposed(x) = count(i == :E for i in x)
infectious(x) = count(i == :I for i in x)
recovered(x) = count(i == :R for i in x)

"""
Constant migration rate, could easily be incorporated into the graph edges
"""
function migrate!(agent, model)
    pid = agent.pos
    if rand() < model.properties["migration_rate"]
        agent.pos = rand(neighbors(model.space.graph, pid))
    end
end

"""
Transmit within the node
"""
function transmit!(agent, model)
    if agent.status == :I
        # agents in the node
        ids = ids_in_position(agent.pos, model)
        # force of infection
        foi = model.properties.β / length(ids)
        # loop over all other agents ($$$$)
        for id in ids
            if model.agents[id].status == :S
                if rand() < foi
                    model.agents[id].status = :E
                end
            end
        end
    end
end

"""
Update status
"""
function update!(agent, model)
    if agent.status == :E
        if rand() < model.properties.σ 
            agent.status = :I
        end        
    elseif agent.status == :I
        if rand() < model.properties.γ 
            agent.status = :R
        end
    end
    return
end

##########################################

function run_model(num_agents::Integer = 10^4)
    # parameters
    γ = 1/10 # recovery rate; 1 / average duration of infection
    β = 2.5 * γ # R_0 / gamma
    σ = 1/3 # incubation rate; 1 / average duration of incubation

    # Initialize
    model = initialize(β, σ, γ; num_agents=num_agents)

    # Data collection
    to_collect = [(:status, f) for f in (susceptible, exposed, infectious, recovered)]
    to_collect = [to_collect; [(:status, f, (a) -> a.pos == node) for node in 1:nv(model.space.graph) for f in (susceptible, exposed, infectious, recovered)]]

    # Run
    data, _ = run!(model, agent_step!, 100; adata=to_collect)

    return data
end

function make_plots(data)
    # Plot
    # 1: total number in each state
    fig1 = Figure()
    ax = Axis(fig1[1,1])
    lines!(ax, data.step, data.susceptible_status)
    lines!(ax, data.step, data.exposed_status)
    lines!(ax, data.step, data.infectious_status)
    lines!(ax, data.step, data.recovered_status)
    save("fig1.png", fig1)

    # 2: Plot in the grid
    ioff = 4 + 1
    fig2 = Figure()
    for i in 1:5
        for j in 1:5
            ax = Axis(fig2[i,j])
            for k in 1:4
                lines!(ax, data.step, data[:, ioff + (i-1)*(j-1)*4 + k])
            end
        end
    end
    fig2
    save("fig2.png",fig2)
end

end