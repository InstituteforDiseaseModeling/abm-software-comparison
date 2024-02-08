from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.space import MultiGrid
from mesa.datacollection import DataCollector
from tqdm import tqdm
from mesa.visualization.modules import CanvasGrid
from mesa.visualization.ModularVisualization import ModularServer



class SEIRSAgent(Agent):
    def __init__(self, unique_id, model):
        super().__init__(unique_id, model)
        self.state = "S"  # Initial state is Susceptible
        self.exposed_duration = 5  # Duration in days an agent remains exposed
        self.infected_duration = 7  # Duration in days an agent remains infectious
        self.days_infected = 0

    def move(self):
        possible_steps = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False
        )
        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)

    def step(self):
        if self.state == "E":
            self.days_infected += 1
            if self.days_infected >= self.exposed_duration:
                self.state = "I"

        elif self.state == "I":
            self.days_infected += 1
            if self.days_infected >= self.infected_duration:
                self.state = "R"
                self.days_infected = 0

    def infect_neighbors(self):
        neighbors = self.model.grid.get_neighbors(self.pos, moore=True, radius=1)
        for neighbor in neighbors:
            if (
                self.random.random() < self.model.infection_rate
                and neighbor.state == "S"
            ):
                neighbor.state = "E"

class SEIRSModel(Model):
    def __init__(self, N, num_nodes, width, height):
        grid = CanvasGrid(agent_portrayal, 10, 10, 500, 500)

        self.num_agents = int(N / num_nodes)
        self.schedule = RandomActivation(self)
        self.grid = MultiGrid(width, height, True)
        self.infection_rate = 0.2
        self.nodes = self.create_nodes(num_nodes)

        for i in range(self.num_agents):
            agent = SEIRSAgent(i, self)
            self.place_agent(agent)

        self.datacollector = DataCollector(
            agent_reporters={"State": "state"}
        )

        self.server = ModularServer(
            SEIRSModel,
            [self.grid],
            "SEIRS Model",
            {"N": 1e5, "num_nodes": 5, "width": 10, "height": 10}
        )

    def create_nodes(self, num_nodes):
        nodes = []
        for i in range(num_nodes):
            node_agents = [SEIRSAgent(j, self) for j in range(i * self.num_agents // num_nodes, (i + 1) * self.num_agents // num_nodes)]
            nodes.append(node_agents)
        return nodes

    def place_agent(self, agent):
        x, y = self.random.randint(0,self.grid.width-1), self.random.randint(0,self.grid.height-1)
        print( f"adding agent at {x},{y}" )
        self.grid.place_agent(agent, (x, y))
        self.schedule.add(agent)

    def step(self):
        self.datacollector.collect(self)
        self.schedule.step()

def agent_portrayal(agent):
    portrayal = {"Shape": "circle", "Color": "red" if agent.state == "I" else "white", "Filled": "true", "Layer": 0, "r": 0.5}
    return portrayal

# Example usage:
model = SEIRSModel(N=1e5, num_nodes=5, width=10, height=10)

model.server.port = 8521  # The default port
model.server.launch()

"""
for day in tqdm(range(365)):  # Run for 1 year (365 days)
    model.step()

# Access data collected during the simulation
agent_states = model.datacollector.get_agent_vars_dataframe().reset_index()
print(agent_states)
"""
