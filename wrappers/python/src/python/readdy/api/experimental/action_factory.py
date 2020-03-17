# todo some argument checking
class ActionFactory(object):
    def __init__(self, sim):
        self._sim = sim

    def integrator_euler_brownian_dynamics(self, time_step):
        return self._sim.create_action_euler_bd(time_step)

    def calculate_forces(self):
        return self._sim.create_action_calculate_forces()

    def create_neighbor_list(self, interaction_distance):
        return self._sim.create_action_create_neighbor_list(interaction_distance)

    def update_neighbor_list(self):
        return self._sim.create_action_update_neighbor_list()

    def clear_neighbor_list(self):
        return self._sim.create_action_clear_neighbor_list()

    def reaction_handler_uncontrolled_approximation(self, time_step):
        return self._sim.create_action_uncontrolled_approximation(time_step)

    def reaction_handler_gillespie(self, time_step):
        return self._sim.create_action_gillespie(time_step)

    def reaction_handler_detailed_balance(self, time_step):
        return self._sim.create_action_detailed_balance(time_step)

    def topology_reaction_handler(self, time_step):
        return self._sim.create_action_evaluate_topology_reactions(time_step)

    def break_bonds(self, time_step, break_config):
        return self._sim.create_action_break_bonds(time_step, break_config)

    def evaluate_observables(self):
        return self._sim.create_action_evaluate_observables()

    def make_checkpoint(self, base_path, max_n_saves):
        return self._sim.create_action_make_checkpoint(base_path, max_n_saves)
