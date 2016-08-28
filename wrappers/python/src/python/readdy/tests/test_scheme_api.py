import unittest

from readdy._internal.api import Simulation


class TestSchemeApi(unittest.TestCase):
    def test_sanity(self):
        simulation = Simulation()
        simulation.set_kernel("SingleCPU")
        configurator = simulation.run_scheme_readdy(False)
        scheme = configurator \
            .with_integrator("EulerBDIntegrator") \
            .include_forces(False) \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .evaluate_observables(False) \
            .configure()
        scheme.run(10)

    def test_sanity_oneliner(self):
        simulation = Simulation()
        simulation.set_kernel("SingleCPU")

        simulation.run_scheme_readdy(False) \
            .with_integrator("EulerBDIntegrator") \
            .include_forces(False) \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .evaluate_observables(False) \
            .configure_and_run(10)

        simulation.run_scheme_readdy(False) \
            .with_integrator("EulerBDIntegrator") \
            .include_forces(False) \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .evaluate_observables(False) \
            .configure().run(10)


if __name__ == '__main__':
    unittest.main()
