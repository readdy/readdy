//
// Created by mho on 6/18/18.
//

#include <pybind11/pybind11.h>
#include <readdy/model/potentials/PotentialOrder1.h>
#include <readdy/model/potentials/PotentialOrder2.h>

namespace py = pybind11;

class TestPotential : public readdy::model::potentials::PotentialOrder1{
public:
    explicit TestPotential(particle_type_type ptype) : PotentialOrder1(ptype) {}

    std::string describe() const override {
        return "Test potential of order 1";
    }

    std::string type() const override {
        return "Test potential";
    }

    readdy::scalar calculateEnergy(const readdy::Vec3 &position) const override {
        return 0;
    }

    void calculateForce(readdy::Vec3 &force, const readdy::Vec3 &position) const override {
        force += {.1, .1, .1};
    }

};

class TestPairPotential : public readdy::model::potentials::PotentialOrder2 {
public:
    TestPairPotential(particle_type_type type1, particle_type_type type2, readdy::scalar cutoff)
            : PotentialOrder2(type1, type2), cutoff(cutoff) {}

    std::string describe() const override {
        return "Test potential of order 2";
    }

    std::string type() const override {
        return "Test pair potential";
    }

    readdy::scalar calculateEnergy(const readdy::Vec3 &x_ij) const override {
        return 2*x_ij * x_ij;
    }

    void calculateForce(readdy::Vec3 &force, const readdy::Vec3 &x_ij) const override {
        const auto distSquared = x_ij*x_ij;
        force -= distSquared * readdy::Vec3{1., 1., 1.};
    }

    readdy::scalar getCutoffRadiusSquared() const override {
        return cutoff*cutoff;
    }
private:
    readdy::scalar cutoff;
};

PYBIND11_MODULE (custom_potential_example, m) {
    py::module::import("readdy");
    py::class_<TestPotential, readdy::model::potentials::PotentialOrder1>(m, "TestPotential")
            .def(py::init<readdy::particle_type_type>());
    py::class_<TestPairPotential, readdy::model::potentials::PotentialOrder2>(m, "TestPairPotential")
            .def(py::init<readdy::particle_type_type, readdy::particle_type_type, readdy::scalar>());
}
