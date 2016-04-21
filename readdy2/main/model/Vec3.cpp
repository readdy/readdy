#include <readdy/model/Vec3.h>

/**
 * << detailed description >>
 *
 * @file Vec3.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.04.16
 */



readdy::model::Vec3 &readdy::model::Vec3::operator+=(const Vec3 &rhs) {
    data[0] += rhs.data[0];
    data[1] += rhs.data[1];
    data[2] += rhs.data[2];
    return *this;
}

readdy::model::Vec3 &readdy::model::Vec3::operator*=(double a) {
    data[0] *= a;
    data[1] *= a;
    data[2] *= a;
    return *this;
}

readdy::model::Vec3::Vec3(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

double readdy::model::Vec3::operator[](const size_t i) const{
    return data[i];
}











