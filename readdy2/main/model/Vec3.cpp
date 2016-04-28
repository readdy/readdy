#include <readdy/model/Vec3.h>
#include <assert.h>

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

double readdy::model::Vec3::operator[](const unsigned int i) const{
    assert(0 <= i && i < 3);
    return data[i];
}

readdy::model::Vec3::Vec3()  : Vec3(0, 0, 0){

}

bool readdy::model::Vec3::operator==(const Vec3 &rhs) {
    return data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2];
}

bool readdy::model::Vec3::operator!=(const Vec3 &rhs) {
    return !(data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2]);
}















