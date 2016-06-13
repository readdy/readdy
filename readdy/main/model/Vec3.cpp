/**
 * << detailed description >>
 *
 * @file Vec3.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.04.16
 */

#include <readdy/model/Vec3.h>
#include <assert.h>
#include <ostream>

namespace readdy {
    namespace model {
        Vec3 &Vec3::operator+=(const Vec3 &rhs) {
            data[0] += rhs.data[0];
            data[1] += rhs.data[1];
            data[2] += rhs.data[2];
            return *this;
        }

        Vec3 &Vec3::operator*=(double a) {
            data[0] *= a;
            data[1] *= a;
            data[2] *= a;
            return *this;
        }

        Vec3::Vec3(const std::array<double, 3> &xyz) {
            data = std::array<double, 3>(xyz);
        }

        Vec3::Vec3(double x, double y, double z) {
            data[0] = x;
            data[1] = y;
            data[2] = z;
        }

        double Vec3::operator[](const unsigned int i) const{
            assert(0 <= i && i < 3);
            return data[i];
        }

        Vec3::Vec3()  : Vec3(0, 0, 0){

        }

        bool Vec3::operator==(const Vec3 &rhs) {
            return data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2];
        }

        bool Vec3::operator!=(const Vec3 &rhs) {
            return !(data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2]);
        }
        Vec3 operator+(const Vec3& lhs, const Vec3& rhs) {
            return {lhs[0]+rhs[0], lhs[1]+rhs[1], lhs[2]+rhs[2]};
        }

        Vec3 operator-(const Vec3& lhs, const Vec3& rhs) {
            return lhs + (-1*rhs);
        }

        std::ostream& operator<< (std::ostream& os, const Vec3& vec) {
            os << "Vec3("<< vec[0] << ", "<< vec[1] << ", " << vec[2] << ")";
            return os;
        }

        bool operator==(const Vec3& lhs, const Vec3& rhs) {
            return lhs[0] == rhs[0] && lhs[1] == rhs[1] && lhs[2] == rhs[2];
        }
    }
}













