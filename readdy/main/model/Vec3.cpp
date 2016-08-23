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

        Vec3 &Vec3::operator*=(const double a) {
            data[0] *= a;
            data[1] *= a;
            data[2] *= a;
            return *this;
        }

        Vec3 &Vec3::operator/=(const double a) {
            data[0] /= a;
            data[1] /= a;
            data[2] /= a;
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

        double Vec3::operator[](const unsigned int i) const {
            assert(0 <= i && i < 3);
            return data[i];
        }

        double &Vec3::operator[](const unsigned int i) {
            assert(0 <= i && i < 3);
            return data[i];
        }

        Vec3::Vec3() : Vec3(0, 0, 0) {

        }

        bool Vec3::operator==(const Vec3 &rhs) const {
            return data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2];
        }

        bool operator<=(const Vec3 &lhs, const Vec3 &rhs) {
            return lhs[0] <= rhs[0] && lhs[1] <= rhs[1] && lhs[2] <= rhs[2];
        }

        bool operator>(const Vec3 &lhs, const Vec3 &rhs) {
            return !(lhs <= rhs);
        }

        bool operator>=(const Vec3 &lhs, const Vec3 &rhs) {
            return lhs.data[0] >= rhs.data[0] && lhs.data[1] >= rhs.data[1] && lhs.data[2] >= rhs.data[2];
        }

        bool operator<(const Vec3 &lhs, const Vec3 &rhs) {
            return !(lhs >= rhs);
        }

        bool Vec3::operator!=(const Vec3 &rhs) const {
            return !(data[0] == rhs[0] && data[1] == rhs[1] && data[2] == rhs[2]);
        }

        Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs) {
            return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
        }

        Vec3 operator/(const Vec3 &lhs, const double rhs) {
            return {lhs[0] / rhs, lhs[1] / rhs, lhs[2] / rhs};
        }

        Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs) {
            return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
        }

        std::ostream &operator<<(std::ostream &os, const Vec3 &vec) {
            os << "Vec3(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
            return os;
        }

        void fixPosition(Vec3 &vec, const std::array<bool, 3> &periodic, const std::array<double, 3> &boxSize) {
            if (periodic[0]) {
                if (periodic[1]) {
                    if (periodic[2]) {
                        fixPosition<true, true, true>(vec, boxSize[0], boxSize[1], boxSize[2]);
                    } else {
                        fixPosition<true, true, false>(vec, boxSize[0], boxSize[1], boxSize[2]);
                    }
                } else {
                    if (periodic[2]) {
                        fixPosition<true, false, true>(vec, boxSize[0], boxSize[1], boxSize[2]);
                    } else {
                        fixPosition<true, false, false>(vec, boxSize[0], boxSize[1], boxSize[2]);
                    }
                }
            } else {
                if (periodic[1]) {
                    if (periodic[2]) {
                        fixPosition<false, true, true>(vec, boxSize[0], boxSize[1], boxSize[2]);
                    } else {
                        fixPosition<false, true, false>(vec, boxSize[0], boxSize[1], boxSize[2]);
                    }
                } else {
                    if (periodic[2]) {
                        fixPosition<false, false, true>(vec, boxSize[0], boxSize[1], boxSize[2]);
                    } else {
                        fixPosition<false, false, false>(vec, boxSize[0], boxSize[1], boxSize[2]);
                    }
                }
            }
        }

        double distSquared(const Vec3 &lhs, const Vec3 &rhs, const std::array<bool, 3> &periodic, const std::array<double, 3> &boxSize) {
            if (periodic[0]) {
                if (periodic[1]) {
                    if (periodic[2]) {
                        return distSquared<true, true, true>(lhs, rhs, boxSize[0], boxSize[1], boxSize[2]);
                    } else {
                        return distSquared<true, true, false>(lhs, rhs, boxSize[0], boxSize[1], boxSize[2]);
                    }
                } else {
                    if (periodic[2]) {
                        return distSquared<true, false, true>(lhs, rhs, boxSize[0], boxSize[1], boxSize[2]);
                    } else {
                        return distSquared<true, false, false>(lhs, rhs, boxSize[0], boxSize[1], boxSize[2]);
                    }
                }
            } else {
                if (periodic[1]) {
                    if (periodic[2]) {
                        return distSquared<false, true, true>(lhs, rhs, boxSize[0], boxSize[1], boxSize[2]);
                    } else {
                        return distSquared<false, true, false>(lhs, rhs, boxSize[0], boxSize[1], boxSize[2]);
                    }
                } else {
                    if (periodic[2]) {
                        return distSquared<false, false, true>(lhs, rhs, boxSize[0], boxSize[1], boxSize[2]);
                    } else {
                        return distSquared<false, false, false>(lhs, rhs, boxSize[0], boxSize[1], boxSize[2]);
                    }
                }
            }
        }
    }
}













