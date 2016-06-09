/**
 * << detailed description >>
 *
 * @file Vec3.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.04.16
 */

#ifndef READDY_MAIN_VEC3_H
#define READDY_MAIN_VEC3_H

#include <array>

namespace readdy {
    namespace model {
        class Vec3 {
        public:
            Vec3();
            Vec3(double x, double y, double z);
            Vec3& operator+=(const Vec3& rhs);
            Vec3& operator*=(const double a);
            double operator[](const unsigned int i) const;
            bool operator==(const Vec3& rhs);
            bool operator!=(const Vec3& rhs);

            friend std::ostream& operator<< (std::ostream&, const Vec3&);
            friend bool operator==(const Vec3& lhs, const Vec3& rhs);
        private:
            std::array<double, 3> data;
        };

        inline double operator*(const Vec3 &lhs, const Vec3 &rhs) {
            return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2];
        }

        inline Vec3 operator*(const Vec3 &lhs, const double &rhs) {
            return Vec3(rhs*lhs[0], rhs*lhs[1], rhs*lhs[2]);
        }

        inline Vec3 operator*(const double &rhs, const Vec3 &lhs) {
            return lhs * rhs;
        }
    }
}

#endif //READDY_MAIN_VEC3_H
