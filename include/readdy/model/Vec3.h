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
            Vec3(const std::array<double, 3> &xyz);
            Vec3& operator+=(const Vec3& rhs);
            Vec3& operator*=(const double a);
            double operator[](const unsigned int i) const;
            double& operator[](const unsigned int i);
            bool operator==(const Vec3& rhs);
            bool operator!=(const Vec3& rhs);

            friend std::ostream& operator<< (std::ostream&, const Vec3&);
            friend bool operator==(const Vec3& lhs, const Vec3& rhs);
            friend Vec3 operator+(const Vec3& lhs, const Vec3& rhs);
            friend Vec3 operator-(const Vec3& lhs, const Vec3& rhs);

        private:
            std::array<double, 3> data;

            // todo: do this here?≠
            //const std::array<double, 3> *const boxSize;
            //const std::array<bool, 3> *const periodic;
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

        template<bool PX, bool PY, bool PZ, typename... Args>
        double dist(const Vec3& lhs, const Vec3& rhs, Args... args);

        template<bool PX, bool PY, bool PZ>
        inline double dist(const Vec3& lhs, const Vec3& rhs, const double& dx, const double &dy, const double &dz) {
            // todo: assuming lhs and rhs are within box already. compare: https://en.wikipedia.org/wiki/Periodic_boundary_conditions
            auto dv = rhs - lhs;
            if(PX) {
                if(dv[0] > dx * .5) dv[0] -= dx;
                if(dv[0] <= dx * .5) dv[0] += dx;
            }
            if(PY) {
                if(dv[1] > dy * .5) dv[1] -= dy;
                if(dv[1] <= dy * .5) dv[1] += dy;
            }
            if(PZ) {
                if(dv[2] > dz * .5) dv[2] -= dz;
                if(dv[2] <= dz * .5) dv[2] += dz;
            }
            return dv * dv;
        };

        template<>
        inline double dist<false, false, false>(const Vec3& lhs, const Vec3 &rhs) {
            return (lhs - rhs)*(lhs-rhs);
        };
    }
}

#endif //READDY_MAIN_VEC3_H
