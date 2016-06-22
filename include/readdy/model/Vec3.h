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
#include <math.h>

namespace readdy {
    namespace model {
        class Vec3 {
        public:
            Vec3();

            Vec3(double x, double y, double z);

            Vec3(const std::array<double, 3> &xyz);

            Vec3 &operator+=(const Vec3 &rhs);

            Vec3 &operator*=(const double a);

            Vec3 &operator/=(const double a);

            double operator[](const unsigned int i) const;

            double &operator[](const unsigned int i);

            bool operator==(const Vec3 &rhs);

            bool operator!=(const Vec3 &rhs);

            friend std::ostream &operator<<(std::ostream &, const Vec3 &);

            friend bool operator==(const Vec3 &lhs, const Vec3 &rhs);

            friend Vec3 operator+(const Vec3 &lhs, const Vec3 &rhs);

            friend Vec3 operator-(const Vec3 &lhs, const Vec3 &rhs);

        private:
            std::array<double, 3> data;
        };

        inline double operator*(const Vec3 &lhs, const Vec3 &rhs) {
            return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
        }

        inline Vec3 operator*(const Vec3 &lhs, const double rhs) {
            return Vec3(rhs * lhs[0], rhs * lhs[1], rhs * lhs[2]);
        }

        inline Vec3 operator*(const double rhs, const Vec3 &lhs) {
            return lhs * rhs;
        }

        template<bool PX, bool PY, bool PZ, typename... Args>
        void fixPosition(Vec3 &vec, Args... args);

        template<bool PX, bool PY, bool PZ>
        inline void fixPosition(Vec3 &vec, const double &dx, const double &dy, const double &dz) {
            if (PX) {
                vec[0] -= floor((vec[0] + .5 * dx) / dx) * dx;
            }
            if (PY) {
                vec[1] -= floor((vec[1] + .5 * dx) / dx) * dx;
            }
            if (PZ) {
                vec[2] -= floor((vec[2] + .5 * dx) / dx) * dx;
            }
        };

        void fixPosition(Vec3 &vec, const std::array<bool, 3> &periodic, const std::array<double, 3> &boxSize);

        template<bool PX, bool PY, bool PZ, typename... Args>
        Vec3 shortestDifference(const Vec3 &lhs, const Vec3 &rhs, Args... args);

        template<bool PX, bool PY, bool PZ>
        inline Vec3 shortestDifference(const Vec3 &lhs, const Vec3 &rhs, const double &dx, const double &dy, const double &dz) {
            auto dv = rhs - lhs;
            if (PX) {
                if (dv[0] > dx * .5) dv[0] -= dx;
                else if (dv[0] <= dx * .5) dv[0] += dx;
            }
            if (PY) {
                if (dv[1] > dy * .5) dv[1] -= dy;
                else if (dv[1] <= dy * .5) dv[1] += dy;
            }
            if (PZ) {
                if (dv[2] > dz * .5) dv[2] -= dz;
                else if (dv[2] <= dz * .5) dv[2] += dz;
            }
            return dv;
        };

        template<>
        inline Vec3 shortestDifference<false, false, false>(const Vec3 &lhs, const Vec3 &rhs) {
            return rhs - lhs;
        };

        template<bool PX, bool PY, bool PZ, typename... Args>
        double distSquared(const Vec3 &lhs, const Vec3 &rhs, Args... args);

        template<bool PX, bool PY, bool PZ>
        inline double distSquared(const Vec3 &lhs, const Vec3 &rhs, const double &dx, const double &dy, const double &dz) {
            auto dv = shortestDifference<PX, PY, PZ>(lhs, rhs, dx, dy, dz);
            return dv * dv;
        };

        template<>
        inline double distSquared<false, false, false>(const Vec3 &lhs, const Vec3 &rhs) {
            return (lhs - rhs) * (lhs - rhs);
        };

        double distSquared(const Vec3 &lhs, const Vec3 &rhs, const std::array<bool, 3> &periodic, const std::array<double, 3> &boxSize);
    }
}

#endif //READDY_MAIN_VEC3_H
