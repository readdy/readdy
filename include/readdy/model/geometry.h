//
// Created by mho on 8/23/21.
//

#pragma once

#include "readdy/common/common.h"

namespace readdy::model::geometry {

template<typename dtype=readdy::scalar>
struct Sphere {

    static constexpr std::string_view name = "Sphere";

    template<bool inclusion>
    [[nodiscard]] bool contains(Vec3 position) const {
        if constexpr(inclusion) {
            return (position - center).normSquared() < radius * radius;
        } else {
            return (position - center).normSquared() > radius * radius;
        }
    }

    template<bool inclusion>
    [[nodiscard]] Vec3 smallestDifference(Vec3 position) const {
        auto delta = position - center;
        auto distToOriginSquared = delta.normSquared();
        if constexpr(inclusion) {
            if (distToOriginSquared < radius * radius) {
                return {0, 0, 0};  // zero vector
            } else {
                auto distToOrigin = std::sqrt(distToOriginSquared);
                return delta * (distToOrigin - radius) / distToOrigin;
            }
        } else {
            if (distToOriginSquared > radius * radius) {
                return {0, 0, 0};  // zero vector
            } else {
                auto distToOrigin = std::sqrt(distToOriginSquared);
                return delta * (distToOrigin - radius) / distToOrigin;
            }
        }
    }

    [[nodiscard]] std::string describe() const {
        return fmt::format("center={}, radius={}", center, radius);
    }

    Vec3 center {};
    dtype radius {1.};
};

template<typename dtype=scalar>
struct Box {

    static constexpr std::string_view name = "Box";

    template<bool inclusion>
    [[nodiscard]] bool contains(Vec3 position) const {
        bool result = true;
        #pragma unroll
        for(std::uint8_t d = 0; d < 3; ++d) {
            if constexpr(inclusion) {
                result &= position[d] > v0[d] && position[d] < v1[d];
            } else {
                result &= position[d] < v0[d] && position[d] > v1[d];
            }
        }
        return result;
    }

    template<bool inclusion>
    [[nodiscard]] Vec3 smallestDifference(Vec3 position) const {
        Vec3 difference {0, 0, 0};

        #pragma unroll
        for(std::uint8_t d = 0; d < 3; ++d) {
            if constexpr(inclusion) {
                if (position[d] < v0[d]) {
                    difference[d] = position[d] - v0[d];
                } else if (position[d] > v1[d]) {
                    difference[d] = position[d] - v1[d];
                }
            } else {
                if (position[d] > v0[d]) {
                    difference[d] = position[d] - v0[d];
                } else if (position[d] < v1[d]) {
                    difference[d] = position[d] - v1[d];
                }
            }
        }

        return difference;
    }

    [[nodiscard]] std::string describe() const {
        return fmt::format("minimum vertex v0={} and maximum vertex v1={}", v0, v1);
    }

    // lower left vertex and upper right vertex
    Vec3 v0 {}, v1 {};
};

template<typename dtype=scalar>
struct Capsule {

    static constexpr std::string_view name = "Capsule";

    template<typename T>
    auto differentSign(T t1, T t2) const {
        return (t1 >= 0 && t2 < 0) || (t1 < 0 && t2 >= 0);
    }

    [[nodiscard]] Sphere<dtype> closestCircleCenter(const Vec3 &position) const {
        auto normalizedDirection = direction / direction.norm();
        // we define the line x(l) = x_0 + l*v, where v is the normalized direction vector
        // then for given position the lambda is (x_0 - p)*v*v:
        auto v = (((position - center) * normalizedDirection) * normalizedDirection);
        auto lambda = v.norm();
        // check if any of the signs in v differ from direction:
        if (differentSign(v.x, normalizedDirection.x) || differentSign(v.y, normalizedDirection.y) || differentSign(v.z, normalizedDirection.z)) {
            lambda = -lambda;
        }
        // clamp lambda to half length of capsule
        lambda = std::clamp(lambda, -length / 2, length / 2);
        // now we obtain point corresponding to lambda
        auto circleCenter = center + lambda * normalizedDirection;
        return {.center = circleCenter, .radius = radius};
    }

    template<bool inclusion>
    [[nodiscard]] bool contains(Vec3 position) const {
        return closestCircleCenter(position).template contains<inclusion>(position);
    }

    template<bool inclusion>
    [[nodiscard]] Vec3 smallestDifference(Vec3 position) const {
        return closestCircleCenter(position).template smallestDifference<inclusion>(position);
    }

    [[nodiscard]] std::string describe() const {
        return fmt::format("center={}, direction={}, radius={}, length={}", center, direction, radius, length);
    }

    Vec3 center {}, direction {};
    dtype radius {1.}, length {1.};
};

}
