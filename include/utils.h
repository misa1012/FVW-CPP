#ifndef FVW_UTILS_H
#define FVW_UTILS_H

#include <vector>
#include <cmath>

namespace fvw
{

    struct Vec3
    {
        double x, y, z;
        Vec3(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
        Vec3 operator+(const Vec3 &other) const { return {x + other.x, y + other.y, z + other.z}; }
        Vec3 operator-(const Vec3 &other) const { return {x - other.x, y - other.y, z - other.z}; }
        Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
        double norm() const { return std::sqrt(x * x + y * y + z * z); }
        double dot(const Vec3 &other) const { return x * other.x + y * other.y + z * other.z; }
        Vec3 cross(const Vec3 &other) const
        {
            return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
        }
    };

} // namespace fvw

#endif // FVW_UTILS_H