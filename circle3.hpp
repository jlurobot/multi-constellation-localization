
#ifndef __CIRCLE3_HPP__
#define __CIRCLE3_HPP__

#include <iostream>
#include <assert.h>
#include <eigen3/Eigen/Dense>
#include <cmath>


struct Point3
{
    double x;
    double y;
    double z;
};

struct Circle3 : Point3
{
    double radius;
};

template <typename T>
static inline auto p2(T x) { return x * x; }

bool __check_circle3(Circle3 c1, Circle3 c2)
{
    double d = sqrt(p2(c1.x - c2.x) + p2(c1.y - c2.y) + p2(c1.z - c2.z));
    return d < c1.radius + c2.radius && d > fabsf(c1.radius - c2.radius);
}
bool __check_circle3(Circle3 c1, Circle3 c2, Circle3 c3)
{
    return __check_circle3(c1, c2) && __check_circle3(c1, c3) && __check_circle3(c2, c3);
}

bool __point_on_circle3(Circle3 c, Point3 p)
{
    double d = sqrt(p2(c.x - p.x) + p2(c.y - p.y) + p2(c.z - p.z));
    return fabs(d - c.radius) < 0.000001;
}

static inline bool __double_nequals_zero(double A)
{
    return A != 0.0;
}

#include <tuple>

bool __solver_check_normalize(Point3 n0)
{
    return __double_nequals_zero(n0.x) || __double_nequals_zero(n0.y) || __double_nequals_zero(n0.z);
}

std::tuple<double, double, bool> __solve_nhle2(double A1, double A2, double B1, double B2, double D1, double D2)
{
    double R = A1 * B2 - A2 * B1;
    if (__double_nequals_zero(R))
        return {
            (D2 * B1 - D1 * B2) / R,
            (D1 * A2 - D2 * A1) / R,
            true};

    return {0.0, 0.0, false};
}

std::tuple<double, double, double, bool> __solve_nhle3(double A1, double A2, double B1, double B2, double C1, double C2, double D1, double D2, Point3 n0)
{
    if (__double_nequals_zero(n0.z))
    {
        auto [R1, R2, ok] = __solve_nhle2(A1, A2, B1, B2, D1, D2);
        return {R1, R2, 0.0, ok};
    }

    if (__double_nequals_zero(n0.y))
    {
        auto [R1, R2, ok] = __solve_nhle2(A1, A2, C1, C2, D1, D2);
        return {R1, 0.0, R2, ok};
    }

    auto [R1, R2, ok] = __solve_nhle2(B1, B2, C1, C2, D1, D2);
    return {0.0, R1, R2, ok};
}

//The implementation of Eq.(15) in the paper
inline void __solve_circle3_delte(Circle3 circle, Point3 x, Point3* delta) {
    double d = sqrt(p2(circle.x - x.x) + p2(circle.y - x.y) + p2(circle.z - x.z));
    double delta_d = d - circle.radius;
    delta->x -= delta_d * (x.x - circle.x) / d;
    delta->y -= delta_d * (x.y - circle.y) / d;
    delta->z -= delta_d * (x.z - circle.z) / d;
}

std::tuple<Point3, bool> __solve_circle3_unsolve(Circle3 c1, Circle3 c2, Circle3 c3) {
    Point3 inital = {
        (c1.x + c2.x + c3.x) / 3.0f,
        (c1.y + c2.y + c3.y) / 3.0f,
        (c1.z + c2.z + c3.z) / 3.0f
    };

    for(int i = 0 ; i < 100; i++) {
        Point3 delta_dir = {0.0f, 0.0f, 0.0f};
        __solve_circle3_delte(c1, inital, &delta_dir);
        __solve_circle3_delte(c2, inital, &delta_dir);
        __solve_circle3_delte(c3, inital, &delta_dir);

        if (fabs(delta_dir.x) < 0.00001 && fabs(delta_dir.y) < 0.00001 && fabs(delta_dir.z) < 0.00001) {
            return { inital, true };
        }

        if(isnan(delta_dir.x) || isnan(delta_dir.y) || isnan(delta_dir.z)) {
            return { inital, false };
        }

        inital.x += delta_dir.x * 0.3f;
        inital.y += delta_dir.y * 0.3f;
        inital.z += delta_dir.z * 0.3f;
    }
    return { inital, true };
}

//When there are three equations (three intersecting spheres), the coordinates of the intersection point are solved.
std::tuple<Point3, Point3, bool> solve_circle3(Circle3 c1, Circle3 c2, Circle3 c3)
{
    constexpr std::tuple<Point3, Point3, bool> __invalid_result = {Point3{0.0f, 0.0f, 0.0f}, Point3{0.0f, 0.0f, 0.0f}, false};
    Point3 n0 = {
        (c1.y - c3.y) * (c2.z - c3.z) - (c1.z - c3.z) * (c2.y - c3.y),
        (c1.z - c3.z) * (c2.x - c3.x) - (c1.x - c3.x) * (c2.z - c3.z),
        (c1.x - c3.x) * (c2.y - c3.y) - (c1.y - c3.y) * (c2.x - c3.x)};

    if (!__solver_check_normalize(n0))
    {
        return __invalid_result;
    }

    Point3 G = {
        p2(c1.x) + p2(c1.y) + p2(c1.z) - p2(c1.radius),
        p2(c2.x) + p2(c2.y) + p2(c2.z) - p2(c2.radius),
        p2(c3.x) + p2(c3.y) + p2(c3.z) - p2(c3.radius)};

    double A1 = 2.0f * (c1.x - c2.x), A2 = 2.0f * (c1.x - c3.x);
    double B1 = 2.0f * (c1.y - c2.y), B2 = 2.0f * (c1.y - c3.y);
    double C1 = 2.0f * (c1.z - c2.z), C2 = 2.0f * (c1.z - c3.z);
    double D1 = G.y - G.x, D2 = G.z - G.x;

    auto [R1, R2, R3, ok] = __solve_nhle3(A1, A2, B1, B2, C1, C2, D1, D2, n0);
    
    if (!ok)
    {
        return __invalid_result;
    } 

    double a = p2(n0.x) + p2(n0.y) + p2(n0.z);
    double b = 2.0f * (n0.x * (R1 - c1.x) + n0.y * (R2 - c1.y) + n0.z * (R3 - c1.z));
    double c = p2(R1 - c1.x) + p2(R2 - c1.y) + p2(R3 - c1.z) - p2(c1.radius);

    double delta = p2(b) - 4.0f * a * c;

    if (delta < 0.0f)
    {
        auto [p, ok] = __solve_circle3_unsolve(c1, c2, c3);
        return {p, p, ok};
    }

    double t1 = (-b + sqrt(delta)) / (2.0f * a);
    double t2 = (-b - sqrt(delta)) / (2.0f * a);

    Point3 p1 = {
        R1 + n0.x * t1,
        R2 + n0.y * t1,
        n0.z * t1};

    Point3 p2 = {
        R1 + n0.x * t2,
        R2 + n0.y * t2,
        n0.z * t2};

    return {p1, p2, true};
}

double __solve_distance(Point3 __p1, Point3 __p2)
{
    return sqrt(p2(__p1.x - __p2.x) + p2(__p1.y - __p2.y) + p2(__p1.z - __p2.z));
}

Point3 __middle(Point3 __p1, Point3 __p2) {
    return {
        (__p1.x + __p2.x) / 2.0,
        (__p1.y + __p2.y) / 2.0,
        (__p1.z + __p2.z) / 2.0,
    };
}


//Solving for the initial values s0, the output serves as the input for Eq.(14), (15), and (16).
/*
Due to space constraints in the article, this implementation was presented in a simplified manner. 
However, in practice, as long as the initial values are not at the sphere centers, 
their influence on the final results is negligible.
*/
std::pair<Point3, bool> __solve_initialize(const Circle3 *c, size_t point_count) {
    auto [P1, P2, ok] = solve_circle3(c[0], c[1], c[2]);
    int i = 1;
    for(; i < int(point_count) - 4; i++) {
        std::tie(P1, P2, ok) = solve_circle3(c[i], c[i + 1], c[i + 2]);
        if(ok)
            break;
    }

    if(!ok) {
        return {Point3{0.0f, 0.0f, 0.0f}, false};
    }

    auto [P3, P4, ok2] = solve_circle3(c[point_count - 1], c[point_count - 2], c[point_count - 3]);

    for(int j = 1; point_count - 1 - j > i; j++) {
        std::tie(P3, P4, ok2) = solve_circle3(c[point_count - 1 - j], c[point_count - 2 - j], c[point_count - 3 - j]);\
        if(ok2)
            break;
    }   

    if(!ok2) {
        return {Point3{0.0f, 0.0f, 0.0f}, false};
    }

    auto R13 = __solve_distance(P1, P3);
    auto R14 = __solve_distance(P1, P4);


    auto R23 = __solve_distance(P2, P3);
    auto R24 = __solve_distance(P2, P4);

    auto M = std::min(std::min(R13, R14), std::min(R23, R24));

    if(M == R13) return {__middle(P1, P3), true};
    if(M == R14) return {__middle(P1, P4), true};
    if(M == R23) return {__middle(P2, P3), true};
    if(M == R24) return {__middle(P2, P4), true};

    return {Point3{0.0f, 0.0f, 0.0f}, false};
}

//Iteratively solve the equation system, corresponding to the Eq.(14), (15), and (16) in the paper.
std::pair<Point3, bool> solve_circles(const Circle3 *observed_circles, size_t point_count, Point3* initial = nullptr)
{
    assert(point_count > 3);

    auto [__initial_point, ok] = __solve_initialize(observed_circles, point_count);
    if (!ok)
    {
        return {__initial_point, false};
    }

    if(initial) {
        initial[0] = __initial_point;
    }

    Point3 __current_point = __initial_point;

    constexpr double rate = 0.1f;
    constexpr size_t iterations = 30;
    constexpr double esp = 0.0005;

    double prev_loss = 100000.0f;

    for (int i = 0; i < iterations; i++)
    {
        double loss = 0.0f;
        Point3 delta = {0.0f, 0.0f, 0.0f};
        for (int j = 0; j < point_count; j++)
        {
            __solve_circle3_delte(observed_circles[j], __initial_point, &delta);
            float d = sqrtf(p2(observed_circles[j].x - __initial_point.x) + p2(observed_circles[j].y - __initial_point.y) + p2(observed_circles[j].z - __initial_point.z));
            loss += p2(d - observed_circles[j].radius);
        }
        
        if(fabsf(loss - prev_loss) < esp && i > 5) {
            break;
        }
        
        float dynamic_rate = 0.9f / point_count;
        __initial_point.x += dynamic_rate * delta.x;
        __initial_point.y += dynamic_rate * delta.y;
        __initial_point.z += dynamic_rate * delta.z;
        prev_loss = loss;
    }

    if(isnan(__initial_point.x) || isnan(__initial_point.y) || isnan(__initial_point.z)) {
        return {__current_point, true};
    }

    return { __initial_point, true};
}

#endif // CIRCLE3_H
