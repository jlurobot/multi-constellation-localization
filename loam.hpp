#ifndef __LOAM_H__
#define __LOAM_H__

#include <algorithm>
#include <vector>
#include <cmath>
#include "adaptor.hpp"
#include <atomic>
#include <chrono>
#include "common.h"

struct smoothness_t{ 
    float value;
    size_t ind;

    bool operator<(const smoothness_t& other) const {
        return value < other.value;
    }
};

struct ranged_points {
    size_t range_offsets[64];

    template<typename value_type>
    array_view<value_type> ring_span(int ring_id, value_type* points) {
        if(ring_id == 0) {
            return array_view(points, range_offsets[0]);
        }
        return array_view(points + range_offsets[ring_id - 1], range_offsets[ring_id] - range_offsets[ring_id - 1]);
    }

};

template<typename point_type>
struct cloud_featured {
    std::vector<point_type> corner_points;
    std::vector<point_type> surface_points;
};

template<typename point_type>
inline void downsample_surf(std::vector<point_type>& surface_points) {
    auto selected = surface_points.begin();
    for(auto i = selected + 1; i != surface_points.end(); i++) {
        if(p2(i->x - selected->x) + p2(i->y - selected->y) + p2(i->z - selected->z) < 4.0f) {
            continue;
        }
        *++selected = *i;
    }
    surface_points.erase(selected + 1, surface_points.end());
}

template<typename point_type>
inline cloud_featured<point_type> get_features(point_type* begin, point_type* end) {

    constexpr size_t H_SCAN = 1800;
    constexpr float edgeThreshold = 1.0;
    constexpr float surfThreshold = 0.1;

    cloud_featured<point_type> features;
    ranged_points points;

    end = std::remove_if(begin, end, [](const point_type& p) {
        return p2(p.x) + p2(p.y) + p2(p.z) < 2.0f;
    });

    std::vector<float> range(std::distance(begin, end));
    std::vector<int> cols(std::distance(begin, end));
    // since points is arranged by time, we can use std::stable_partition to split the points into rings

    point_type* ring_start = begin;
    size_t ring_id = 0;
    while(ring_start != end && ring_id < 64) {
        point_type* ring = std::stable_partition(ring_start, end, [ring_id](point_type& p) { return p.ring == ring_id; });
        points.range_offsets[ring_id] = ring - begin;
        ring_id ++;
        ring_start = ring;
    }

    // now we have the points arranged by rings, we can calculate the features 
    // for each ring
    std::vector<float> curvature(std::distance(begin, end));
    std::vector<smoothness_t> smoothness(std::distance(begin, end));

    for(point_type* i = begin; i != end; i++) {
        float d = i->x * i->x + i->y * i->y + i->z * i->z;
        range[i - begin] = std::sqrt(d);

        float angle = std::atan2(i->x, i->y) * 180.0f / M_PI;
        int columnIdn = -round((angle-90.0f)/(360.0f / H_SCAN)) + H_SCAN/2;
        if (columnIdn >= H_SCAN)
            columnIdn -= H_SCAN;
        
        if (columnIdn < 0 || columnIdn >= H_SCAN)
            columnIdn = 0;
        cols[i - begin] = columnIdn;

    }

    std::vector<bool> neighbor_picked(std::distance(begin, end), false);
    std::vector<bool> flag(std::distance(begin, end), false);

    size_t cloudSize = curvature.size();
    for (int i = 5; i < cloudSize - 5; i++)
    {
        float curv = range[i-5] + range[i-4]
                        + range[i-3] + range[i-2]
                        + range[i-1] - range[i] * 10
                        + range[i+1] + range[i+2]
                        + range[i+3] + range[i+4]
                        + range[i+5];            

        curvature[i] = curv*curv;
        smoothness[i].value = curvature[i];
        smoothness[i].ind = i;
    }

    // mark occluded points and parallel beam points
    for (size_t i = 5; i < cloudSize - 6; ++i)
    {
        // occluded points
        float depth1 = range[i];
        float depth2 = range[i+1];
        int diff = std::abs(cols[i+1] - cols[i]);

        if (diff < 10){
            // 10 pixel diff in range image
            if (depth1 - depth2 > 0.3){
                neighbor_picked[i - 5] = true;
                neighbor_picked[i - 4] = true;
                neighbor_picked[i - 3] = true;
                neighbor_picked[i - 2] = true;
                neighbor_picked[i - 1] = true;
                neighbor_picked[i] = true;
            }else if (depth2 - depth1 > 0.3){
                neighbor_picked[i + 1] = true;
                neighbor_picked[i + 2] = true;
                neighbor_picked[i + 3] = true;
                neighbor_picked[i + 4] = true;
                neighbor_picked[i + 5] = true;
                neighbor_picked[i + 6] = true;
            }
        }
        // parallel beam
        float diff1 = std::abs(range[i-1] - range[i]);
        float diff2 = std::abs(range[i+1] - range[i]);

        if (diff1 > 0.02 * range[i] && diff2 > 0.02 * range[i])
            neighbor_picked[i] = true;
    }

    for (int i = 0; i < ring_id; i++)
    {
        auto cloud_span = points.ring_span(i, begin);
        auto smoothness_span = points.ring_span(i, smoothness.data());

        for (int j = 0; j < 6; j++)
        {

            int sp = (cloud_span.size() * j) / 6;
            int ep = (cloud_span.size() * (j + 1)) / 6 - 1;

            if (sp >= ep)
                continue;

            std::sort(smoothness_span.begin() + sp, smoothness_span.begin() + ep);

            int largestPickedNum = 0;
            for (int k = ep; k >= sp; k--)
            {
                int ind = smoothness_span[k].ind;
                if (neighbor_picked[ind] == false && curvature[ind] > edgeThreshold && range[ind] > 2.0f)
                {
                    largestPickedNum++;
                    if (largestPickedNum <= 20){
                        flag[ind] = true;
                        features.corner_points.push_back(begin[ind]);
                    } else {
                        break;
                    }

                    neighbor_picked[ind] = true;
                    for (int l = 1; l <= 5; l++)
                    {
                        int columnDiff = std::abs(int(cols[ind + l] - cols[ind + l - 1]));
                        if (columnDiff > 10)
                            break;
                        neighbor_picked[ind + l] = true;
                    }
                    for (int l = -1; l >= -5; l--)
                    {
                        int columnDiff = std::abs(int(cols[ind + l] - cols[ind + l + 1]));
                        if (columnDiff > 10)
                            break;
                        neighbor_picked[ind + l] = true;
                    }
                }
            }

            for (int k = sp; k <= ep; k++)
            {
                int ind = smoothness_span[k].ind;
                if (neighbor_picked[ind] == false && curvature[ind] < surfThreshold && range[ind] > 2.0f)
                {

                    flag[ind] = false;
                    neighbor_picked[ind] = true;

                    for (int l = 1; l <= 5; l++) {

                        int columnDiff = std::abs(cols[ind + l] - cols[ind + l - 1]);
                        if (columnDiff > 10)
                            break;

                        neighbor_picked[ind + l] = true;
                    }
                    for (int l = -1; l >= -5; l--) {

                        int columnDiff = std::abs(cols[ind + l] - cols[ind + l + 1]);
                        if (columnDiff > 10)
                            break;

                        neighbor_picked[ind + l] = true;
                    }
                }
            }

            for (int k = sp; k <= ep; k++)
            {
                if (!flag[smoothness_span[k].ind]){
                    features.surface_points.push_back(cloud_span[k]);
                }
            }
        }

    }

    downsample_surf(features.surface_points);
    return features;
}


struct searched_line {
    float nx, ny, nz;
    float cx, cy, cz;
    bool ok;
};

template<typename point_type, typename query_point_type>
inline searched_line search_line(const array_adaptor<point_type>& tree, query_point_type point) {
    
    size_t pointSearchInd[5];
    float pointSearchSqDis[5];

    tree.query(point, 5, pointSearchInd, pointSearchSqDis);
    Eigen::Matrix3f A1 = Eigen::Matrix3f::Zero();

    if (pointSearchSqDis[4] < 1.0) {
        float cx = 0, cy = 0, cz = 0;
        for (int j = 0; j < 5; j++) {
            cx += tree.m_data[pointSearchInd[j]].x;
            cy += tree.m_data[pointSearchInd[j]].y;
            cz += tree.m_data[pointSearchInd[j]].z;
        }
        cx /= 5; cy /= 5;  cz /= 5;

        float a11 = 0, a12 = 0, a13 = 0, a22 = 0, a23 = 0, a33 = 0;
        for (int j = 0; j < 5; j++) {
            float ax = tree.m_data[pointSearchInd[j]].x - cx;
            float ay = tree.m_data[pointSearchInd[j]].y - cy;
            float az = tree.m_data[pointSearchInd[j]].z - cz;

            a11 += ax * ax; a12 += ax * ay; a13 += ax * az;
            a22 += ay * ay; a23 += ay * az;
            a33 += az * az;
        }
        a11 /= 5; a12 /= 5; a13 /= 5; a22 /= 5; a23 /= 5; a33 /= 5;

        A1(0, 0) = a11; A1(0, 1) = a12; A1(0, 2) = a13;
        A1(1, 0) = a12; A1(1, 1) = a22; A1(1, 2) = a23;
        A1(2, 0) = a13; A1(2, 1) = a23; A1(2, 2) = a33;

        Eigen::EigenSolver<Eigen::Matrix3f> es(A1);
        Eigen::Vector3f D1 = es.eigenvalues().real();
        Eigen::Matrix3f V1 = es.eigenvectors().real();

        if (D1(0) > 3 * D1(1)) {
            searched_line line;
            line.nx = V1(0, 0);
            line.ny = V1(1, 0);
            line.nz = V1(2, 0);
            line.cx = cx;
            line.cy = cy;
            line.cz = cz;
            line.ok = true;
            return line;
        }
    }
    searched_line line;
    line.ok = false;
    return line;
}

struct coeff {
    float px, py, pz;
    float x, y, z;
    float b;
    float s;
};


template<typename point_type>
inline coeff line_coeff(const searched_line& line, const point_type& p) {
    coeff c;
    if(!line.ok) {
        c.px = p.x;
        c.py = p.y;
        c.pz = p.z;
        c.x = 0;
        c.y = 0;
        c.z = 0;
        c.b = 0;
        c.s = 0;
        return c;
    }

    float x0 = p.x;
    float y0 = p.y;
    float z0 = p.z;

    float x1 = line.cx + 0.1 * line.nx;
    float y1 = line.cy + 0.1 * line.ny;
    float z1 = line.cz + 0.1 * line.nz;
    float x2 = line.cx - 0.1 * line.nx;
    float y2 = line.cy - 0.1 * line.ny;
    float z2 = line.cz - 0.1 * line.nz;

    float a012 = sqrt(((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1)) * ((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1)) 
                    + ((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1)) * ((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1)) 
                    + ((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1)) * ((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1)));

    float l12 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));

    float la = ((y1 - y2)*((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1)) 
                + (z1 - z2)*((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1))) / a012 / l12;

    float lb = -((x1 - x2)*((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1)) 
                - (z1 - z2)*((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1))) / a012 / l12;

    float lc = -((x1 - x2)*((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1)) 
                + (y1 - y2)*((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1))) / a012 / l12;

    float ld2 = a012 / l12;

    float s = 1 - 0.9 * fabs(ld2);

    c.x = s * la;
    c.y = s * lb;
    c.z = s * lc;
    c.b = s * ld2;
    c.s = s;

    c.px = p.x;
    c.py = p.y;
    c.pz = p.z;
    
    return c;
}

struct plane {
    float a, b, c, d;
    bool ok;
};

template<typename point_type>
inline plane search_plane(const array_adaptor<point_type>& tree, const point_type& p)
{
    size_t pointSearchInd[5];
    float pointSearchSqDis[5];

    tree.query(p, 5, pointSearchInd, pointSearchSqDis);

    Eigen::Matrix<float, 5, 3> matA0;
    Eigen::Matrix<float, 5, 1> matB0;
    Eigen::Vector3f matX0;

    matA0.setZero();
    matB0.fill(-1);
    matX0.setZero();

    if (pointSearchSqDis[4] < 1.0) {
        for (int j = 0; j < 5; j++) {
            matA0(j, 0) = tree.m_data[pointSearchInd[j]].x;
            matA0(j, 1) = tree.m_data[pointSearchInd[j]].y;
            matA0(j, 2) = tree.m_data[pointSearchInd[j]].z;
        }

        matX0 = matA0.colPivHouseholderQr().solve(matB0);

        float pa = matX0(0, 0);
        float pb = matX0(1, 0);
        float pc = matX0(2, 0);
        float pd = 1;

        float ps = sqrt(pa * pa + pb * pb + pc * pc);
        pa /= ps; pb /= ps; pc /= ps; pd /= ps;

        bool planeValid = true;
        for (int j = 0; j < 5; j++) {
            if (fabs(pa * tree.m_data[pointSearchInd[j]].x +
                        pb * tree.m_data[pointSearchInd[j]].y +
                        pc * tree.m_data[pointSearchInd[j]].z + pd) > 0.2) {
                planeValid = false;
                break;
            }
        }

        plane pl;
        pl.a = pa;
        pl.b = pb;
        pl.c = pc;
        pl.d = pd;
        
        pl.ok = planeValid;
        return pl;
    }

    plane pl;
    pl.ok = false;
    return pl;
}

//Calculating the point-to-plane distance, referring to Eq.(8) 
template<typename point_type>
inline coeff plane_coeff(const plane& pl, const point_type& p) {
    if (pl.ok) {
        float pd2 = pl.a * p.x + pl.b * p.y + pl.c * p.z + pl.d;

        float s = 1 - 0.9 * fabs(pd2) / sqrt(sqrt(p.x * p.x + p.y * p.y + p.z * p.z));
        coeff c;

        c.x = s * pl.a;
        c.y = s * pl.b;
        c.z = s * pl.c;
        c.b = s * pd2;  
        c.s = s;
        c.px = p.x;
        c.py = p.y;
        c.pz = p.z;

        return c;
    }

    coeff c;
    c.s = 0;
    return c;
}

template<typename point_type>
inline void transform_point(point_type& p, const Eigen::Matrix4d& t) {
    Eigen::Vector4d v(p.x, p.y, p.z, 1);
    v = t * v;
    p.x = v.x();
    p.y = v.y();
    p.z = v.z();
}

struct jacobian {
    double j[6];
};

struct jacobian_g {
    double srx, crx, sry, cry, srz, crz;
};

inline void init_jacobian_g(jacobian_g& g, const Transform& t) {
    g.srx = sin(t.roll);
    g.crx = cos(t.roll);
    g.sry = sin(t.pitch);
    g.cry = cos(t.pitch);
    g.srz = sin(t.yaw);
    g.crz = cos(t.yaw);
}

//Compute the Jacobian matrix
inline jacobian J(const coeff& c, const jacobian_g& g) {
    jacobian j;
    j.j[3] = (c.py * (g.srx*g.srz + g.sry*g.crx*g.crz) + c.pz *(g.srz*g.crx - g.srx*g.sry*g.crz)) * c.x + 
            (c.py*(g.sry*g.srx*g.crx - g.srx*g.crz) - c.pz*(g.srx*g.sry*g.srz + g.crx * g.crz)) * c.y +
            (c.py * g.crx*g.cry - c.pz * g.srx*g.cry) * c.z;
    
    j.j[4] = (-c.px * g.sry * g.crz + c.py * g.srx * g.cry * g.crz + c.pz * g.crx * g.cry * g.crz) * c.x +
            (-c.px*g.sry*g.srz + c.py*g.srx*g.srz*g.cry + c.pz*g.srz*g.crx*g.cry) * c.y +
             (-c.px*g.cry - c.py*g.srx*g.sry - c.pz*g.sry*g.crx) * c.z;
             
    j.j[5] = (-c.px*g.srz*g.cry + c.py*(-g.sry*g.srx*g.srz - g.crx*g.crz) + c.pz*(-g.sry*g.srz*g.crx + g.srx*g.crz))* c.x +
             (c.px*g.cry*g.crz + c.py*(g.sry*g.srx*g.crz - g.srz*g.crx) + c.pz*(g.sry*g.crx*g.crz + g.srx*g.srz))*c.y;

    j.j[0] = c.x;
    j.j[1] = c.y;
    j.j[2] = c.z;
    return j;
}


// Gauss-Newton method for optimization, referring to Eq.(10)
template<typename source_point_type, typename target_point_type>
inline Transform __LM_iteration(const cloud_featured<source_point_type>& source, array_adaptor<target_point_type>& corner,
    array_adaptor<target_point_type>& surf, Eigen::MatrixXf& A, Eigen::VectorXf& b, const Transform& initial_guess, float* loss = nullptr) {
    
    Eigen::Matrix4d transform = to_eigen(initial_guess);
    jacobian_g g;
    init_jacobian_g(g, initial_guess);


    std::atomic<int> index = 0;

    size_t corner_size = source.corner_points.size();
    size_t surf_size = source.surface_points.size();
    if(loss) {
        *loss = 0;
    }

#pragma omp parallel for
    for(size_t i = 0; i < corner_size + surf_size; ++i) {
        coeff c;
        c.s = 0;
        if(i < corner_size) {
            source_point_type p2 = source.corner_points[i];
            transform_point(p2, transform);
        
            searched_line sl = search_line(corner, p2);
            if(sl.ok) {
                c = line_coeff(sl, p2);
            }
        } else {
            int idx = i - corner_size;
            source_point_type p2 = source.surface_points[idx];
            transform_point(p2, transform);
        
            plane sp = search_plane(surf, p2);
            if(sp.ok) {
                c = plane_coeff(sp, p2);
            }
        }

        if(c.s < 0.1f)
            continue;
        jacobian j = J(c, g);
        int idx = index.fetch_add(1);
        A(idx, 0) = j.j[0];
        A(idx, 1) = j.j[1];
        A(idx, 2) = j.j[2];
        A(idx, 3) = j.j[3];
        A(idx, 4) = j.j[4];
        A(idx, 5) = j.j[5];
        b(idx) = -c.b;

        if(loss) {
            loss[0] += c.b;
        }
    }

    int idx = index.load();
    Eigen::Matrix<float, 6, 6> ATA = A.topRows(idx).transpose() * A.topRows(idx);
    Eigen::Matrix<float, 6, 1> ATb = A.topRows(idx).transpose() * b.topRows(idx);
    Eigen::Matrix<float, 6, 1> x = ATA.householderQr().solve(ATb);

    Transform delta;
    delta.x = initial_guess.x + x(0, 0);
    delta.y = initial_guess.y + x(1, 0);
    delta.z = initial_guess.z + x(2, 0);
    delta.roll = initial_guess.roll + x(3, 0);
    delta.pitch = initial_guess.pitch + x(4, 0);
    delta.yaw = initial_guess.yaw + x(5, 0);

    if(loss) {
        loss[0] /= index;
    }
    return delta;
}

template<typename _Call>
struct Caller {
    _Call call;
    Caller(_Call&& call) : call(std::forward<_Call>(call)) {}
    template<typename... Args>
    auto operator()(Args&&... args) {
        return call(std::forward<Args>(args)...);
    }
};

template<>
struct Caller<nullptr_t> {
    Caller(nullptr_t) {}
    template<typename... Args>
    void operator()(Args&&... args) {}
};

template<typename source_point_type, typename target_point_type, typename _Call = nullptr_t>
Transform LM(const cloud_featured<source_point_type>& source, const cloud_featured<target_point_type>& target, 
            const Transform& initial_guess = Transform(), float* loss = nullptr, _Call&& call = nullptr) {
    array_adaptor<target_point_type> corner(target.corner_points.data(), target.corner_points.size());
    array_adaptor<target_point_type> surf(target.surface_points.data(), target.surface_points.size());
    Eigen::MatrixXf A(source.corner_points.size() + source.surface_points.size(), 6);
    Eigen::VectorXf b(source.corner_points.size() + source.surface_points.size());

    Transform result = initial_guess;
    Caller<_Call> caller(std::forward<_Call>(call));
    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < 30; ++iter) {
        Transform u = __LM_iteration(source, corner, surf, A, b, result, loss);
        float deltaR = sqrtf(p2(u.roll - result.roll) + p2(u.pitch - result.pitch) + p2(u.yaw - result.yaw));
        float deltaT = sqrtf(p2(u.x - result.x) + p2(u.y - result.y) + p2(u.z - result.z));
        // printf("iter: %d, deltaR: %f, deltaT: %f\r\n", iter, deltaR, deltaT);
        result = u;
        if (deltaR < 0.0005 && deltaT < 0.0005) {
            break;
        }
        caller(result);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    return result;
}

#endif