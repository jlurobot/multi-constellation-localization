#ifndef __COMMON_H__
#define __COMMON_H__

#include <Eigen/Dense>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <chrono>

struct XYZIRT {
    PCL_ADD_POINT4D;
    PCL_ADD_INTENSITY;
    uint16_t ring;
    double time;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} EIGEN_ALIGN16;

POINT_CLOUD_REGISTER_POINT_STRUCT(XYZIRT, (float, x, x)(float, y, y)(float, z, z)(float, intensity, intensity)(std::uint16_t, ring, ring)(double, time, time))


template<typename T>
struct array_view {
    T* __ptr;
    size_t __size;

    array_view(T* begin, size_t size) : __ptr(begin), __size(size) {}

    T& operator[](size_t i) {
        return __ptr[i];
    }

    const T& operator[](size_t i) const {
        return __ptr[i];
    }

    T* begin() {
        return __ptr;
    }

    T* end() {
        return __ptr + __size;
    }

    const T* begin() const {
        return __ptr;
    }

    const T* end() const {
        return __ptr + __size;
    }

    size_t size() const {
        return __size;
    }

    bool empty() const {
        return __size == 0;
    }

    T* data() {
        return __ptr;
    }
};

template<typename P1, typename P2>
auto distance2(const P1& p1, const P2& p2) {
    return p2(p1.x - p2.x) + p2(p1.y - p2.y) + p2(p1.z - p2.z);
}

struct Transform {
    double x = 0.0f ,y = 0.0f ,z = 0.0f;
    double roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
};


static inline Eigen::Matrix4d to_eigen(const Transform& tr) {
    // T = s.Matrix([[1, 0, 0, tx], [0, 1, 0, ty], [0, 0, 1, tz], [0, 0, 0, 1]])
    // Rx = s.Matrix([[1, 0, 0, 0], [0, s.cos(rx), -s.sin(rx), 0], [0, s.sin(rx), s.cos(rx), 0], [0, 0, 0, 1]])
    // Ry = s.Matrix([[s.cos(ry), 0, s.sin(ry), 0], [0, 1, 0, 0], [-s.sin(ry), 0, s.cos(ry), 0], [0, 0, 0, 1]])
    // Rz = s.Matrix([[s.cos(rz), -s.sin(rz), 0, 0], [s.sin(rz), s.cos(rz), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

    // Transformation matrix
    // M = T * Rz * Ry * Rx
    
    Eigen::Matrix4d m = Eigen::Matrix4d::Identity();
    m(0, 3) = tr.x;
    m(1, 3) = tr.y;
    m(2, 3) = tr.z;

    double cr = cos(tr.roll);
    double sr = sin(tr.roll);
    double cp = cos(tr.pitch);
    double sp = sin(tr.pitch);
    double cy = cos(tr.yaw);
    double sy = sin(tr.yaw);

    m(0, 0) = cp * cy;
    m(0, 1) = cy*sp*sr - cr*sy;
    m(0, 2) = sr*sy + cr*cy*sp;
    m(1, 0) = cp*sy;
    m(1, 1) = cr*cy + sp*sr*sy;
    m(1, 2) = cr*sp*sy - cy*sr;
    m(2, 0) = -sp;
    m(2, 1) = cp*sr;
    m(2, 2) = cr*cp;

    return m;
}

static inline Transform from_eigen(const Eigen::Matrix4d& m) {
    Transform tr;
    tr.x = m(0, 3);
    tr.y = m(1, 3);
    tr.z = m(2, 3);
    tr.pitch = atan2(-m(2, 0), sqrt(m(0, 0)*m(0, 0) + m(1, 0)*m(1, 0)));

    double c = cos(tr.pitch);
    tr.yaw = atan2(m(1, 0) / c, m(0, 0) / c);
    tr.roll = atan2(m(2, 1) / c, m(2, 2) / c);

    return tr;
}


std::vector<std::string> list_files(const char* path, const char* ext);

struct progress_bar {
    int previous_percent = 0;
    size_t total = 0;
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    
    progress_bar(size_t total);
    void print_progress(size_t value);
    void done();

private:
    void do_print(size_t markers, size_t spaces, size_t value); 
};

#endif