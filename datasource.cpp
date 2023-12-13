#include "datasource.h"

#include <filesystem>
#include <pcl/io/pcd_io.h>

#include "common.h"
struct trace_item {
    double time;
    float x;
    float y;
    float z;
    float roll, pitch, yaw;
};

//Interpolation is applied to the GPS data, ensuring that each point cloud is associated with a corresponding GPS location.
bool gpslerp(const std::vector<trace_item>& traces, double t, trace_item* result) {
    if(traces.size() == 0 || t < traces[0].time || t > traces.back().time) {
        return false;
    }

    auto it = std::lower_bound(traces.begin(), traces.end(), t, [](const trace_item& trace, double t) {
        return trace.time < t;
    });

    if(it == traces.begin()) {
        *result = *it;
        return false;
    }

    auto prev = it - 1;

    double t0 = prev->time;
    double t1 = it->time;

    double w0 = (t1 - t) / (t1 - t0);

    result->time = t;
    result->x = prev->x * w0 + it->x * (1.0 - w0);
    result->y = prev->y * w0 + it->y * (1.0 - w0);
    result->z = prev->z * w0 + it->z * (1.0 - w0);

    result->roll = prev->roll * w0 + it->roll * (1.0 - w0);
    result->pitch = prev->pitch * w0 + it->pitch * (1.0 - w0);
    result->yaw = prev->yaw * w0 + it->yaw * (1.0 - w0);

    return true;
}

#include "files.h"

//Reading timestamps from self-collected data.
std::vector<double> read_timestamp(const char* file) {
    FILE* fp = fopen(file, "r");
    if(!fp) {
        printf("Failed to open %s(%s)\n", file, strerror(errno));
        return {};
    }
    int frame_id= 0;
    double timestamp;
    std::vector<double> result;
    while(fscanf(fp, "%d %lf", &frame_id, &timestamp) == 2) {

        result.push_back(timestamp / 1.e9);
    }

    fclose(fp);
    return result;
}

//Reading trajectory from self-collected dataset.
void load_trace9(const char* filename, std::vector<trace_item>& items) {
    double cov_x, cov_y, cov_z;
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("Failed to open %s(%s)\n", filename, strerror(errno));
        return;
    }

    char buf[1024];
    while (fgets(buf, sizeof(buf), fp)) {
        trace_item item;
        if (sscanf(buf, "%lf,%f,%f,%f,%f,%f,%f,%lf,%lf,%lf", &item.time, &item.x, &item.y, &item.z, &item.roll, &item.pitch, &item.yaw, &cov_x, &cov_y, &cov_z) != 10) {
            printf("Failed to parse %s\n", buf);
            continue;
        }
        items.push_back(item);
    }

    fclose(fp);
}

//Searching for the position corresponding to the current timestamp among all timestamps and trajectories.
std::vector<trace_item> downsample_traces(std::vector<trace_item>& input, const std::vector<double>& timestamp){
    std::vector<trace_item> result;
    for(auto val : timestamp) {
        trace_item item;
        if(gpslerp(input, val, &item)) {
            result.push_back(item);
        } else {
            if(val < input[0].time) {
                result.push_back(input[0]);
            } else {
                result.push_back(input.back());
            }
        }
    }

    return result;
}

static bool pcd_loader(const char* path, pcl::PointCloud<XYZIRT>& cloud) {
    return pcl::io::loadPCDFile(path, cloud) >= 0;
}

DataSource createPCDSource(const char* path, const char* timestamp, const char* traces) {
    DataSource source = new DataSourceStruct();
    source->paths = list_files(path, ".pcd");
    source->timestamps = read_timestamp(timestamp);
    printf("Found %zd pcd files\n", source->paths.size());
    printf("Found %zd timestamps\n", source->timestamps.size());
    
    std::vector<trace_item> trace;
    load_trace9(traces, trace);
    

    auto strace = downsample_traces(trace, source->timestamps);

    for(auto& var : strace) {
        Transform t;
        t.x = var.x;
        t.y = var.y;
        t.z = var.z;
        t.roll = var.roll;
        t.pitch = var.pitch;
        t.yaw = var.yaw;

        source->traces.push_back(to_eigen(t));
    }
    
    printf("Found %zd traces\n", trace.size());

    source->reader = &pcd_loader;
    return source;
}

struct origin_info {
    double ecef_x;
    double ecef_y;
    double ecef_z;
};


//Converting GPS coordinates to ECEF coordinate system.
inline void wgs84_to_ecef(double latitude, double longitude, double altitude, double* x, double* y, double* z) {
    constexpr double a = 6378137.0;
    constexpr double b = 6356752.3142;
    constexpr double e2 = (a * a - b * b) / (a * a);

    latitude = latitude / 180.0 * M_PI;
    longitude = longitude / 180.0 * M_PI;

    double N = a / sqrt(1 - e2 * sin(latitude) * sin(latitude));
    *x = (N + altitude) * cos(latitude) * cos(longitude);
    *y = (N + altitude) * cos(latitude) * sin(longitude);
    *z = (N * (1 - e2) + altitude) * sin(latitude);
}

//Converting GPS coordinates to a Cartesian coordinate system.
inline void wgs84_to_local_tangent(const origin_info* o, double latitude, double longitude, double altitude, double* local_x, double* local_y, double* local_z) {
    double ecef_x, ecef_y, ecef_z;
    wgs84_to_ecef(latitude, longitude, altitude, &ecef_x, &ecef_y, &ecef_z);

    double dx = o->ecef_x - ecef_x;
    double dy = o->ecef_y - ecef_y;
    double dz = o->ecef_z - ecef_z;

    double sin_lat = sin(latitude / 180.0 * M_PI);
    double cos_lat = cos(latitude / 180.0 * M_PI);
    double sin_lon = sin(longitude / 180.0 * M_PI);
    double cos_lon = cos(longitude / 180.0 * M_PI);

    *local_x = -sin_lon * dx + cos_lon * dy;
    *local_y = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
    *local_z = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;
}

//Reading trajectory from KITTI dataset.
static std::vector<Eigen::Matrix4d> load_traceX(const char* folder) {
    std::vector<std::string> files = list_files(folder, ".txt");
    std::vector<Eigen::Matrix4d> result;
    bool first_frame = true;
    origin_info o;
    
    for(auto& file : files) {
        FILE* fp = fopen(file.c_str(), "r");
        if(!fp) {
            printf("Failed to open %s(%s)\n", file.c_str(), strerror(errno));
            continue;
        }

        char buf[1024];
        Eigen::Matrix4d mat = Eigen::Matrix4d::Identity();

        while(fgets(buf, sizeof(buf), fp)) {
            double lat, lon, alt;
            Transform t;
            if(sscanf(buf, "%lf %lf %lf %lf %lf %lf", &lat, &lon, &alt, &t.roll, &t.pitch, &t.yaw) != 6) {
                printf("Failed to parse %s\n", buf);
                continue;
            }

            if(first_frame) {
                wgs84_to_ecef(lat, lon, alt, &o.ecef_x, &o.ecef_y, &o.ecef_z);
                first_frame = false;
            }
            wgs84_to_local_tangent(&o, lat, lon, alt, &t.x, &t.y, &t.z);
            result.push_back(to_eigen(t));
        }
        fclose(fp);
    }

    return result;
};
 
struct XYZI {
    float x;
    float y;
    float z;
    float intensity;
};

inline size_t get_rings(const XYZI* points, size_t count) {
        bool neg = false;
        for (size_t i = 0; i < count; i++) {
                float yaw = atan2f(points[i].y, points[i].x);
                if (yaw < 0)
                        neg = true;
                if (yaw > 0 && neg)
                        return i;
        }
        return count;
}

//Reading point cloud data from KITTI dataset.
static bool BIN_reader(const char* filename, pcl::PointCloud<XYZIRT>& cloud) {
    FILE* fp = fopen(filename, "rb");
    if(!fp) {
        printf("Failed to open %s(%s)\n", filename, strerror(errno));
        return false;
    }

    fseek(fp, 0, SEEK_END);
    size_t size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    size_t count = size / sizeof(XYZI);
    cloud.resize(count);

    std::vector<XYZI> buffer(count);

    if(fread(buffer.data(), sizeof(XYZI), count, fp) != count) {
        printf("Failed to read %s\n", filename);
        fclose(fp);
        return false;
    }

    fclose(fp);

    int ring_index = 0;
    size_t offset = 0;
    for(int i = 0 ;i < 64; i++) {
        size_t next_ring = get_rings(buffer.data() + offset, count - offset);
        for(size_t j = offset; j < offset + next_ring; j++) {
            cloud[j].x = buffer[j].x;
            cloud[j].y = buffer[j].y;
            cloud[j].z = buffer[j].z;
            cloud[j].intensity = buffer[j].intensity;
            cloud[j].ring = ring_index;
            cloud[j].time = 0.99 / count * j;
        }
        offset += next_ring;
        ring_index++;
    }

    if(offset != count) {
        printf("Failed to parse %s\n", filename);
        return false;
    }
    
    return true;
}

//Reading timestamp from KITTI dataset.
static double parse_timestamp(const char* line) {
    struct tm tm;
    memset(&tm, 0, sizeof(tm));
    int nano;
    if(sscanf(line, "%d-%d-%d %d:%d:%d.%d", &tm.tm_year, &tm.tm_mon, &tm.tm_mday, &tm.tm_hour, &tm.tm_min, &tm.tm_sec, &nano) != 7) {
        printf("Failed to parse %s\n", line);
        return 0;
    }

    tm.tm_year -= 1900;
    tm.tm_mon -= 1;

    return (double)mktime(&tm) + nano / 1e9;
}

//Reading the timestamps corresponding to the KITTI point clouds.
static std::vector<double> load_timestamp_Kitti(const char* filename) {
    std::vector<double> result;

    FILE* fp = fopen(filename, "r");
    if(!fp) {
        printf("Failed to open %s(%s)\n", filename, strerror(errno));
        return result;
    }

    char buf[1024];
    while(fgets(buf, sizeof(buf), fp)) {
        double t = parse_timestamp(buf);
        if(t == 0) {
            printf("Failed to parse %s\n", buf);
            continue;
        }

        result.push_back(t);
    }

    fclose(fp);

    return result;
}

DataSource createBINSource(const char* path, const char* timestamp, const char* traces) {
    DataSource source = new DataSourceStruct();
    source->timestamps = load_timestamp_Kitti(timestamp);
    source->traces = load_traceX(traces);
    source->paths = list_files(path, ".bin");
    source->reader = &BIN_reader;
    return source;
}

DataSource createDataSource(DataSourceType type, const char* path, const char* timestamp, const char* traces)
{
    DataSource result = nullptr;
    switch(type) {
        case DataSourceType::BIN_Directory:
            result = createBINSource(path, timestamp, traces);
            break;
        case DataSourceType::PCD_Directory:
            result = createPCDSource(path, timestamp, traces);
            break;
        default:
            printf("Unknown data source type\n");
            return nullptr;
    }

    if(result == nullptr) {
        return nullptr;
    }

    if(result->paths.size() != result->timestamps.size()) {
        printf("Timestamps and files count mismatch\n");
        destroyDataSource(result);
        return nullptr;
    }

    if(result->traces.size() != result->timestamps.size()) {
        printf("Timestamps and traces count mismatch\n");
        destroyDataSource(result);
        return nullptr;
    }

    return result;
}

void destroyDataSource(DataSource source) {
    if(source == nullptr) {
        return;
    }

    delete source;
}
