#ifndef __DATA_SOURCE_H__
#define __DATA_SOURCE_H__

#include <pcl/point_cloud.h>
#include <string>
#include "common.h"

using Reader = bool(*)(const char* path, pcl::PointCloud<XYZIRT>& cloud);

struct DataSourceStruct {
    Reader reader;
    std::vector<std::string> paths;
    std::vector<double> timestamps;
    std::vector<Eigen::Matrix4d> traces;
};

enum class DataSourceType {
    PCD_Directory,
    BIN_Directory,
};

typedef struct DataSourceStruct* DataSource;

DataSource createDataSource(DataSourceType type, const char* path, const char* timestamp, const char* traces);

void destroyDataSource(DataSource source);

#endif // __DATA_SOURCE_H__
