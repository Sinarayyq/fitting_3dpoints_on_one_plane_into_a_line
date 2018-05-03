#pragma once 

#ifndef VISUALIZER_H
#define VISUALIZER_H

#include <pcl/visualization/cloud_viewer.h>
#include <pcl/point_types.h>

enum camera_position { xy, yz, xz };

void visualizePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2, std::string label_viewer_window, camera_position camera_pos);

void visualizePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, std::string viewer_window_label, camera_position camera_pos);




#endif