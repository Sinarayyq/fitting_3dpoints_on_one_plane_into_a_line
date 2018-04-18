#include <visualize.h>


boost::shared_ptr<pcl::visualization::PCLVisualizer> visCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_to_visualize, std::string window_label, camera_position camera_pos)
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer(window_label));
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addPointCloud(cloud_to_visualize, "sample cloud");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "sample cloud");
	//viewer->addCoordinateSystem(1, "global");
	viewer->addCoordinateSystem(500, 0, 0, 0, xy);

	// Set camera position according to different input.
	switch (camera_pos)
	{
	case xy:
		viewer->setCameraPosition(
			0, 0, 1,//double pos_x, double pos_y, double pos_z,
			0, 0, -1,//double view_x, double view_y, double view_z,
			1, 0, 0);//double up_x, double up_y, double up_z, int viewport = 0
	case yz:
		viewer->setCameraPosition(
			1, 0, 0,//double pos_x, double pos_y, double pos_z,                                    
			-1, 0, 0,//double view_x, double view_y, double view_z,
			0, 1, 0);//double up_x, double up_y, double up_z, int viewport = 0
	case xz:
		viewer->setCameraPosition(
			0, -1, 0,//double pos_x, double pos_y, double pos_z,                                    
			0, 1, 0,//double view_x, double view_y, double view_z,
			0, 0, 1);//double up_x, double up_y, double up_z, int viewport = 0
	}

	return (viewer);
}

void visualizePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, std::string viewer_window_label, camera_position camera_pos)
{
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;

	// Make use of the 
	viewer = visCloud(cloud, viewer_window_label, camera_pos);

	// Reset camera according to the input data. Zoom out so that all data points can be viewed.
	viewer->resetCamera();

	while (!viewer->wasStopped())
	{
		viewer->spinOnce();
	}
	viewer->close();

}



void visualizePointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2, std::string label_viewer_window, camera_position camera_pos)
{
	// Window setup
	pcl::visualization::PCLVisualizer viewer(label_viewer_window);
	viewer.setBackgroundColor(0, 0, 0);
	//viewer.addCoordinateSystem(1);
	viewer.addCoordinateSystem(500, 0, 0, 0, xy);
	viewer.setCameraPosition(0, 0, 1,//double pos_x, double pos_y, double pos_z,                                    
		0, 0, 0,//double view_x, double view_y, double view_z,
		0, 1, 0);//double up_x, double up_y, double up_z, int viewport = 0);

				 // Add the first cloud: size 2, white(default)
	viewer.addPointCloud(cloud1, "cloud1");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloud1");
	// Define color "red"
	pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> red(cloud2, 255, 255, 0);
	// Add the second cloud: size 3, red (defined)
	viewer.addPointCloud<pcl::PointXYZ>(cloud2, red, "cloud2");
	viewer.setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "cloud2");

	//Auto recenters the view.
	viewer.resetCamera();

	while (!viewer.wasStopped())
	{
		viewer.spinOnce();
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
	viewer.close();
}
