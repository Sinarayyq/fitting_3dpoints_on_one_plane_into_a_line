////iterative_ransac
//#include <iostream>
//#include <vector>
//#include <time.h>
//
//#include <pcl/common/geometry.h>
//#include <pcl/filters/extract_indices.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/io/io.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/io/vtk_lib_io.h>
//#include <pcl/ModelCoefficients.h>
//#include <pcl/point_cloud.h>
//#include <pcl/point_types.h>
//#include <pcl/ros/conversions.h>
//#include <pcl/sample_consensus/method_types.h>
//#include <pcl/sample_consensus/model_types.h>
//#include <pcl/segmentation/sac_segmentation.h>
//#include <pcl/filters/voxel_grid.h>
//#include <pcl/filters/extract_indices.h>
//
//
//#include <io.h>
//#include <visualize.h>
//
//using namespace std;
//
//class PatchType
//{
//public:
//	bool plane;      //0 = not a plane, 1 = is a plane
//	bool cylinder;
//	bool cone;
//
//	pcl::ModelCoefficients::Ptr coefficients_plane;
//	pcl::ModelCoefficients::Ptr coefficients_cylinder;
//	pcl::ModelCoefficients::Ptr coefficients_cone;
//
//	PatchType();
//	~PatchType() {};
//};
//
//PatchType::PatchType():plane(0),cylinder(0),cone(0), coefficients_plane(new pcl::ModelCoefficients()), 
//                       coefficients_cylinder(new pcl::ModelCoefficients()), coefficients_cone(new pcl::ModelCoefficients())
//{
//	//coefficients_plane = NULL;
//	//coefficients_cylinder = NULL;
//	//coefficients_cone = NULL;
//}
//
//
//pcl::PointCloud<pcl::PointXYZ>::Ptr BuildPointCloudStructure(std::vector<pcl::PointXYZ> subsampled_points, int j, int i, std::vector<int> number)
//{
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	int size_cloud = i - j + 1;
//
//	cloud->width = size_cloud;
//	cloud->height = 1;
//	cloud->is_dense = false;
//	cloud->points.resize(size_cloud);
//	for (size_t k = 0; k < size_cloud; ++k)
//	{
//		cloud->points[k].x = subsampled_points[number[k+j]].x;
//		cloud->points[k].y = subsampled_points[number[k+j]].y;
//		cloud->points[k].z = subsampled_points[number[k+j]].z;
//	}
//
//	return cloud;
//	
//}
//
//pcl::PointCloud<pcl::Normal>::Ptr BuildPointNormalStructure(std::vector<pcl::Normal> subsample_points_normals, int j, int i, std::vector<int> number)
//{
//	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
//	int size_cloud = i - j + 1;
//
//	normals->width = size_cloud;
//	normals->height = 1;
//	normals->is_dense = false;
//	normals->points.resize(size_cloud);
//	for (size_t k = 0; k < size_cloud; ++k)
//	{
//		normals->points[k].normal_x = subsample_points_normals[number[k+j]].normal_x;
//		normals->points[k].normal_y = subsample_points_normals[number[k+j]].normal_y;
//		normals->points[k].normal_z = subsample_points_normals[number[k+j]].normal_z;
//	}
//
//	return normals;
//
//}
//
//
//void setSegmentationParametersForPlane(pcl::SACSegmentation<pcl::PointXYZ>& seg_plane)
//{
//	// Create the segmentation object
//	//pcl::SACSegmentation<pcl::PointXYZ> seg;
//	// Optional
//	seg_plane.setOptimizeCoefficients(true);
//	// Mandatory
//	seg_plane.setModelType(pcl::SACMODEL_PLANE);
//	seg_plane.setMethodType(PLANE_METHOD_TYPE);
//	seg_plane.setMaxIterations(PLANE_MAX_NUM_ITER);
//	seg_plane.setDistanceThreshold(PLANE_TOL);
//}
//
//void setSegmentationParametersForCylinder(pcl::SACSegmentationFromNormals<pcl::PointXYZ, pcl::Normal>& seg_cyl)
//{
//	// Set all the parameters for cylinder segmentation object
//	seg_cyl.setOptimizeCoefficients(true);
//	seg_cyl.setModelType(pcl::SACMODEL_CYLINDER);
//	seg_cyl.setMethodType(CYL_METHOD_TYPE);
//	seg_cyl.setNormalDistanceWeight(CYL_WEIGHT_NORMAL_DISTANCE);
//	seg_cyl.setMaxIterations(CYL_MAX_NUM_ITER);
//	seg_cyl.setDistanceThreshold(CYL_TOL);
//	//double a, b;
//	//seg_cyl.getRadiusLimits(a, b);
//	//std::cout << "a = " << a << "  b =" << b << std::endl;
//	seg_cyl.setRadiusLimits(CYL_MIN_RADIUS_LIMIT, CYL_MAX_RADIUS_LIMIT);
//}
//
//void setSegmentationParametersForCone(pcl::SACSegmentationFromNormals<pcl::PointXYZ, pcl::Normal>& seg_cone)
//{
//	// Set all the parameters for cone segmentation object
//	seg_cone.setOptimizeCoefficients(true);
//	seg_cone.setModelType(pcl::SACMODEL_CONE);
//	seg_cone.setMethodType(CONE_METHOD_TYPE);
//	seg_cone.setMinMaxOpeningAngle(CONE_MIN_OPENING_ANGLE / 180.0  * M_PI, CONE_MAX_OPENING_ANGLE / 180.0 * M_PI); //it is in radiants; min=5degree, max=80degree
//	seg_cone.setNormalDistanceWeight(CONE_WEIGHT_NORMAL_DISTANCE);
//	seg_cone.setMaxIterations(CONE_MAX_NUM_ITER);
//	seg_cone.setDistanceThreshold(CONE_TOL);
//	//double a, b;
//	//seg_cone.getRadiusLimits(a, b);
//	//std::cout << "a = " << a << "  b =" << b << std::endl;
//	seg_cone.setRadiusLimits(CONE_MIN_RADIUS_LIMIT, CONE_MAX_RADIUS_LIMIT);
//}
//
//
//bool IteraticeRansac(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, pcl::PointCloud<pcl::Normal>::Ptr normals, PatchType *patchtype_temp)
//{
//	size_t size_cloud = cloud->size();
//	bool flag = 0;
//
//	//plane
//	pcl::SACSegmentation<pcl::PointXYZ> seg_plane;
//	pcl::PointIndices::Ptr inliers_plane(new pcl::PointIndices());
//	pcl::ModelCoefficients::Ptr coefficients_plane(new pcl::ModelCoefficients());
//	setSegmentationParametersForPlane(seg_plane);
//	seg_plane.setInputCloud(cloud);
//	seg_plane.segment(*inliers_plane, *coefficients_plane);
//	if (((inliers_plane->indices.size()) / size_cloud) > PERCENTAGE)
//	{
//		flag = 1;
//		(*patchtype_temp).plane = 1;
//		*((*patchtype_temp).coefficients_plane) = *coefficients_plane;
//	}
//
//		
//
//	//cylinder
//	pcl::SACSegmentationFromNormals<pcl::PointXYZ, pcl::Normal> seg_cylinder;
//	pcl::PointIndices::Ptr inliers_cylinder(new pcl::PointIndices());
//	pcl::ModelCoefficients::Ptr coefficients_cylinder(new pcl::ModelCoefficients());
//	setSegmentationParametersForCylinder(seg_cylinder);
//	seg_cylinder.setInputCloud(cloud);
//	seg_cylinder.setInputNormals(normals);
//	seg_cylinder.segment(*inliers_cylinder, *coefficients_cylinder);
//	if (((inliers_cylinder->indices.size()) / size_cloud) > PERCENTAGE)
//	{
//		flag = 1;
//		(*patchtype_temp).cylinder = 1;
//		*((*patchtype_temp).coefficients_cylinder) = *coefficients_cylinder;
//	}
//
//			
//
//	//cone
//	pcl::SACSegmentationFromNormals<pcl::PointXYZ, pcl::Normal> seg_cone;
//	pcl::PointIndices::Ptr inliers_cone(new pcl::PointIndices());
//	pcl::ModelCoefficients::Ptr coefficients_cone(new pcl::ModelCoefficients());
//	setSegmentationParametersForCone(seg_cone);
//	seg_cone.setInputCloud(cloud);
//	seg_cone.setInputNormals(normals);
//	seg_cone.segment(*inliers_cone, *coefficients_cone);
//	if (((inliers_cone->indices.size()) / size_cloud) > PERCENTAGE)
//	{
//		flag = 1;
//		(*patchtype_temp).cone = 1;
//		*((*patchtype_temp).coefficients_cone) = *coefficients_cone;
//	}
//
//		
//	
//	
//
//	return flag;
//}
//
//
//
//
//int main(int argc, char** argv)
//{
//	readParameterFile("C:\\Development\\iterative_ransac\\source\\input extract indices.txt");
//	std::vector<pcl::PointXYZ> subsampled_points = ReadTXTFileAsPointCloud("C:\\Users\\sinara\\Desktop\\subsampled_points.txt");
//	std::vector<pcl::Normal> subsample_points_normals = ReadTXTFilePointNormal("C:\\Users\\sinara\\Desktop\\subsample_points_normals.txt");
//	
//
//
//	time_t start = clock();
//	int size_subsampled_points = subsampled_points.size();
//	
//	if (size_subsampled_points < 1)
//	{
//		std::cout << "No subsampled points" << std::endl;
//		return -1;
//	}
//
//	int j = 0; // start
//	std::vector<int> number = ReadTXTNumber("C:\\Users\\sinara\\Desktop\\distance\\0.txt");
//	pcl::PointCloud<pcl::PointXYZ>::Ptr subsampled_point_cloud = BuildPointCloudStructure(subsampled_points, 0, size_subsampled_points - 1, number);
//	std::vector<PatchType> patchtype;
//	PatchType patchtype_last;
//	for (int i = 40; i < size_subsampled_points; i++)
//	{
//		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = BuildPointCloudStructure(subsampled_points, j, i, number);            //need to be changed when adding generatrix 
//		pcl::PointCloud<pcl::Normal>::Ptr normals = BuildPointNormalStructure(subsample_points_normals, j, i, number);    //need to be changed when adding generatrix
//		PatchType patchtype_temp;
//		if (IteraticeRansac(cloud, normals, &patchtype_temp))   //1 = points can be recognized as at least one type; 0 = points cannot be recognized as any type
//		{
//			patchtype_last = patchtype_temp;
//			if (i == size_subsampled_points - 1)
//			{
//				visualizePointCloud(subsampled_point_cloud, cloud, "patch", xy);
//			}
//			continue;
//		}
//		else
//		{
//			patchtype.push_back(patchtype_last);
//			visualizePointCloud(subsampled_point_cloud, cloud, "patch", xy);
//			j = i;
//			//number = ReadTXTNumber("C:\\Users\\sinara\\Desktop\\distance\\"+std::to_string(i)+".txt");
//			//subsampled_point_cloud = BuildPointCloudStructure(subsampled_points, 0, size_subsampled_points - 1, number);
//			i += 40;
//		}
//		
//	}
//	patchtype.push_back(patchtype_last);
//	time_t end = clock();
//	std::cout << (double)(end - start) / CLOCKS_PER_SEC * 1000.0 << std::endl;
//
//
//	std::system("pause");
//	return 0;
//}
//
//
//




////pcl subsampling
//#include <iostream>
//#include <visualize.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/point_types.h>
//#include <pcl/filters/voxel_grid.h>
//
//#include <pcl/features/normal_3d.h>
//#include <pcl/features/pfh.h>
//#include <pcl/features/fpfh.h>
//#include <pcl/features/fpfh_omp.h>
//#include <pcl/features/vfh.h>
//#include <pcl/features/gfpfh.h>
//
//#include <pcl/point_cloud.h>
//#include <pcl/kdtree/kdtree_flann.h>
//#include <vector>
//
//#include <pcl/sample_consensus/ransac.h>
//#include <pcl/sample_consensus/sac_model_plane.h>
//#include <pcl/sample_consensus/sac_model_sphere.h>
//#include <pcl/filters/random_sample.h>
//#include <pcl/filters/uniform_sampling.h>
//#include <pcl/filters/normal_space.h>
//
//#include <pcl/kdtree/kdtree_flann.h>
//
//int main(int argc, char** argv)
//{
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
//
//	// Fill in the cloud data
//	pcl::PCDReader reader;
//	// Replace the path below with the path where you saved your file
//	std::string load_file;
//	std::getline(std::cin, load_file);
//	reader.read(load_file, *cloud); // Remember to download the file first!
//	//std::cout << cloud->points[0].x << " " << cloud->points[0].y << " " << cloud->points[0].z << std::endl;
//	/*std::cerr << "PointCloud before filtering: " << cloud->width * cloud->height
//		<< " data points (" << pcl::getFieldsList(*cloud) << ").";*/
//
//	//pcl::UniformSampling<pcl::PointXYZ> random_sample;
//	//random_sample.setInputCloud(cloud);      //the size of it 6500 
//	//random_sample.setRadiusSearch(120); //5000 
//	////random_sample.setSeed(rand());
//	//random_sample.filter(*cloud_filtered);
//	//random_sample.
//	// Create the filtering object
//	//pcl::<pcl::PointXYZ> sor;
//	/*sor.setInputCloud(cloud);
//	sor.setLeafSize(120, 120, 120);
//	sor.filter(*cloud_filtered);*/
//	pcl::VoxelGrid<pcl::PointXYZ> sample;
//	sample.setInputCloud(cloud);
//	sample.setLeafSize(120, 120, 120);
//	sample.filter(*cloud_filtered);
//	//visualizePointCloud(cloud, cloud_filtered, "patch", xy);
//
//	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
//	kdtree.setInputCloud(cloud);
//
//	int K = 1;
//
//	size_t size_cloud_filtered = cloud_filtered->size();
//	std::vector<int> index(K);
//	std::vector<double> d(K);
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_final(new pcl::PointCloud<pcl::PointXYZ>);
//	for (size_t i = 0; i < size_cloud_filtered; i++)
//	{
//		if (kdtree.nearestKSearch(cloud_filtered->points[i], K, index, d) > 0)
//		{
//			size_t no = rand();
//			cloud_final->push_back(pcl::PointXYZ(cloud->points[index[K-1]].x, cloud->points[index[K-1]].y, cloud->points[index[K-1]].z));
//		}
//	}
//
//
//
//
//
//
//
//	visualizePointCloud(cloud, cloud_final, "patch", xy);
//	/*std::cerr << "PointCloud after filtering: " << cloud_filtered->width * cloud_filtered->height
//		<< " data points (" << pcl::getFieldsList(*cloud_filtered) << ").";*/
//
//	pcl::io::savePCDFileASCII("C:\\Users\\sinara\\Desktop\\grid_test_output.pcd", *cloud_final);
//	/*pcl::PCDWriter writer;
//	writer.write("C:\\Users\\sinara\\Desktop\\table_scene_lms400_downsampled.pcd", *cloud_filtered,
//		Eigen::Vector4f::Zero(), Eigen::Quaternionf::Identity(), false);*/
//
//	return (0);
//}






#include <iostream>
#include <vector>
#include <time.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <math.h>
#include <pcl/io/pcd_io.h>
#include <pcl/PointIndices.h>

#include <io.h>
#include <visualize.h>


Eigen::Matrix3f buildTransformMatrixFromAxis(Eigen::Vector3f x_axis, Eigen::Vector3f y_axis, Eigen::Vector3f z_axis)
{
	using namespace Eigen;

	Matrix3f transform_matrix;
	transform_matrix << x_axis, y_axis, z_axis;
	return transform_matrix;
}

Eigen::Matrix3f getTransformMatrixForAlignmentWithNormalToPlane(pcl::PointXYZ origin, pcl::PointXYZ point_on_x_axis, std::vector<double> plane_param)
{
	Eigen::Vector3f x_axis;
	Eigen::Vector3f y_axis;
	Eigen::Vector3f z_axis;

	//z axis definition	
	z_axis[0] = plane_param[0] ;
	z_axis[1] = plane_param[1] ;
	z_axis[2] = plane_param[2] ;
	double norm_z = z_axis.norm();
	for (std::size_t i = 0; i < 3; ++i)
	{
		z_axis[i] = z_axis[i] / norm_z;
	}

	//x axis
	x_axis[0] = point_on_x_axis.x - origin.x;
	x_axis[1] = point_on_x_axis.y - origin.y;
	x_axis[2] = point_on_x_axis.z - origin.z;
	//std::cout << x_axis << std::endl;
	double norm_x = x_axis.norm();
	for (std::size_t i = 0; i < 3; ++i)
	{
		x_axis[i] = x_axis[i] / norm_x;
	}

	//y axis definition
	y_axis = z_axis.cross(x_axis);
	double norm_y = y_axis.norm();
	for (std::size_t i = 0; i < 3; ++i)
	{
		y_axis[i] = y_axis[i] / norm_y;
	}
	/*std::cerr << "x_axis = (" << x_axis[0] << "," << x_axis[1] << "," << x_axis[2] << ")" << std::endl;
	std::cerr << "y_axis = (" << y_axis[0] << "," << y_axis[1] << "," << y_axis[2] << ")" << std::endl;
	std::cerr << "z_axis = (" << z_axis[0] << "," << z_axis[1] << "," << z_axis[2] << ")" << std::endl;*/

	//Get the transformation matrix referred to the new coordinate system

	Eigen::Matrix3f trasf_matrix = buildTransformMatrixFromAxis(x_axis, y_axis, z_axis);
	Eigen::Matrix3f trasf_matrix_inv = trasf_matrix.inverse();

	return trasf_matrix_inv;
	//return trasf_matrix;
}


pcl::PointCloud<pcl::PointXYZ>::Ptr transformCloudByMatrix(pcl::PointCloud<pcl::PointXYZ>::Ptr patch_cloud, Eigen::Matrix3f transform_matrix, pcl::PointXYZ origin)
{
	const int num_nodes = patch_cloud->points.size();
	Eigen::Vector3f point, transformed_point;
	pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_points(new pcl::PointCloud<pcl::PointXYZ>);
	Eigen::Vector3f appex(origin.x, origin.y, origin.z);

	for (std::size_t i = 0; i < num_nodes; ++i)
	{
		point = Eigen::Vector3f(patch_cloud->at(i).getArray3fMap());
		/*std::cout << point << std::endl;
		std::cout << appex << std::endl;
		std::cout << transform_matrix << std::endl;
		std::cout << point - appex << std::endl;*/
		transformed_point = transform_matrix*(point-appex);
		
		//std::cout << transformed_point << std::endl;
		(*transformed_points).push_back(pcl::PointXYZ(transformed_point[0], transformed_point[1], transformed_point[2]));
	}

	return transformed_points;
}

int LinearRegression(pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_points, double &a, double &b)
{
	double xsum = 0, x2sum = 0, ysum = 0, xysum = 0;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
	int number = (*transformed_points).size();
	for (int i = 0; i < number; i++)
	{
		xsum = xsum + (*transformed_points)[i].x;                        //calculate sigma(xi)
		ysum = ysum + (*transformed_points)[i].y;                        //calculate sigma(yi)
		x2sum = x2sum + pow((*transformed_points)[i].x, 2);                //calculate sigma(x^2i)
		xysum = xysum + (*transformed_points)[i].x * (*transformed_points)[i].y;                    //calculate sigma(xi*yi)
	}
	a = (number * xysum - xsum * ysum) / (number * x2sum - xsum * xsum);            //calculate slope
	b = (x2sum * ysum - xsum * xysum) / (x2sum * number - xsum * xsum);            //calculate intercept

	return 0;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr ProjectPointsOntoLine(pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_points, double a, double b, double &error)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr projected_points(new pcl::PointCloud<pcl::PointXYZ>);
	int size = (*transformed_points).size();
	double sum = 0;
	for (int i = 0; i < size; i++)
	{
		double x0 = (*transformed_points)[i].x, y0 = (*transformed_points)[i].y;
		double x = (x0 + a*y0 - a*b) / (1 + a*a);
		double y = a*x + b;
		(*projected_points).push_back(pcl::PointXYZ(x, y, 0));
		sum += sqrt(pow(x - x0, 2) + pow(y - y0, 2));
	}
	error = sum;
	return projected_points;
}

int FindHeadAndTail(pcl::PointCloud<pcl::PointXYZ>::Ptr projected_points, pcl::PointCloud<pcl::PointXYZ>::Ptr &transformed_head_tail)
{
	int size = (*projected_points).size();
	double x0 = (*projected_points)[0].x, y0 = (*projected_points)[0].y;
	double x = (*projected_points)[1].x - x0;
	int head_no = 0, tail_no = 0;
	double head_dist = 0, tail_dist = 0;

	for (int i = 1; i < size; i++)
	{
		double xn = (*projected_points)[i].x, yn = (*projected_points)[i].y;
		double dist = sqrt(pow(xn - x0, 2) + pow(yn - y0, 2));;
		if ((xn - x0)*x >= 0)
		{
			if (dist > tail_dist)
			{
				tail_dist = dist;
				tail_no = i;
			}
		}
		else
		{
			if (dist > head_dist)
			{
				head_dist = dist;
				head_no = i;
			}
		}
	}
	(*transformed_head_tail).push_back((*projected_points)[head_no]);
	(*transformed_head_tail).push_back((*projected_points)[tail_no]);
	return 0;
}

int Fitting3DPointsToLine(pcl::PointCloud<pcl::PointXYZ>::Ptr &points_list, std::vector<double> plane_normal, double &error, pcl::PointXYZ &head, pcl::PointXYZ &tail, bool flag)
{
	visualizePointCloud(points_list, "points_list", xy);
	Eigen::Matrix3f transform_matrix1 = getTransformMatrixForAlignmentWithNormalToPlane((*points_list).points[0], (*points_list).points[1], plane_normal);
	pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_points = transformCloudByMatrix(points_list, transform_matrix1, (*points_list).points[0]);
	//pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_points = transformPlanarPatchPoints(points_list, plane_normal);
	
	//visualizePointCloud(transformed_points, "transformed_points", xy);
	/*visualizePointCloud(points_list, transformed_points, "transformed_points", xy);
	std::cout << sqrt(pow((*points_list)[0].x - (*points_list)[1].x,2) + pow((*points_list)[0].y - (*points_list)[1].y,2) +pow((*points_list)[0].z - (*points_list)[1].z,2)) << std::endl;
	std::cout << sqrt(pow((*transformed_points)[0].x - (*transformed_points)[1].x, 2) + pow((*transformed_points)[0].y - (*transformed_points)[1].y, 2) + pow((*transformed_points)[0].z - (*transformed_points)[1].z,2)) << std::endl;*/
	double a, b;
	LinearRegression(transformed_points, a, b);
	//pcl::PointCloud<pcl::PointXYZ>::Ptr line(new pcl::PointCloud<pcl::PointXYZ>);
	//(*line).push_back(pcl::PointXYZ(0, b, 0));
	//(*line).push_back(pcl::PointXYZ(1000, 1000*a+b, 0));

	pcl::PointCloud<pcl::PointXYZ>::Ptr projected_points = ProjectPointsOntoLine(transformed_points, a, b, error);
	pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_head_tail(new pcl::PointCloud<pcl::PointXYZ>);
	if (flag == 0)
	{
		(*transformed_head_tail).push_back(projected_points->points[0]);
		(*transformed_head_tail).push_back(projected_points->points[1]);
	}
	else
	{
		FindHeadAndTail(projected_points, transformed_head_tail);
	}
	
	//visualizePointCloud(transformed_points, transformed_head_tail, "transformed_head_tail", xy);




	pcl::PointCloud<pcl::PointXYZ>::Ptr coordinate_system((new pcl::PointCloud<pcl::PointXYZ>));
	(*coordinate_system).push_back(pcl::PointXYZ(0, 0, 0));
	(*coordinate_system).push_back(pcl::PointXYZ(1, 0, 0));   //x
	(*coordinate_system).push_back(pcl::PointXYZ(0, 0, 1));   //z
	pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_coordinate_system = transformCloudByMatrix(coordinate_system, transform_matrix1, (*points_list).points[0]);
	//visualizePointCloud(points_list, coordinate_system, "coordinate_system", xy);
	//visualizePointCloud(transformed_points, transformed_coordinate_system, "transformed_coordinate_system", xy);
	//std::vector<double> transformed_coordinate_normal;
	std::vector<double> normal;
	normal.push_back((*transformed_coordinate_system)[2].x - (*transformed_coordinate_system)[0].x);
	normal.push_back((*transformed_coordinate_system)[2].y - (*transformed_coordinate_system)[0].y);
	normal.push_back((*transformed_coordinate_system)[2].z - (*transformed_coordinate_system)[0].z);





	Eigen::Matrix3f transform_matrix2 = getTransformMatrixForAlignmentWithNormalToPlane((*transformed_coordinate_system)[0], (*transformed_coordinate_system)[1], normal);
	pcl::PointCloud<pcl::PointXYZ>::Ptr transformed_transformed_points = transformCloudByMatrix(transformed_points, transform_matrix2, (*transformed_coordinate_system)[0]);
	visualizePointCloud(points_list, transformed_transformed_points, "transformed_transformed_points", xy);
	pcl::PointCloud<pcl::PointXYZ>::Ptr head_tail = transformCloudByMatrix(transformed_head_tail, transform_matrix2, (*transformed_coordinate_system)[0]);

	visualizePointCloud(points_list, head_tail, "head_tail", xy);
	head = (*head_tail)[0];
	tail = (*head_tail)[1];


	return 0;
}

int main()
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr points_list(new pcl::PointCloud<pcl::PointXYZ>);
	std::vector<double> plane_normal;
	std::string load_file;
	std::getline(std::cin, load_file);
	ReadTXTFileAsPointCloud(load_file.c_str(), points_list, plane_normal);
	double error;
	pcl::PointXYZ head, tail;
	bool flag = 0;

	Fitting3DPointsToLine(points_list, plane_normal, error, head, tail, flag);

	return 0;
}