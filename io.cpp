#include "io.h"

//std::string PATH_PCD_DOC;
int PLANE_METHOD_TYPE;
int CYL_METHOD_TYPE;
int CONE_METHOD_TYPE;
double PLANE_TOL;
double CYL_TOL;
double CONE_TOL;
double MIN_INLIERS;
double MST_RATIO;
double MST_RATIO2;
int PLANE_MAX_NUM_ITER;
int CYL_MAX_NUM_ITER;
int CONE_MAX_NUM_ITER;
int K_FOR_NORMAL_SEARCH;
double CYL_WEIGHT_NORMAL_DISTANCE;
double CONE_WEIGHT_NORMAL_DISTANCE;
double CYL_MIN_RADIUS_LIMIT;
double CYL_MAX_RADIUS_LIMIT;
double CONE_MIN_RADIUS_LIMIT;
double CONE_MAX_RADIUS_LIMIT;
double CONE_MIN_OPENING_ANGLE;
double CONE_MAX_OPENING_ANGLE;
double TOLERANCE_FOR_ADJACENT_CHAINS;
double TOLERANCE_FOR_APPEX;
double POISSON_DISK_SAMPLING_RADIUS;
int POISSON_DISK_SAMPLING_NUMBER;
std::string PATH_TEMPORARY_FILE;
double PLANE_AREA_PERCENTAGE;
double CYL_AREA_PERCENTAGE;
double CONE_AREA_PERCENTAGE;
double PERCENTAGE;

const int INF = std::numeric_limits<int>::max();

void setParameters(std::vector<std::string>& parameters)
{
	for (std::size_t i = 0; i < parameters.size(); ++i)
	{
		switch (i)
		{
			/*case 0: {PATH_PCD_DOC = parameters[i]; }
			break;*/
		case 0: {
			if (std::stoi(parameters[i]) == 0)
			{
				PLANE_METHOD_TYPE = pcl::SAC_RANSAC;
			}
			else if (std::stoi(parameters[i]) == 1)
			{
				PLANE_METHOD_TYPE = pcl::SAC_PROSAC;
			}
			else if (std::stoi(parameters[i]) == 2)
			{
				PLANE_METHOD_TYPE = pcl::SAC_LMEDS;
			}
			else if (std::stoi(parameters[i]) == 3)
			{
				PLANE_METHOD_TYPE = pcl::SAC_MSAC;
			}
			else if (std::stoi(parameters[i]) == 4)
			{
				PLANE_METHOD_TYPE = pcl::SAC_RRANSAC;
			}
			else if (std::stoi(parameters[i]) == 5)
			{
				PLANE_METHOD_TYPE = pcl::SAC_RMSAC;
			}
			else if (std::stoi(parameters[i]) == 6)
			{
				PLANE_METHOD_TYPE = pcl::SAC_MLESAC;
			}
			else
			{
				PLANE_METHOD_TYPE = INF;
			}
		}
				break;
		case 1: {
			if (std::stoi(parameters[i]) == 0)
			{
				CYL_METHOD_TYPE = pcl::SAC_RANSAC;
			}
			else if (std::stoi(parameters[i]) == 1)
			{
				CYL_METHOD_TYPE = pcl::SAC_PROSAC;
			}
			else if (std::stoi(parameters[i]) == 2)
			{
				CYL_METHOD_TYPE = pcl::SAC_LMEDS;
			}
			else if (std::stoi(parameters[i]) == 3)
			{
				CYL_METHOD_TYPE = pcl::SAC_MSAC;
			}
			else if (std::stoi(parameters[i]) == 4)
			{
				CYL_METHOD_TYPE = pcl::SAC_RRANSAC;
			}
			else if (std::stoi(parameters[i]) == 5)
			{
				CYL_METHOD_TYPE = pcl::SAC_RMSAC;
			}
			else if (std::stoi(parameters[i]) == 6)
			{
				CYL_METHOD_TYPE = pcl::SAC_MLESAC;
			}
			else
			{
				CYL_METHOD_TYPE = INF;
			}
		}
				break;
		case 2: {
			if (std::stoi(parameters[i]) == 0)
			{
				CONE_METHOD_TYPE = pcl::SAC_RANSAC;
			}
			else if (std::stoi(parameters[i]) == 1)
			{
				CONE_METHOD_TYPE = pcl::SAC_PROSAC;
			}
			else if (std::stoi(parameters[i]) == 2)
			{
				CONE_METHOD_TYPE = pcl::SAC_LMEDS;
			}
			else if (std::stoi(parameters[i]) == 3)
			{
				CONE_METHOD_TYPE = pcl::SAC_MSAC;
			}
			else if (std::stoi(parameters[i]) == 4)
			{
				CONE_METHOD_TYPE = pcl::SAC_RRANSAC;
			}
			else if (std::stoi(parameters[i]) == 5)
			{
				CONE_METHOD_TYPE = pcl::SAC_RMSAC;
			}
			else if (std::stoi(parameters[i]) == 6)
			{
				CONE_METHOD_TYPE = pcl::SAC_MLESAC;
			}
			else
			{
				CONE_METHOD_TYPE = INF;
			}
		}
				break;
		case 3: {PLANE_TOL = std::stod(parameters[i]); }
				break;
		case 4: {CYL_TOL = std::stod(parameters[i]); }
				break;
		case 5: {CONE_TOL = std::stod(parameters[i]); }
				break;
		case 6: {MIN_INLIERS = std::stod(parameters[i]); }
				break;
		case 7: {MST_RATIO = std::stod(parameters[i]); }
				break;
		case 8: {MST_RATIO2 = std::stod(parameters[i]); }
				break;
		case 9: {PLANE_MAX_NUM_ITER = std::stoi(parameters[i]); }
				break;
		case 10: {CYL_MAX_NUM_ITER = std::stoi(parameters[i]); }
				 break;
		case 11: {CONE_MAX_NUM_ITER = std::stoi(parameters[i]); }
				 break;
		case 12: {K_FOR_NORMAL_SEARCH = std::stoi(parameters[i]); }
				 break;
		case 13: {CYL_WEIGHT_NORMAL_DISTANCE = std::stod(parameters[i]); }
				 break;
		case 14: {CONE_WEIGHT_NORMAL_DISTANCE = std::stod(parameters[i]); }
				 break;
		case 15: {CYL_MIN_RADIUS_LIMIT = std::stod(parameters[i]); }
				 break;
		case 16: {CYL_MAX_RADIUS_LIMIT = std::stod(parameters[i]); }
				 break;
		case 17: {CONE_MIN_RADIUS_LIMIT = std::stod(parameters[i]); }
				 break;
		case 18: {CONE_MAX_RADIUS_LIMIT = std::stod(parameters[i]); }
				 break;
		case 19: {CONE_MIN_OPENING_ANGLE = std::stod(parameters[i]); }
				 break;
		case 20: {CONE_MAX_OPENING_ANGLE = std::stod(parameters[i]); }
				 break;
		case 21: {TOLERANCE_FOR_ADJACENT_CHAINS = std::stod(parameters[i]); }
				 break;
		case 22: {TOLERANCE_FOR_APPEX = std::stod(parameters[i]); }
				 break;
		case 23: {POISSON_DISK_SAMPLING_RADIUS = std::stoi(parameters[i]); }
				 break;
		case 24: {PATH_TEMPORARY_FILE = parameters[i]; }
				 break;
		case 25: {PLANE_AREA_PERCENTAGE = std::stod(parameters[i]); }
				 break;
		case 26: {CYL_AREA_PERCENTAGE = std::stod(parameters[i]); }
				 break;
		case 27: {CONE_AREA_PERCENTAGE = std::stod(parameters[i]); }
				 break;
		case 28: {PERCENTAGE = std::stod(parameters[i]); }
				 break;

		}
	}
}

bool readParameterFile(std::string parameterTxtFilePath)
{
	std::ifstream infile(parameterTxtFilePath);
	std::string line;
	std::vector<std::string> parameters;
	const size_t num_parameters = 29;
	//size_t n_lines = 0;

	while (std::getline(infile, line))
	{
		if (line[0] != '#')
		{
			parameters.push_back(line);
			//std::istringstream iss(line);
			//int a, b;
			//if (!(iss >> a >> b)) { break; } // error

			// process pair (a,b)
		}
		//++n_lines;
	}

	////Print parameters in cmd window:
	//for (std::size_t i = 0; i < parameters.size(); ++i)
	//{
	//	std::cerr << parameters[i] << std::endl;
	//}
	//std::cerr << "numero di paramentri: " << parameters.size() << std::endl;

	//Check if the parameter file contains the right number of parameters
	if (parameters.size() != num_parameters)
	{
		std::cerr << "ERROR in the text file of parameters!" << std::endl;
		return false;
	}
	//std::cerr << "The parameter file is correct." << std::endl;
	setParameters(parameters);
	return true;
}

std::vector<pcl::PointXYZ> ReadTXTFileAsPointCloud(const char *cfilename)
{
	std::vector<pcl::PointXYZ> subsampled_points;
	char a[100];
	double x = 0, y = 0, z = 0;
	int count = 0;

	std::ifstream in(cfilename, std::ifstream::in);
	do
	{
		in.getline(a, 100, '\n');

		if (sscanf(a, "%lf%lf%lf", &x, &y, &z) == 3)
		{
			subsampled_points.push_back(pcl::PointXYZ(x, y, z));
		}
		count++;
	} while (!in.eof());

	return subsampled_points;
}

int ReadTXTFileAsPointCloud(const char *cfilename, pcl::PointCloud<pcl::PointXYZ>::Ptr &points_list, std::vector<double> &plane_normal)
{
	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	char a[100];
	double x = 0, y = 0, z = 0;
	int count = 0;

	std::ifstream in(cfilename, std::ifstream::in);

	in.getline(a, 100, '\n');

	if (sscanf(a, "%lf%lf%lf", &x, &y, &z) == 3)
	{
		plane_normal.push_back(x);
		plane_normal.push_back(y);
		plane_normal.push_back(z);
	}


	do
	{
		in.getline(a, 100, '\n');

		if (sscanf(a, "%lf%lf%lf", &x, &y, &z) == 3)
		{
			(*points_list).push_back(pcl::PointXYZ(x, y, z));
		}
		count++;
	} while (!in.eof());

	return 0;
}


std::vector<pcl::Normal> ReadTXTFilePointNormal(const char *cfilename)
{
	std::vector<pcl::Normal> subsampled_points;
	char a[100];
	double x = 0, y = 0, z = 0;
	int count = 0;

	std::ifstream in(cfilename, std::ifstream::in);
	do
	{
		in.getline(a, 100, '\n');

		if (sscanf(a, "%lf%lf%lf", &x, &y, &z) == 3)
		{
			subsampled_points.push_back(pcl::Normal(x, y, z));
		}
		count++;
	} while (!in.eof());

	return subsampled_points;
}

std::vector<int> ReadTXTNumber(std::string cfilename)
{
	std::vector<int> number;
	char a[100];
	double x = 0, y = 0, z = 0;
	int count = 0;

	std::ifstream in(cfilename, std::ifstream::in);
	do
	{
		in.getline(a, 100, '\n');

		if (sscanf(a, "%lf", &x) == 1)
		{
			number.push_back(x);
		}
		count++;
	} while (!in.eof());

	return number;
}