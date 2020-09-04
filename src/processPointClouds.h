// PCL lib Functions for processing point clouds 

#ifndef PROCESSPOINTCLOUDS_H_
#define PROCESSPOINTCLOUDS_H_
#pragma once

#include <pcl/io/pcd_io.h>
#include <pcl/common/common.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/crop_box.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/common/transforms.h>
#include <iostream> 
#include <string>  
#include <vector>
#include <ctime>
#include <chrono>
#include "render/box.h"
#include <unordered_set>
#include "quiz/cluster/kdtree.h"
template<typename PointT>
class ProcessPointClouds {
public:

    //constructor
    ProcessPointClouds();
    //deconstructor
    ~ProcessPointClouds();

    void numPoints(typename pcl::PointCloud<PointT>::Ptr cloud);

    typename pcl::PointCloud<PointT>::Ptr FilterCloud(typename pcl::PointCloud<PointT>::Ptr cloud, float filterRes, Eigen::Vector4f minPoint, Eigen::Vector4f maxPoint);

    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> SeparateClouds(pcl::PointIndices::Ptr inliers, typename pcl::PointCloud<PointT>::Ptr cloud);

    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> SegmentPlane(typename pcl::PointCloud<PointT>::Ptr cloud, int maxIterations, float distanceThreshold);

    
    
    std::vector<typename pcl::PointCloud<PointT>::Ptr> Clustering(typename pcl::PointCloud<PointT>::Ptr cloud, float clusterTolerance, int minSize, int maxSize);

    Box BoundingBox(typename pcl::PointCloud<PointT>::Ptr cluster);

    void savePcd(typename pcl::PointCloud<PointT>::Ptr cloud, std::string file);

    typename pcl::PointCloud<PointT>::Ptr loadPcd(std::string file);

    std::vector<boost::filesystem::path> streamPcd(std::string dataPath);

  
    
    //own implementation
    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> SegmentPlane2(typename pcl::PointCloud<PointT>::Ptr cloud, int maxIterations, float distanceThreshold);
    //own implementation of clustering
    std::vector<typename pcl::PointCloud<PointT>::Ptr> Clustering2(typename pcl::PointCloud<PointT>::Ptr cloud, float clusterTolerance, int minSize, int maxSize);
    private:
    

    //cluster helper of eucledian clustering
    void clusterHelper(int indice, const std::vector<std::vector<float>>& points, std::vector<int>& cluster,std::vector<bool>& processed, KdTree* tree, float distanceTol)
    {
        processed[indice]=true;
        cluster.push_back(indice);
        std::vector<int> nearest = tree->search(points[indice],distanceTol);

        for(int id: nearest)
        {
            if(!processed[id])
                clusterHelper(id,points, cluster, processed, tree, distanceTol);
        }
    }
    //own implementation of eucledian clustering
    std::vector<std::vector<int>> euclideanCluster(const std::vector<std::vector<float>>& points, KdTree* tree, float distanceTol)
    {

        // TODO: Fill out this function to return list of indices for each cluster

        std::vector<std::vector<int>> clusters;
        std::vector<bool> processed(points.size(),false);

        int i=0;
        while(i<points.size())
        {
            if(processed[i])
            {
                i++;
                continue;
            }

            std::vector<int> cluster;
            clusterHelper(i,points,cluster,processed,tree,distanceTol);
            clusters.push_back(cluster);
            i++;

        }
        return clusters;

    }


    //own implementation of ransac

    std::unordered_set<int> Ransac2(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, int maxIterations, float distanceTol)
    {
        std::unordered_set<int> inliersResult;
        srand(time(NULL));
        
        while(maxIterations--)
        {
            std::unordered_set<int> inliers;
            while(inliers.size()<3)
                inliers.insert(rand()%(cloud->points.size()));
            float x1,y1,z1,x2,y2,z2,x3,y3,z3;

            auto itr = inliers.begin();
            x1=cloud->points[*itr].x;
            y1=cloud->points[*itr].y;
            z1=cloud->points[*itr].z;
            itr++;
            x2=cloud->points[*itr].x;
            y2=cloud->points[*itr].y;
            z2=cloud->points[*itr].z;
            itr++;
            x3=cloud->points[*itr].x;
            y3=cloud->points[*itr].y;
            z3=cloud->points[*itr].z;
            //vector v1 travels from point1 to point 2
            //vector v2 travels from point 1 to point 3
            // float v1[3],v2[3];
            float v1[]={x2-x1,y2-y1,z2-z1};
            float v2[]={x3-x1,y3-y1,z3-z1};
            // v1[0]=x2-x1;
            // v1[1]=y2-y1;
            // v1[2]=z2-z1;
            // v2[0]=x3-x1;
            // v2[1]=y3-y1;
            // v2[2]=z3-z1;
            //normal vector of v1,v2
            float normal[]={(v1[1]*v2[2])-(v1[2]*v2[1]),(v1[2]*v2[0])-(v1[0]*v2[2]),(v1[0]*v2[1])-(v1[1]*v2[0])};

            float a = normal[0];
            float b = normal[1];
            float c = normal[2];
            float d = -(a*x1+b*y1+c*z1);

            for(int index = 0; index<cloud->points.size();index++)
            {
                if(inliers.count(index)>0)
                    continue;

                //pcl::PointXYZI point = cloud->points[index];
                float x4 = cloud->points[index].x;
                float y4 = cloud->points[index].y;
                float z4 = cloud->points[index].z;
                float dist = fabs(a*x4+b*y4+c*z4+d)/sqrt(a*a+b*b+c*c);
                if(dist<=distanceTol)
                    inliers.insert(index);


            }
            if(inliers.size()>inliersResult.size())
            {
                inliersResult = inliers;
            }
        }
        return inliersResult;

    }
};
#endif /* PROCESSPOINTCLOUDS_H_ */