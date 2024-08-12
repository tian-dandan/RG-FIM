#include <vector>
#include <string>
#include <map>
#include <array>
#include <queue>
#include <iostream>

using namespace std;
class Region {
public:
    Region();
    int ID;
    std::vector<std::pair<int, int>> pixels;
    std::vector<std::pair<int, int>> boundary;
    int TYPE;
    int value;
    std::string growType;
    std::string shapeType;
    double accessibility;

    // Planimetric attributes
    long pixel_area;
    long area;
    int centroid_X;
    int centroid_Y;
    int perimeter;
    float thickness;
    float length;
    float width;

    // Shape attributes
    float shapeIndex;
    float compactness;
    float elongatedness;
    float asymmetry;
    float orientation;
    float fractal;
    float rectangularity;
    float ellipticity;
    float triangularity;

    // Surface attributes
    float mean_elv_1, mean_slp_1, mean_asp_1, mean_cur_1, mean_curvpro_1, mean_curvpl_1;
    float mean_elv_2, mean_slp_2, mean_asp_2, mean_cur_2, mean_curvpro_2, mean_curvpl_2;

    // Volumetric attributes
    float mean_depth, max_depth, sdv_depth, volume;

    // Thematic attributes

    // Drawing attributes
    // Rectangle
    float tl_x, tl_y, tr_x, tr_y, bl_x, bl_y, br_x, br_y;
    // Ellipse
    double a, b;

    double max_distance_centr;
    double mean_distance_centr;
    double mean_dist_bnd;

    // Statistic descriptors
    double min_x, min_y, max_x, max_y;
    double var_x, var_y, covar_xy;
    double lamda1, lamda2;

    // Moment descriptors
    double** m; // moments
    // double u[,];//central moments
    // double n[,];//normalized central moments
    // double phi_1,phi_2,phi_3,phi_4; //invariant moments
    void buildGeometry(const std::vector<std::vector<short>>& binaryImage);
    void print();
    std::string toJsonString() const; 
    void writeToImage(std::vector<std::vector<int>>& image) const;

    
};


class RegionManagement {
    
    std::vector<std::vector<short>> binaryImage;
    std::vector<std::vector<bool>> visited;
    int numRows;
    int numCols;
    int minRow,maxRow,minCol,maxCol;//row and column numbers for exploreBlob to use

    std::map<int, Region> blobs, blobs_back;
    void exploreBlob(int row, int col, Region& blob);
    bool isOnBoundary(int row, int col) const {
        return row == minRow || row == maxRow - 1 || col == minCol || col == maxCol - 1;
    }
    bool isWithinBoundary(int row, int col) const{
        return row >= minRow && row < maxRow && col >= minCol && col < maxCol;
    }

    std::vector<std::pair<int, int>> neighborDefinitions = {
        {-1, 0},  // Up
        {1, 0},   // Down
        {0, -1},  // Left
        {0, 1}    // Right
        // You can add more definitions based on your needs
    };    
public:
    RegionManagement():numRows(0),numCols(0){};
    RegionManagement(const std::vector<std::vector<short>>& binaryImage);
    string ExtractBoundary(int index);

    void printRegion(int index)
    {
        std::string msg;
        Region region = blobs[index];
        region.print();
        return;
    }
    void printAllRegions()
    {
        for(auto it = blobs.begin(); it != blobs.end();it++){
            it->second.print();
        }
    }
    std::vector<std::vector<int>> writeImage()const;
    void regionGeneration();
    void regionGeneration2(int nR, int nC);// the algorithm that splits the image to subregions and generate objects within each region

    void removeSmallObjs(int thresh_foreground, int thresh_background);
    void removeSmallObjs2(int thresh_foreground, int thresh_background, int nR, int nC);

    std::string writeRegionsToJSON() const;
};