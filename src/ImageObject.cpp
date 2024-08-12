#include "ImageObject.h"
#include <queue>
#include <stdexcept>
#include <sstream>
#include <cfloat>
#include <math.h>
#include <limits.h>
#define ID_BACKGROUND 0
using namespace std;

void RegionManagement::regionGeneration() {
    if(numRows == 0 || numCols == 0){
		cout<<"The class is not initialized properly\n";
        throw exception();
    }
	minRow = minCol = 0;
	maxRow = numRows;
	maxCol = numCols;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            // Modify individual bits using the [] operator
            visited[i][j] = false;
        }
    }

    blobs.clear();
	blobs_back.clear();
	cout<<"Scanning for foreground objects"<<endl;

    int currentId = 1;  // Start blob IDs from 1
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            if (binaryImage[i][j] != ID_BACKGROUND && !visited[i][j]) {
                Region blob;
                blob.ID = currentId;
                blob.value = binaryImage[i][j];
				//cout<<"exploreBlob "<<currentId<<endl;

                exploreBlob(i, j, blob);
                blob.buildGeometry(binaryImage);
                blobs[currentId] = blob;
                currentId++;
            }
        }
    }
    cout<<"Number of object regions created: "<<blobs.size()<<endl;

	cout<<"Scanning for background objects"<<endl;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            // Modify individual bits using the [] operator
            visited[i][j] = false;
        }
    }
	//scan for background objects
	for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            if (binaryImage[i][j] == ID_BACKGROUND && !visited[i][j]) {
                Region blob;
                blob.ID = currentId;
                blob.value = binaryImage[i][j];
                exploreBlob(i, j, blob);
                blob.buildGeometry(binaryImage);
                //blob.print();
                blobs_back[currentId] = blob;
                currentId++;
            }
        }
    }
	cout<<"Number of background regions created: "<<blobs_back.size()<<endl;
}
void RegionManagement::regionGeneration2(int nR, int nC) {
    if(numRows == 0 || numCols == 0){
		cout<<"The class is not initialized properly\n";
        throw exception();
    }
	minRow = minCol = 0;
	maxRow = numRows;
	maxCol = numCols;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            // Modify individual bits using the [] operator
            visited[i][j] = false;
        }
    }

    blobs.clear();
	blobs_back.clear();
	cout<<"Scanning for foreground objects"<<endl;
	int r_interval = numRows / nR;
	int c_interval = numCols / nC;

    int currentId = 1;  // Start blob IDs from 1
	for(int r = 0; r < nR; r++){
		for(int c = 0; c < nC; c++){
			minRow = r * r_interval;
			minCol = c * c_interval;
			maxRow = min(minRow + r_interval,numRows);
			maxCol = min(minCol + c_interval,numCols);
			for (int i = minRow; i < maxRow; ++i) {
				for (int j = minCol; j < maxCol; ++j) {
					if (binaryImage[i][j] != ID_BACKGROUND && !visited[i][j]) {
						Region blob;
						blob.ID = currentId;
						blob.value = binaryImage[i][j];
						exploreBlob(i, j, blob);
						blob.buildGeometry(binaryImage);
						blobs[currentId] = blob;
						currentId++;
					}
				}
			}
		}
	}
    cout<<"Number of object regions created: "<<blobs.size()<<endl;

	cout<<"Scanning for background objects"<<endl;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            // Modify individual bits using the [] operator
            visited[i][j] = false;
        }
    }
	//scan for background objects
	for(int r = 0; r < nR; r++){
		for(int c = 0; c < nC; c++){
			minRow = r * r_interval;
			minCol = c * c_interval;
			maxRow = min(minRow + r_interval,numRows);
			maxCol = min(minCol + c_interval,numCols);
			for (int i = minRow; i < maxRow; ++i) {
				for (int j = minCol; j < maxCol; ++j) {
					if (binaryImage[i][j] == ID_BACKGROUND && !visited[i][j]) {
						Region blob;
						blob.ID = currentId;
						blob.value = binaryImage[i][j];
						exploreBlob(i, j, blob);
						blob.buildGeometry(binaryImage);
						//blob.print();
						blobs_back[currentId] = blob;
						currentId++;
					}
				}
			}
		}
	}
	cout<<"Number of background regions created: "<<blobs_back.size()<<endl;
}

void RegionManagement::exploreBlob(int startRow, int startCol, Region& blob) {
        std::queue<std::pair<int,int>> pixelQueue;
        pixelQueue.push({startRow, startCol});
		int target_value = blob.value;
		visited[startRow][startCol] = true;

		//std::vector<std::pair<int, int>> neighbors = {
        //{-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1}};
        while (!pixelQueue.empty()) {
			auto currentPixel = pixelQueue.front();
			pixelQueue.pop();

			int row = currentPixel.first;
			int col = currentPixel.second;
			blob.pixels.emplace_back(row, col);
			bool isBound = false;
			// Enqueue neighbors
			for (const auto& neighbor : neighborDefinitions) {
				int nRow = row + neighbor.first;
				int nCol = col + neighbor.second;
				if (isWithinBoundary(nRow,nCol) && !visited[nRow][nCol] && binaryImage[nRow][nCol] == target_value){
					visited[nRow][nCol] = true;
					pixelQueue.push({nRow, nCol});
				}
				else if(!isWithinBoundary(nRow,nCol) ||binaryImage[nRow][nCol] != target_value){//this neighbor is different
					isBound = true;
				}
			}
			if(isBound) blob.boundary.emplace_back(row, col);

			//cout<<"Queue size: "<<pixelQueue.size()<<endl;

        }
        return;
}


RegionManagement::RegionManagement(const std::vector<std::vector<short>>& binaryImage) : binaryImage(binaryImage) {
    cout<<"Initlaizing the blob algorithm"<<endl;

    try{
        numRows = binaryImage.size();
        numCols = binaryImage[0].size();
        std::cout<<numRows<<","<<numCols<<std::endl;
        visited.resize(numRows, std::vector<bool>(numCols, false));
    }
    catch(std::exception &e){

        cout<<"something is wrong preventing RegionManagement constructor to complete"<<endl;
		
    }
    cout<<"Finished initlaizing the blob algorithm"<<endl;
}
std::vector<std::vector<int>> RegionManagement::writeImage()const
{

    vector< vector<int> > outImage;
    outImage.resize(numRows, std::vector<int>(numCols, 0));
    for(auto it = blobs.begin(); it != blobs.end();it++){
        it->second.writeToImage(outImage);
    }
    return outImage;
}

void RegionManagement::removeSmallObjs2(int thresh_foreground, int thresh_background, int nR, int nC)
{
	//initialize a new image with the background values



	//write the blobs pixels to a new image as foreground values if the size exceeds the threshold
	bool changed = true;
	while(changed){
		changed = false;
		for (auto& row : binaryImage) {
			for (auto& pixel : row) {
				pixel = 0;
			}
		}

		for(auto it = blobs.begin(); it != blobs.end();it++){
			Region r = it->second;
			if(r.pixel_area >= thresh_foreground)
			{
				for(auto it_pixel = r.pixels.begin(); it_pixel != r.pixels.end(); it_pixel++){
					int row = it_pixel->first;
					int col = it_pixel->second;
					binaryImage[row][col] = 1;//foreground, regardless the original values
				}
			}
			else{
				//cout<<"region id: "<<r.ID<<" was removed because the size "<<r.pixel_area<<" is less than "<<thresh_foreground<<endl;
				changed = true;
			}
		}
	
	

		//write the background blobs pixels to a new image as foreground values if the size is less than the threshold

		for(auto it = blobs_back.begin(); it != blobs_back.end();it++){
			Region r = it->second;
			if(r.pixel_area < thresh_background)
			{
				//cout<<"region id: "<<r.ID<<" was removed because the size "<<r.pixel_area<<" is less than "<<thresh_background<<endl;
				changed = true;
				for(auto it_pixel = r.pixels.begin(); it_pixel != r.pixels.end(); it_pixel++){
					int row = it_pixel->first;
					int col = it_pixel->second;
					binaryImage[row][col] = 1;//foreground, regardless the original values
				}
			}
		}
		//rescan the image to form blobs and blobs_back
		this->regionGeneration2(nR, nC);
	}
}

void RegionManagement::removeSmallObjs(int thresh_foreground, int thresh_background)
{
	//initialize a new image with the background values



	//write the blobs pixels to a new image as foreground values if the size exceeds the threshold
	bool changed = true;
	while(changed){
		changed = false;
		for (auto& row : binaryImage) {
			for (auto& pixel : row) {
				pixel = 0;
			}
		}

		for(auto it = blobs.begin(); it != blobs.end();it++){
			Region r = it->second;
			if(r.pixel_area >= thresh_foreground)
			{
				for(auto it_pixel = r.pixels.begin(); it_pixel != r.pixels.end(); it_pixel++){
					int row = it_pixel->first;
					int col = it_pixel->second;
					binaryImage[row][col] = 1;//foreground, regardless the original values
				}
			}
			else{
				//cout<<"region id: "<<r.ID<<" was removed because the size "<<r.pixel_area<<" is less than "<<thresh_foreground<<endl;
				changed = true;
			}
		}
	
	

		//write the background blobs pixels to a new image as foreground values if the size is less than the threshold

		for(auto it = blobs_back.begin(); it != blobs_back.end();it++){
			Region r = it->second;
			if(r.pixel_area < thresh_background)
			{
				//cout<<"region id: "<<r.ID<<" was removed because the size "<<r.pixel_area<<" is less than "<<thresh_background<<endl;
				changed = true;
				for(auto it_pixel = r.pixels.begin(); it_pixel != r.pixels.end(); it_pixel++){
					int row = it_pixel->first;
					int col = it_pixel->second;
					binaryImage[row][col] = 1;//foreground, regardless the original values
				}
			}
		}
		//rescan the image to form blobs and blobs_back
		this->regionGeneration();
	}
}
string RegionManagement::ExtractBoundary(int id){
	Region r = blobs[id];
	int i = 0;

	stringstream ss;
	ss<<"{";
    for (const auto &pair : r.boundary) {
		ss<<"\""<<i<<"\":"<<"{\"col\":"<<pair.second<<","<<"\"row\":"<<pair.first<<"},";
		i++;
    }
    ss.seekp(-1,ss.cur);

	ss<<"}";
    return ss.str();
}

std::string RegionManagement::writeRegionsToJSON() const
{
    std::stringstream ss;
    ss<<"{";

    for(auto it = blobs.begin(); it != blobs.end();it++){

        ss << "\"" << it->second.ID <<"\":" <<it->second.toJsonString()<< ",";
    }
    ss.seekp(-1,ss.cur);
    ss<<"}";
    return ss.str();
}
void Region::buildGeometry(const std::vector<std::vector<short>>& binaryImage)
{

    //1. Planimetric Attributes
//PIXEL AREA
	area = pixel_area = pixels.size();
    perimeter = boundary.size();

//1.1 CENTROID
	long sumX = 0;
	long sumY = 0;
    int i,j;
    if(boundary.size() == 0){
		cout<<"Boundary of the region is empty\n";
        throw new std::exception();
    }
    for (auto pixelIt = boundary.begin(); pixelIt != boundary.end(); ++pixelIt) {
        int row = pixelIt->first;
        int col = pixelIt->second;
        sumX = sumX + col;	
		sumY = sumY + row;

    }

	centroid_X = sumX/boundary.size();
	centroid_Y = sumY/boundary.size();

	////Statistic Descriptors
	////VAR_X VAR_Y COV_XY
	var_x = 0;
	var_y = 0;
	covar_xy = 0;
    for (auto pixelIt = pixels.begin(); pixelIt != pixels.end(); ++pixelIt) {
        int row = pixelIt->first;
        int col = pixelIt->second;
        var_x += (col - centroid_X)*(col - centroid_X);
        var_y += (row - centroid_Y)*(row - centroid_Y);
        covar_xy += (col - centroid_X) * (row - centroid_Y);

    }
	////EIGEN VALUE
	if(pixels.size() == 1)
	{
		lamda1 = 0.5;
		lamda2 = 0.5;
	}
	else
	{
		lamda1 = 0.5*(var_x + var_y + sqrt((var_x - var_y) * (var_x - var_y) + 4.0 * covar_xy * covar_xy));
		lamda2 = 0.5*(var_x + var_y - sqrt((var_x - var_y) * (var_x - var_y) + 4.0 * covar_xy * covar_xy));
	}

 //   //1.5 Length Width //WRONG!!
	length = lamda1;
	width = lamda2;

//Moment
    //binMoment
	array< array<double, 4>, 4> m;
	for(i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			double m_pq = 0;
			for (auto pixelIt = pixels.begin(); pixelIt != pixels.end(); ++pixelIt) {
                int row = pixelIt->first;
                int col = pixelIt->second;
				m_pq += pow((double)col,(double)i) * pow((double)row,(double)j);
			}
			m[i][j] = m_pq;
		}
	}
    //binCentralMoment
    	array< array<double, 4>, 4> u;
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			double u_pq = 0;
			for (auto pixelIt = pixels.begin(); pixelIt != pixels.end(); ++pixelIt) {
                int row = pixelIt->first;
                int col = pixelIt->second;		
				u_pq += pow((double)col - centroid_X, (double)i) * pow((double)row - centroid_Y, (double)j);
			}
			u[i][j] = u_pq;
		}		
	}



//1.5 Length and Width
	length = 2*sqrt((2.0*(u[2][0] + u[0][2] + sqrt((u[2][0] - u[0][2]) * (u[2][0] - u[0][2]) + 4.0 * u[1][1] * u[1][1])))/u[0][0]);
	width = 2*sqrt((2.0*(u[2][0] + u[0][2] - sqrt((u[2][0] - u[0][2]) * (u[2][0] - u[0][2]) + 4.0 * u[1][1] * u[1][1])))/u[0][0]);

	//I1
	double I1 = (u[2][0]*u[0][2]-u[1][1]*u[1][1])/pow(u[0][0],4);

//2. Shape attributes
const double PI  =3.141592653589793238463;

//SHAPE_INDEX-   SQRT(4*PI*A)/P
	shapeIndex = sqrt(4*PI*(area))/(perimeter);
//2.1 Compactness_index(CI) - 4*PI*A/pow(P,2)
	compactness = 4.0*PI*(area)/(perimeter * perimeter);
//2.2 Elongatedness;     /* ratio of length to width*/
	if(width !=0 )
	{
		elongatedness = length/ width;
	}
	else
	{
		//Assume single-pixel-width straight line. So length/width should be equal to size
		elongatedness = pixels.size();
	}
//2.3 ASYMMETRY
	if(length == 0){
		asymmetry = 0;
	}
	else{
		asymmetry = 1.0 - width/ length;

	}
//2.4 orientation - primary ORIENTATION of the image object
	if((u[1][1]<0.000001) && (u[1][1] > -0.000001))  //covariance = u[1][1]
	{
		orientation=0.0;
	}
	else
	{
		orientation = 0.5 * atan2(2*u[1][1],(u[2][0]-u[0][2]));
		//orientation = Math::Atan2(lamda1 - var_x, covar_xy) * 180.0 / PI;
	}


//2.6 rectangularity  /* A/(l*w) */
	if(width !=0 && length!=0)
	{
		rectangularity = area/(length * width);
	}
	else
	{
		rectangularity = 0;
	}

//2.7 ellipticity
	if(I1<=1/(16*PI*PI))
	{
		ellipticity = 16*PI*PI*I1;
	}
	else
	{
		ellipticity = 1/(16*PI*PI*I1);
	}
//2.8 triangularity;
	if(I1<=1.0/108.0)
	{
		triangularity = 108*I1;
	}
	else
	{
		triangularity = 1.0/(108*I1);
	}

//5.Draw bounding rectangle - vertex
	//Step 1: find edge point
	double V_width,V_length;
	double minWidth,maxWidth,minLength,maxLength;
	maxWidth = 0;	  maxLength = 0;
	minWidth = DBL_MAX;   minLength = DBL_MAX;
	tl_x = 0; tl_y = 0;
	tr_x = 0; tr_y = 0;
	bl_x = 0; bl_y = 0;
	br_x = 0; br_y = 0;
	int x1,x2,x3,x4;
	int y1,y2,y3,y4;
	double tan1 = tan(orientation);
	double cot;
	if(tan1!=0){cot = 1/tan1;}
	else{cot = DBL_MAX;}
    for (auto pixelIt = boundary.begin(); pixelIt != boundary.end(); ++pixelIt) {
        int thisRow = pixelIt->first;
        int thisCol = pixelIt->second;
		V_width = (thisRow - centroid_Y) - tan1*(thisCol - centroid_X);
		if(V_width>=0)
		{
			if(V_width>maxWidth)
			{
				maxWidth = V_width;
				x1 = thisCol;
				y1 = thisRow;
			}
		}
		else
		{
			if(V_width<minWidth)
			{
				minWidth = V_width;
				x2 = thisCol;
				y2 = thisRow;
			}
		}
					
		V_length = (thisRow - centroid_Y) + cot*(thisCol - centroid_X);
		if(V_length>=0)
		{
			if(V_length>maxLength)
			{
				maxLength = V_length;
				x3 = thisCol;
				y3 = thisRow;
			}
		}
		else
		{
			if(V_length<minLength)
			{
				minLength= V_length;
				x4 = thisCol;
				y4 = thisRow;
			}
		}

	}//End for boundary point

	//Step2: vertices of boundary rectangle
	tl_x = (x1*tan1 + x3*cot + y3 - y1)/(tan1 + cot);
	tl_y = (y1*cot + y3*tan1 + x3 - x1)/(tan1 + cot);
	tr_x = (x1*tan1 + x4*cot + y4 - y1)/(tan1 + cot);
	tr_y = (y1*cot + y4*tan1 + x4 - x1)/(tan1 + cot);
	bl_x = (x2*tan1 + x3*cot + y3 - y2)/(tan1 + cot);
	bl_y = (y2*cot + y3*tan1 + x3 - x2)/(tan1 + cot);
	br_x = (x2*tan1 + x4*cot + y4 - y2)/(tan1 + cot);
	br_y = (y2*cot + y4*tan1 + x4 - x2)/(tan1 + cot);
	


}//build geometry


Region::Region(){
    ID = 0;
    pixels.clear();
    boundary.clear();
}

void Region::print()
{
    cout<<toJsonString()<<endl;

}

std::string Region::toJsonString() const {
    std::stringstream ss;
    ss << "{"
    <<"\"ID\":" <<ID<<","
    << "\"pixel_area\":" << pixel_area << ","
    << "\"area\":" << area << ","
    << "\"centroid_X\":" << centroid_X << ","
    << "\"centroid_Y\":" << centroid_Y << ","
    << "\"perimeter\":" << perimeter << ","
    << "\"thickness\":" << thickness << ","
    << "\"length\":" << length << ","
    << "\"width\":" << width << ","
    << "\"shapeIndex\":" << shapeIndex << ","
    << "\"compactness\":" << compactness << ","
    << "\"elongatedness\":" << elongatedness << ","
    << "\"asymmetry\":" << asymmetry << ","
    << "\"orientation\":" << orientation << ","
    << "\"fractal\":" << fractal << ","
    << "\"rectangularity\":" << rectangularity << ","
    << "\"ellipticity\":" << ellipticity << ","
    << "\"triangularity\":" << triangularity          
    // Add the remaining properties similarly
    << "}";

    return ss.str();
}
void Region::writeToImage(std::vector<std::vector<int>>& image) const {
    int regionValue = value;

    for (const auto& pixel : pixels) {
        int row = pixel.first;
        int col = pixel.second;

        // Ensure the pixel coordinates are within the image bounds
        if (row >= 0 && row < image.size() && col >= 0 && col < image[0].size()) {
            image[row][col] = ID;
        }
    }
}


