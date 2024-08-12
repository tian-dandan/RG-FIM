
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Python.h>

#include <pybind11/numpy.h>
#include <iostream>
#include "ImageObject.h"

using namespace std;

namespace py = pybind11;

template <typename T>
vector< vector<short> > createShort2DArrayFromNumpy(py::array_t<T> input_array)
{
	// Get the shape of the input array
	auto shape = input_array.shape();
	size_t rows = shape[0];
	size_t cols = shape[1];

	// Create a new short array with the same shape
    vector< vector<short> > result;
	result.resize(rows,vector<short>(cols));
	// Get pointers to the data of both arrays
	T *input_ptr = static_cast<T *>(input_array.request().ptr);

	// Copy the data from the input array to the short array
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			result[i][j] = input_ptr[i * cols + j];
		}
	}

	return result;
}
template <typename T>
vector< vector<int> > createInt2DArrayFromNumpy(py::array_t<T> input_array)
{
	// Get the shape of the input array
	auto shape = input_array.shape();
	size_t rows = shape[0];
	size_t cols = shape[1];

	// Create a new short array with the same shape
    vector< vector<int> > result;
	result.resize(rows,vector<int>(cols));
	// Get pointers to the data of both arrays
	T *input_ptr = static_cast<T *>(input_array.request().ptr);

	// Copy the data from the input array to the short array
	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < cols; ++j)
		{
			result[i][j] = input_ptr[i * cols + j];
		}
	}

	return result;
}
py::array_t<int> imageToPyArray(const std::vector<std::vector<short>>& inputVector) {
    // Determine the shape of the array
    int numRows = static_cast<int>(inputVector.size());
    int numCols = numRows > 0 ? static_cast<int>(inputVector[0].size()) : 0;

    // Create a Pybind11 array with the given shape
    auto resultArray = py::array_t<int>({numRows, numCols});

    // Get a pointer to the underlying data
    auto resultPtr = resultArray.mutable_unchecked<2>();

    // Copy the data from the vector to the Pybind11 array
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            resultPtr(i, j) = static_cast<int>(inputVector[i][j]);
			//cout<<static_cast<int>(inputVector[i][j]);
        }
    }

    return resultArray;
}


PYBIND11_MODULE(ImageObject, m)
{
	m.def("Create2DArrayFromNumpy",&createInt2DArrayFromNumpy<int>,"Create 2d array of the integer type for c++ processing",py::arg("input_image"));
	m.def("Create2DArrayFromNumpy",&createInt2DArrayFromNumpy<unsigned char>,"Create 2d array of the byte type for c++ processing",py::arg("input_image"));
	m.def("imageToPyArray",&imageToPyArray,"Put the output image back to pyarray",py::arg("out_image"));
	//m.def("CreateBlobs", &CreateBlobs<int>,"Create blobs from a binary image", py::arg("input_image"));
	//m.def("CreateBlobs", &CreateBlobs<unsigned char>,"Create blobs from a binary image", py::arg("input_image"));
	m.def("foo", []() {
        return "Hello, World!";
    });
	py::class_<RegionManagement>(m, "RegionManagement")
		.def(py::init<>())
		.def(py::init<const std::vector<std::vector<short>>&>())
		.def("regionGeneration",&RegionManagement::regionGeneration)
		.def("regionGeneration2",&RegionManagement::regionGeneration2)
		.def("removeSmallObjs",&RegionManagement::removeSmallObjs,py::arg("threshold_forground"),py::arg("threshold_background"))
		.def("removeSmallObjs2",&RegionManagement::removeSmallObjs2)
		.def("writeJSON",&RegionManagement::writeRegionsToJSON)
		.def("extractBoundary",&RegionManagement::ExtractBoundary,py::arg("ID"))
		.def("writeImage",&RegionManagement::writeImage);

	#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
	#else
		m.attr("__version__") = "dev";
	#endif
}