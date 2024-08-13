//
// Created by Arkadiy Simonov on 08.08.24.
//

#include "IntensityMap.h"



IntensityMap ReadHDF5(string filename)
{
    H5File file( filename, H5F_ACC_RDONLY );
    DataSet dataset = file.openDataSet( "data" );

    DataSpace dataspace = dataset.getSpace();

    hsize_t dims_out[3];
    int no_dimensions = dataspace.getSimpleExtentDims( dims_out, NULL);

    //\TODO: assert that there is not too many dimensions in input vector

    for(int i=no_dimensions; i<3; i++)
        dims_out[i]=1; //fill additional dimensions if input matrix has less then 3 dimensions

    double * temp_map_holder = (double*) malloc(dims_out[0]*dims_out[1]*dims_out[2] * sizeof(double));

    dataset.read(temp_map_holder, PredType::NATIVE_DOUBLE);

    IntensityMap result(dims_out[0], dims_out[1], dims_out[2]); ///\TODO: copy data directly into versa, not like this, with additional array
    for(int i=0; i<dims_out[0]*dims_out[1]*dims_out[2]; i++)
        result.at(i) = temp_map_holder[i];

    free(temp_map_holder);

    return result;
}

template<typename T>
DataType getH5Type() {}

template<>
DataType getH5Type<int> () {
    return PredType::NATIVE_INT;
}
template<>
DataType getH5Type<bool> () {
    return PredType::NATIVE_HBOOL;
}
template<>
DataType getH5Type<double> () {
    return PredType::NATIVE_DOUBLE;
}

template <typename T>
void creadeAndWriteDataset(H5File& file, string datasetName, T* data, hsize_t n, hsize_t* dims) {
    DataSpace dataspace( n, dims );
    DataSet dataset = file.createDataSet( datasetName, getH5Type<T>(), dataspace );
    dataset.write( data, getH5Type<T>() );
}

template <typename T>
void writeConstant(H5File& file, string datasetName, T data) {
    H5::DataSet ds = file.createDataSet(datasetName, getH5Type<T>(), H5::DataSpace(H5S_SCALAR));
    ds.write(&data, getH5Type<T>());
}

template <typename T, std::size_t N>
void writeVector(H5File& file, string datasetName, af::tiny_plain<T,N> v) {
    T dataAsPlaneArray[N];
    for(int i=0; i<N; ++i)
        dataAsPlaneArray[i]=v[i];

    hsize_t sz[N]={N};
    creadeAndWriteDataset<T>(file, datasetName, dataAsPlaneArray, 1, sz);
}

void writeFormatString(H5File& file) {
    string format = "Yell 1.0";

    H5::StrType h5stringType(H5::PredType::C_S1, H5T_VARIABLE); // + 1 for trailing zero
    H5::DataSet ds = file.createDataSet("format", h5stringType, H5::DataSpace(H5S_SCALAR));
    ds.write(format, h5stringType);
}

void WriteHDF5(string filename, IntensityMap& input)
{
    H5File file( filename, H5F_ACC_TRUNC );

    hsize_t dimsf[3];
    //Fill it from input
    vec3<int> dims = input.size();
    for(int i=0; i<3; i++)
        dimsf[i] = dims[i];

    double * temp_map_holder = (double*) malloc(dimsf[0]*dimsf[1]*dimsf[2] * sizeof(double));

    for(int i=0; i<dimsf[0]*dimsf[1]*dimsf[2]; i++)
        temp_map_holder[i] = input.at(i);

    creadeAndWriteDataset<double>(file, "data", temp_map_holder, 3, dimsf);

    free(temp_map_holder);

    Grid& g=input.grid;
    // Write out the grid settings
    //creadeAndWriteDataset<double>(file, "lower_limits", {g.lower_limits[0], g.lower_limits[1], g.lower_limits[2]}, 1, {3});
    writeVector(file, "lower_limits", g.lower_limits);
    writeConstant(file, "is_direct", !g.reciprocal_flag);
    writeVector(file, "step_sizes", g.grid_steps);
    cctbx::uctbx::unit_cell cell = g.cell;
    if(!g.reciprocal_flag)
        cell=cell.reciprocal();

    writeVector(file, "unit_cell", cell.parameters());
    writeFormatString(file);

}