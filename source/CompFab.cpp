//
//  CompFab.cpp
//  voxelizer
//
//
//
#include <fstream>
#include <iostream>
#include "../include/CompFab.h"


using namespace CompFab;





CompFab::Vec3Struct::Vec3Struct()
{
    m_x = m_y = m_z = 0.0;
}

CompFab::Vec3Struct::Vec3Struct(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

void CompFab::Vec3Struct::normalize() {
    
    double magnitude = sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
    
    if(magnitude > EPSILON)
    {
        m_x /= magnitude;
        m_y /= magnitude;
        m_z /= magnitude;
    }
}

//Data Types
CompFab::Vec3iStruct::Vec3iStruct()
{
    m_x = m_y = m_z = 0.0;
}

CompFab::Vec3iStruct::Vec3iStruct(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

CompFab::Vec2fStruct::Vec2fStruct()
{
    m_x = m_y = 0.0;
}

CompFab::Vec2fStruct::Vec2fStruct(double x, double y)
{
    m_x = x;
    m_y = y;
}

CompFab::RayStruct::RayStruct()
{
    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;
    m_direction[0] = 1.0;
    m_direction[1] = m_direction[2] = 0.0;
}

CompFab::RayStruct::RayStruct(Vec3 &origin, Vec3 &direction)
{
    m_origin = origin;
    m_direction = direction;
}

CompFab::TriangleStruct::TriangleStruct(Vec3 &v1, Vec3 &v2,Vec3 &v3)
{
    m_v1 = v1;
    m_v2 = v2;
    m_v3 = v3;
}

CompFab::Vec3 CompFab::operator-(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[0] - v2[0];
    v3[1] = v1[1] - v2[1];
    v3[2] = v1[2] - v2[2];

    return v3;
}

CompFab::Vec3 CompFab::operator+(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[0] + v2[0];
    v3[1] = v1[1] + v2[1];
    v3[2] = v1[2] + v2[2];
    
    return v3;
}


//Cross Product
Vec3 CompFab::operator%(const Vec3 &v1, const Vec3 &v2)
{
    Vec3 v3;
    v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v3[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return v3;
}

//Dot Product
double CompFab::operator*(const Vec3 &v1, const Vec3 &v2)
{
    return v1.m_x*v2.m_x + v1.m_y*v2.m_y+v1.m_z*v2.m_z;
}

//Comparison
/*int operator<(const Vec3 &v1, const Vec3 &v2)
{
    if (v1.m_x > v2.m_x) {
        return 1;
    } else if (v1.m_x < v2.m_x) {
        return -1;
    } else if (v1.m_y > v2.m_y) {
        return 1;
    } else if (v1.m_y < v2.m_y) {
        return -1;
    } else if (v1.m_z > v2.m_z) {
        return 1;
    } else if (v1.m_z < v2.m_z) {
        return -1;
    } else {
        return 0;
    }
}*/


//Grid structure for Voxels
CompFab::VoxelGridStruct::VoxelGridStruct(Vec3 lowerLeft, unsigned int dimX, unsigned int dimY, unsigned int dimZ, double spacing)
{
    m_lowerLeft = lowerLeft;
    m_dimX = dimX;
    m_dimY = dimY;
    m_dimZ = dimZ;
    m_size = dimX*dimY*dimZ;
    m_spacing = spacing;
    
    //Allocate Memory
    m_insideArray = new bool[m_size];
    m_surfaceArray = new bool[m_size];
    m_pieceNumArray = new unsigned int[m_size];
    m_color_array = new double[m_size*3];

    for(unsigned int ii=0; ii<m_size; ++ii) {
        m_insideArray[ii] = false;
        m_surfaceArray[ii] = false;
        m_pieceNumArray[ii] = 0;

        // colors
        m_color_array[ii*3] = 0.0;
        m_color_array[ii*3 + 1] = 0.0;
        m_color_array[ii*3 + 2] = 0.0;
    }
    
}

CompFab::VoxelGridStruct::~VoxelGridStruct()
{
    delete[] m_insideArray;
    delete[] m_surfaceArray;    
    delete[] m_pieceNumArray;
}

//Class structure for puzzle pieces
CompFab::PuzzlePieceStruct::PuzzlePieceStruct(unsigned int id) {
    m_id = id;
}

CompFab::PuzzlePieceStruct::~PuzzlePieceStruct() {
}


std::vector<CompFab::PuzzlePiece> CompFab::PuzzlePieceStruct::check_contiguous(){
    std::vector<CompFab::PuzzlePiece> output;
    return output;
}


//Class structure for puzzles
CompFab::PuzzleStruct::PuzzleStruct(unsigned int dimX, unsigned int dimY, unsigned int dimZ) {
    m_dimX = dimX;
    m_dimY = dimY;
    m_dimZ = dimZ;
    m_size = dimX*dimY*dimZ;
}

CompFab::PuzzleStruct::~PuzzleStruct() {
}

void saveGrid(std::ostream & out, CompFab::VoxelGrid & voxelGrid) {
    // write bool per line
    for (int i=0; i < voxelGrid.m_size; i++) {
        bool isInside = voxelGrid.m_insideArray[i];
        if (isInside) {
            out<<"1\n";
        } else {
            out<<"0\n";
        }
    }

    out<<"#end\n";
}

void readGrid(std::istream & f, CompFab::VoxelGrid & voxelGrid) {
    std::string line;
    std::string endLine("#end");
    std::string trueLine("1");
    std::string falseLine("0");

    int i = 0;
    while(true) {
        std::getline(f,line);

        // if we hit the end line, break out
        if(std::string::npos!=line.find(endLine)) {
            break;
        }

        if(std::string::npos!=line.find(trueLine)) {
            voxelGrid.m_insideArray[i] = true;
        }

        if(std::string::npos!=line.find(falseLine)) {
            voxelGrid.m_insideArray[i] = false;
        }

        i++;
    }

    std::cout<<"lines read from VOX file: "<<i<<"\n";

}

void loadVoxelGrid(const char * filename, CompFab::VoxelGrid & voxelGrid) {
    std::cout<<"loading VOX file "<<filename<<"\n";
    std::ifstream f;
    f.open(filename);
    if(!f.good()) {
        std::cout<<"Error: cannot open VOX file "<<filename<<"\n";
        return;
    }


    readGrid(f, voxelGrid);
    f.close();
    std::cout<<"done loading grid cache file\n";
}

void saveVoxelGrid(const char * filename, CompFab::VoxelGrid & voxelGrid) {
    std::ofstream out(filename);
    if(!out.good()){
        std::cout<<"cannot open voxel grid output file"<<filename<<"\n";
        return;
    }

    saveGrid(out, voxelGrid);
    out.close();
    std::cout<<"done saving grid cache file\n";
}

bool isVoxelCacheFilePresent(const char * filename) {
    std::ifstream f;
    f.open(filename);
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    } 
}

