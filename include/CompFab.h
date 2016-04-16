//
//  CompFab.h
//  voxelizer
//
//
//

#ifndef voxelizer_CompFab_h
#define voxelizer_CompFab_h

#define EPSILON 1e-9

#include <cmath>

namespace CompFab
{
    
    //Data Types
    typedef struct Vec3Struct
    {
        
        Vec3Struct();
        Vec3Struct(double x, double y, double z);

        union
        {
            double m_pos[3];
            struct { double m_x,m_y,m_z; };
        };
        
        inline double & operator[](unsigned int index) { return m_pos[index]; }
        inline const double & operator[](unsigned int index) const { return m_pos[index]; }
        inline void operator+=(const Vec3Struct &a)
        {
            m_x += a.m_x;
            m_y += a.m_y;
            m_z += a.m_z;
        }
        
        void normalize();
        
    }Vec3;

    //Data Types
    typedef struct Vec3iStruct
    {
        
        Vec3iStruct();
        Vec3iStruct(double x, double y, double z);
        union
        {
            int m_pos[3];
            struct {int m_x,m_y,m_z;};
        };
        
        inline int & operator[](unsigned int index) { return m_pos[index]; }
        inline const int & operator[](unsigned int index) const { return m_pos[index]; }
        
    }Vec3i;

    //Data Types
    typedef struct Vec2fStruct
    {
        
        Vec2fStruct();
        Vec2fStruct(double x, double y);
        
        union
        {
            float m_pos[2];
            struct { float m_x,m_y; };
        };
        
        inline float & operator[](unsigned int index) { return m_pos[index]; }
        inline const float & operator[](unsigned int index) const { return m_pos[index]; }
        
    }Vec2f;

    
    //NOTE: Ray direction must be normalized
    typedef struct RayStruct
    {
        
        RayStruct();
        RayStruct(Vec3 &origin, Vec3 &direction);
        
        Vec3 m_origin;
        Vec3 m_direction;
        
    } Ray;
    
    typedef struct TriangleStruct
    {
        
        TriangleStruct(Vec3 &v1, Vec3 &v2,Vec3 &v3);
        
        Vec3 m_v1, m_v2, m_v3;
        
    }Triangle;
    
    //Some useful operations
    //Compute v1 - v2
    Vec3 operator-(const Vec3 &v1, const Vec3 &v2);
    
    Vec3 operator+(const Vec3 &v1, const Vec3 &v2);
    
    //Cross Product
    Vec3 operator%(const Vec3 &v1, const Vec3 &v2);
    
    //Dot Product
    double operator*(const Vec3 &v1, const Vec3 &v2);
    
    
    //Grid structure for Voxels
    typedef struct VoxelGridStruct {
        //Square voxels only
        VoxelGridStruct(Vec3 lowerLeft, unsigned int dimX, unsigned int dimY, unsigned int dimZ, double spacing);
        ~VoxelGridStruct();

        inline bool & isInside(int ii, int jj, int kk) {
            unsigned int i = (unsigned int) ii;
            unsigned int j = (unsigned int) jj;
            unsigned int k = (unsigned int) kk;
            
            return m_insideArray[k*(m_dimX*m_dimY)+j*m_dimY + i];
        }  

        inline void setIsInside(unsigned int i, unsigned int j, unsigned int k) {
            m_insideArray[k*(m_dimX*m_dimY)+j*m_dimY + i] = true;
        }      

        inline bool & isOnSurface(unsigned int i, unsigned int j, unsigned int k) {
            return m_surfaceArray[k*(m_dimX*m_dimY)+j*m_dimY + i];
        }

        inline void setOnSurface(unsigned int i, unsigned int j, unsigned int k) {
            m_surfaceArray[k*(m_dimX*m_dimY)+j*m_dimY + i] = true;
        }  

        inline unsigned int & getPieceNum(unsigned int i, unsigned int j, unsigned int k) {
            return m_pieceNumArray[k*(m_dimX*m_dimY)+j*m_dimY + i];
        }

        inline void setPieceNum(unsigned int i, unsigned int j, unsigned int k, unsigned int pieceNum) {
            m_pieceNumArray[k*(m_dimX*m_dimY)+j*m_dimY + i] = pieceNum;
        }  

        inline int getArrayIndexByCoordinate(unsigned int i, unsigned int j, unsigned int k) {
            return k*(m_dimX*m_dimY)+j*m_dimY + i;
        }

        // colors
        inline double getVoxelColor_r(unsigned int i, unsigned int j, unsigned int k) {
            return m_color_array[getArrayIndexByCoordinate(i,j,k)*3];
        }
        inline double getVoxelColor_g(unsigned int i, unsigned int j, unsigned int k) {
            return m_color_array[getArrayIndexByCoordinate(i,j,k)*3 + 1];
        }
        inline double getVoxelColor_b(unsigned int i, unsigned int j, unsigned int k) {
            return m_color_array[getArrayIndexByCoordinate(i,j,k)*3 + 2];
        }

        inline void setVoxelColor(unsigned int i, unsigned int j, unsigned int k, double r, double g, double b) {
            m_color_array[getArrayIndexByCoordinate(i,j,k)*3] = r;
            m_color_array[getArrayIndexByCoordinate(i,j,k)*3 + 1] = g;
            m_color_array[getArrayIndexByCoordinate(i,j,k)*3 + 2] = b;
        }

        bool *m_insideArray;
        bool *m_surfaceArray;
        unsigned int m_dimX, m_dimY, m_dimZ, m_size;
        double m_spacing;
        Vec3 m_lowerLeft;

        // piece numbers
        unsigned int *m_pieceNumArray;

        // color info
        // 3* # vertices for (r,g,b)
        double *m_color_array;
        
    } VoxelGrid;
}



#endif
