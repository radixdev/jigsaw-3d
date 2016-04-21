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
#include <map>
#include <vector>
#include <algorithm>

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

    //int operator<(const Vec3 &v1, const Vec3 &v2);
    
    
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
        // updatePiecesFromPuzzle(Puzzle) // should change pieceNumArray to be the same as the puzzle structure.

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

    // A single piece
    typedef struct PuzzlePieceStruct {
        //Square voxels only
        PuzzlePieceStruct(unsigned int m_id);
        ~PuzzlePieceStruct();

        // get id 
        inline unsigned int getID() {
            return m_id;
        }

        // get_voxels()
        inline std::vector<int> getVoxels() {
            return m_voxels;
        } 

        // add_voxels(std::vector<int>)
        // takes in a list of new Voxels to add, and adds them if they are no already in the set.
        inline std::vector<int> add_voxels(std::vector<int> newVoxels) {
            for (unsigned int i = 0; i < newVoxels.size(); i++) {
                if (std::find(m_voxels.begin(), m_voxels.end(), newVoxels[i]) == m_voxels.end()) {
                    m_voxels.push_back(newVoxels[i]);
                }
            }   

            return m_voxels;        
        }

        // remove_voxels(std::vector<int>)
        inline std::vector<int> remove_voxels(std::vector<int> badVoxels) {
            std::vector<int>::iterator found;
            std::vector<std::vector<int>::iterator> toErase;

            for (unsigned int i = 0; i < badVoxels.size(); i++) {
                found = std::find(m_voxels.begin(), m_voxels.end(), badVoxels[i]);
                if (found != m_voxels.end()) {
                    toErase.push_back(found);
                }
            } 

            for (unsigned int j = 0; j < toErase.size(); j++) {
                m_voxels.erase(toErase[j]);
            }

            return m_voxels;         
        }

        // based on location, should be newz*(g_voxelGrid -> m_dimX/PIECE_SIZE*g_voxelGrid ->m_dimY/PIECE_SIZE)+newy*g_voxelGrid ->m_dimY/PIECE_SIZE + newx
        unsigned int m_id;
        // keeps track of the voxels inside, k*(m_dimX*m_dimY)+j*m_dimY + i
        std::vector<int> m_voxels;
        
    } PuzzlePiece;

    // A puzzle with a lot of pieces
    typedef struct PuzzleStruct {
        //Square voxels only
        PuzzleStruct(unsigned int dimX, unsigned int dimY, unsigned int dimZ);
        ~PuzzleStruct();

        // get_pieces()
        inline std::vector<PuzzlePiece> get_pieces() {
            std::vector<PuzzlePiece> output;
            for (std::map<unsigned int, PuzzlePiece>::iterator it = m_pieceList.begin(); it != m_pieceList.end(); it++) {
                output.push_back(it -> second);
            }
        }

        // add_piece(id, PuzzlePiece)
        inline void add_or_update_piece(unsigned int id, PuzzlePiece piece) {
            if (m_pieceList.count(id) <= 0) {
                // not in the map yet 
                m_pieceList.insert(std::map<unsigned int, PuzzlePiece>::value_type(id, piece));
            } else {
                m_pieceList[id] = piece;
            }
            return;
        }

        // remove_piece(id)
        inline void remove_piece(unsigned int id) {
            if (m_pieceList.count(id) > 0) {
                // not in the map yet 
                m_pieceList.erase(id);
            }
            return;
        }

        // has_piece_at(i, j, k)
        inline bool has_piece_at(int i, int j, int k) {
            return m_pieceList.count(k*(m_dimX*m_dimY)+j*m_dimY + i) > 0;
        }

        // merge_pieces(id1, id2)

        unsigned int m_dimX, m_dimY, m_dimZ, m_size;
        std::map<unsigned int, PuzzlePiece> m_pieceList;
        
    } Puzzle;
}

void loadVoxelGrid(const char * filename, CompFab::VoxelGrid & voxelGrid);
void saveVoxelGrid(const char * filename, CompFab::VoxelGrid & voxelGrid);

// if the VOX file is present, return true
bool isVoxelCacheFilePresent(const char * filename);

#endif
