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
#include <iostream>

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
        void scale(double val);
        double distanceSquaredTo(Vec3Struct* x);
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
    
    // A single piece
    typedef struct PuzzlePieceStruct {
        PuzzlePieceStruct(unsigned int id);
        PuzzlePieceStruct(unsigned int id, std::vector<Vec3i> voxels, unsigned int dimX, unsigned int dimY, unsigned int dimZ);
        ~PuzzlePieceStruct();

        // get id 
        inline unsigned int getID() {
            return m_id;
        }

        // get_voxels()
        inline std::vector<int>* getVoxels() {
            return &m_voxels;
        } 

        inline int size() {
            return m_voxels.size();
        }

        // // add_voxels(std::vector<int>)
        // // takes in a list of new Voxels to add, and adds them if they are no already in the set.
        // inline std::vector<int> add_voxels(std::vector<int> newVoxels) {
        //     std::cout << "there are " << newVoxels.size() << std::endl;
        //     for (unsigned int i = 0; i < newVoxels.size(); i++) {
        //         if (std::find(m_voxels.begin(), m_voxels.end(), newVoxels[i]) == m_voxels.end()) {
        //             m_voxels.push_back(newVoxels[i]);
        //         }
        //     }   

        //     return m_voxels;        
        // }

        inline void add_voxel(int voxel) {
            m_voxels.push_back(voxel);
        }

        inline std::vector<int> remove_voxel(int badVoxel) {
            std::vector<int>::iterator found;
            found = std::find(m_voxels.begin(), m_voxels.end(), badVoxel);
            if (found != m_voxels.end()) {
                m_voxels.erase(found);
            }
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

        // splits itself into puzzle pieces that are all contiguous
        // returns a vector of puzzle pieces
        // self, if it is all one piece, or many if there are separate chunks.  
        std::vector<PuzzlePieceStruct> check_contiguous(unsigned int gridX, unsigned int gridY, unsigned int gridZ);

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

        // get all the puzzle pieces
        inline std::vector<PuzzlePiece*> get_pieces() {
            std::vector<PuzzlePiece*> output;
            for (std::map<unsigned int, PuzzlePiece*>::iterator it = m_pieceList.begin(); it != m_pieceList.end(); it++) {
                output.push_back((it->second));
            }

            return output;
        }

        // adds a piece into the puzzle by id
        inline void add_piece(unsigned int id) {
            if (m_pieceList.count(id) <= 0) {
                // not in the map yet 
                m_pieceList.insert(std::make_pair<unsigned int, PuzzlePiece*>(id, new PuzzlePiece(id)));
            }
        }

        // remove_piece(id)
        inline void remove_piece(unsigned int id) {
            if (m_pieceList.count(id) > 0) {
                // not in the map yet 
                m_pieceList.erase(id);
            }
        }

        // has_piece_at(i, j, k)
        inline bool has_piece_at(int i, int j, int k) {
            return m_pieceList.count(k*(m_dimX*m_dimY)+j*m_dimY + i) > 0;
        }

        inline bool has_piece_at(int index) {
            return m_pieceList.count(index) > 0;
        }

        inline std::map<unsigned int, PuzzlePiece*>::iterator get_piece_at(unsigned int id) {
            // the iterator is the object of key->value that the map actually stores
            // the value is thus found by get_piece_at(ID_HERE)->second
            return m_pieceList.find(id);
        }

        // inline PuzzlePiece* get_piece_at(unsigned int id) {
        //     // the iterator is the object of key->value that the map actually stores
        //     // the value is thus found by get_piece_at(ID_HERE)->second
        //     return m_pieceList.find(id)->second;
        // }

        inline void check_contiguous(unsigned int gridX, unsigned int gridY, unsigned int gridZ) {
            // loop through all pieces,
            std::vector<PuzzlePiece> newPieces;
            std::vector<PuzzlePiece> intermediatePieces;

            for (std::map<unsigned int, PuzzlePiece*>::iterator p = m_pieceList.begin(); p != m_pieceList.end(); p++ ) {
                intermediatePieces = p->second->check_contiguous(gridX, gridY, gridZ);
                for (int i = 0; i < intermediatePieces.size(); i++) {
                    newPieces.push_back(intermediatePieces[i]);
                }

            }
            // update m_pieceList
            std::cout << "there are " << newPieces.size() << " pieces" << std::endl;
            m_pieceList.clear();
            for (int i = 0; i < newPieces.size(); i++) {
                // let it be known this is a bad idea
                // copy the fields from the stack piece to the heap piece
                PuzzlePiece newPiece = newPieces[i];
                PuzzlePiece* piece = new PuzzlePiece(newPiece.getID());
                piece->m_voxels = newPiece.m_voxels;
                piece->m_id = newPiece.m_id;

                // m_pieceList.insert(std::make_pair<unsigned int, PuzzlePiece>(newPieces[i].getID(), *newPieces[i]));
                m_pieceList[newPieces[i].getID()] = piece;
            }
        }

        // merge_pieces(id1, id2)

        unsigned int m_dimX, m_dimY, m_dimZ, m_size;

        // key is id. value is the piece itself
        std::map<unsigned int, PuzzlePiece*> m_pieceList;
        
    } Puzzle;

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

        inline bool & isOnSurface(int index) {
            return m_surfaceArray[index];
        }

        inline void setOnSurface(unsigned int i, unsigned int j, unsigned int k) {
            m_surfaceArray[k*(m_dimX*m_dimY)+j*m_dimY + i] = true;
        }  

        inline void setOnSurface(int index) {
            m_surfaceArray[index] = true;
        }

        inline unsigned int & getPieceNum(unsigned int i, unsigned int j, unsigned int k) {
            return m_pieceNumArray[k*(m_dimX*m_dimY)+j*m_dimY + i];
        }

        inline unsigned int & getPieceNum(int index) {
            return m_pieceNumArray[index];
        }

        inline void setPieceNum(unsigned int i, unsigned int j, unsigned int k, unsigned int pieceNum) {
            m_pieceNumArray[k*(m_dimX*m_dimY)+j*m_dimY + i] = pieceNum;
        }  

        inline void setPieceNum(int index, unsigned int pieceNum) {
            m_pieceNumArray[index] = pieceNum;
        }

        inline int getArrayIndexByCoordinate(unsigned int i, unsigned int j, unsigned int k) {
            return k*(m_dimX*m_dimY)+j*m_dimY + i;
        }

        // getting another voxel index via translation
        inline int getArrayIndexFromAnotherIndexViaTranslation(int index, unsigned int i, unsigned int j, unsigned int k) {
            return index + getArrayIndexByCoordinate(i,j,k);
        }

        // // getting another voxel index via translation in the x direction (i)
        // inline int getArrayIndexFromAnotherIndexViaTranslation_x(int index, unsigned int x_offset) {
        //     return index + x_offset;
        // }

        // // getting another voxel index via translation in the y direction (j)
        // inline int getArrayIndexFromAnotherIndexViaTranslation_y(int index, unsigned int y_offset) {
        //     return index + y_offset*m_dimY;
        // }

        // // getting another voxel index via translation in the z direction (k)
        // inline int getArrayIndexFromAnotherIndexViaTranslation_z(int index, unsigned int z_offset) {
        //     return index + z_offset*(m_dimX*m_dimY);
        // }

        // get one of the 6 neighbors
        // index is the index of the voxel we're neighboring from
        // offset is in range [0,5]
        inline int getNeighborVoxelIndexFromOffset(int index, int offset) {
            // check its neighbors 
            int neighbor_voxel_index;
            switch (offset) {
                case 0: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, -1, -1); break; 
                case 1: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, -1, 0); break; 
                case 2: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, -1, 1); break; 
                case 3: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, 0, -1); break; 
                case 4: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, 0, 0); break; 
                case 5: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, 0, 1); break; 
                case 6: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, 1, -1); break; 
                case 7: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, 1, 0); break; 
                case 8: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, -1, 1, 1); break; 
                case 9: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 0, -1, -1); break; 
                case 10: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 0, -1, 0); break; 
                case 11: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 0, -1, 1); break; 
                case 12: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 0, 0, -1); break; 
                case 13: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 0, 0, 1); break; 
                case 14: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 0, 1, -1); break; 
                case 15: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 0, 1, 0); break; 
                case 16: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 0, 1, 1); break; 
                case 17: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, -1, -1); break; 
                case 18: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, -1, 0); break; 
                case 19: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, -1, 1); break; 
                case 20: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, 0, -1); break; 
                case 21: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, 0, 0); break; 
                case 22: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, 0, 1); break; 
                case 23: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, 1, -1); break; 
                case 24: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, 1, 0); break; 
                case 25: 
                    neighbor_voxel_index = getArrayIndexFromAnotherIndexViaTranslation(index, 1, 1, 1); break; 
            }

            // check bounds
            if (neighbor_voxel_index < 0 || neighbor_voxel_index >= m_size) {
                // out of bounds, return bad
                return -1;
            }

            return neighbor_voxel_index;
        }

        inline Vec3* getVoxelLocationByIndex(int index) {
            return m_voxelToLocationMapping[index];
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
        
        inline void updatePiecesFromPuzzle(PuzzleStruct* puzzle) {
        // should change pieceNumArray to be the same as the puzzle structure.
            for (int i = 0; i < m_size; i++) {
                m_pieceNumArray[i] = 0;
            }

            // for each constituent voxel in the piece,
            // reset its piece number to match what's in the puzzle struct
            for (std::map<unsigned int, PuzzlePiece*>::iterator p = puzzle->m_pieceList.begin(); p != puzzle->m_pieceList.end(); p++ ) {
                PuzzlePiece* piece = p->second;
                std::vector<int>* constituentVoxelIdsFromPiece = piece->getVoxels();

                for (int i = 0; i < constituentVoxelIdsFromPiece->size(); i++) {
                    int voxelID = (constituentVoxelIdsFromPiece)->at(i);

                    int pieceID = piece->getID();

                    m_pieceNumArray[voxelID] = pieceID;
                }
            }
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

        // voxel -> i,j,k mapping
        Vec3* *m_voxelToLocationMapping;
        
    } VoxelGrid;

}

void loadVoxelGrid(const char * filename, CompFab::VoxelGrid & voxelGrid);
void saveVoxelGrid(const char * filename, CompFab::VoxelGrid & voxelGrid);

// if the VOX file is present, return true
bool isVoxelCacheFilePresent(const char * filename);

#endif
