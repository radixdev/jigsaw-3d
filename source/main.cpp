#include <iostream>
#include <vector>
#include <algorithm>    // std::min
#include <string>
#include <sstream>
#include <stdlib.h>
#include <ctime>
#include "../include/CompFab.h"
#include "../include/Mesh.h"
        
//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;
CompFab::Puzzle *g_puzzle;

double getRandomColor() {
    int resolution =  100;
    return ((double)(rand() % resolution)) / ((double) resolution);
}

std::string getFileSuffixFromArgs(char **argv) {
    return std::string(argv[2]).substr(0, std::string(argv[2]).size()-4);
}

// //Ray-Triangle Intersection
// //Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle, 
     * 0 otherwise */

    // first find if the ray intersections the plane.
    CompFab::Vec3 e1 = triangle.m_v2 - triangle.m_v1;
    CompFab::Vec3 e2 = triangle.m_v3 - triangle.m_v1;

    // dot normal of triangle with ray direction
    CompFab::Vec3 triangle_norm = e1 % e2;
    double check = triangle_norm * ray.m_direction;
    if (std::abs(check) < 1e-6 ) {
        return 0; //this means it is perpendicular to the normal and will not intersect
    }

    // plane of triangle is P dot N + d = 0
    double d = -1.0 * (triangle.m_v2 * triangle_norm);

    // ray is P = origin + t*dir
    // therefor (origin + t*dir) dot N + d = 0
    double t = -1.0 * (ray.m_origin * triangle_norm + d) / (ray.m_direction * triangle_norm);
    if (t < 0) {
        return 0;
    }

    // point of intersection between ray and plane of triangle
    CompFab::Vec3 point = ray.m_origin + CompFab::Vec3(ray.m_direction[0]*t, ray.m_direction[1]*t, ray.m_direction[2]*t);

    // use barycentric coordinates http://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    e1 = triangle.m_v2 - triangle.m_v1;
    e2 = triangle.m_v3 - triangle.m_v1;
    CompFab::Vec3 e3 = point - triangle.m_v1;

    double dot11 = e1*e1;
    double dot12 = e1*e2;
    double dot22 = e2*e2;
    double dot31 = e3*e1;
    double dot32 = e3*e2;

    double denom = 1.0/(dot11 * dot22 - dot12 * dot12);
    double beta = (dot22 * dot31 - dot12 * dot32) *denom;
    double gamma = (dot11 * dot32 - dot12 * dot31) *denom; 
    
    if (beta >= 0 && gamma >= 0 && gamma+beta <1) {
        return 1;
    } else {
        return 0;
    }

}

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    
    unsigned int numHits = 0;
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */
    for (int i = 0; i < g_triangleList.size(); i++) {
        CompFab::Ray r(voxelPos, dir);
        numHits += rayTriangleIntersection(r, g_triangleList[i]);
    }
    return numHits;
}
    
bool loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    
    return true;
   
}


/**
saveWithColors -> save to the outfile using the voxel's set color parameter
**/
void saveVoxelsToObj(const char * outfile, const bool checkSurface, const bool saveWithColors) {
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }

                if(checkSurface && !g_voxelGrid -> isOnSurface(ii, jj, kk)) {
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                if (saveWithColors) {
                    makeCube(box, box0, box1, g_voxelGrid->getVoxelColor_r(ii,jj,kk), g_voxelGrid->getVoxelColor_g(ii,jj,kk), g_voxelGrid->getVoxelColor_b(ii,jj,kk));
                } else {
                    // default drab grey if not saving with colors
                    makeCube(box, box0, box1, 0.5, 0.5, 0.5);
                }
                mout.append(box);
            }
        }
    }

    std::cout << "saving: " << outfile << std::endl;
    mout.save_obj(outfile);
}

void initializeVoxelGrid(char **argv, unsigned int dim) {
    // if the VOX file is present, read from it.
    // else, do the actual voxelization steps

    std::stringstream ss;
    ss<<dim;

    std::string voxelSaveFile = std::string(argv[1]).substr(0, std::string(argv[1]).size()-4) +"_dim"+ss.str()+ "_VOXELS.vox";

    if (isVoxelCacheFilePresent(voxelSaveFile.c_str())) {
        std::cout << "reading from stored VOX file" << std::endl;
        loadVoxelGrid(voxelSaveFile.c_str(), *g_voxelGrid);
    } else {
        std::cout << "generating voxels and VOX file" << std::endl;
        //Cast ray, check if voxel is inside or outside 
        //even number of surface intersections = outside (OUT then IN then OUT)
        // odd number = inside (IN then OUT)
        CompFab::Vec3 voxelPos;
        CompFab::Vec3 direction(1.0,0.0,0.0);
        CompFab::Vec3 direction2(0.0, 1.0, 0.0);
        CompFab::Vec3 direction3(0.0, 0.0, 1.0);

        std::clock_t start;
        start = std::clock();

        if (argv[3] && std::string(argv[3]) == "--multidir") {
            for (int k = 0; k < g_voxelGrid -> m_dimZ; k++) {
                for (int j = 0; j < g_voxelGrid -> m_dimY; j++) {
                    for (int i = 0; i < g_voxelGrid -> m_dimX; i++) {
                        voxelPos = CompFab::Vec3(g_voxelGrid -> m_lowerLeft[0] + g_voxelGrid -> m_spacing*i,
                                                 g_voxelGrid -> m_lowerLeft[1] + g_voxelGrid -> m_spacing*j,
                                                 g_voxelGrid -> m_lowerLeft[2] + g_voxelGrid -> m_spacing*k);
                        int intersectionNum = numSurfaceIntersections(voxelPos, direction);
                        int intersectionNum2 = numSurfaceIntersections(voxelPos, direction2);
                        int intersectionNum3 = numSurfaceIntersections(voxelPos, direction3);

                        if (intersectionNum % 2 == 1 && intersectionNum2 % 2 == 1 && intersectionNum3 % 2 == 1) {
                            // odd -> inside
                            g_voxelGrid->setIsInside(i,j,k);
                        }
                    }
                }
            }

        } else {
            /********* ASSIGNMENT ********
            /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
             * surface defined by the triangles in g_triangleList */

            for (int k = 0; k < g_voxelGrid -> m_dimZ; k++) {
                for (int j = 0; j < g_voxelGrid -> m_dimY; j++) {
                    for (int i = 0; i < g_voxelGrid -> m_dimX; i++) {
                        voxelPos = CompFab::Vec3(g_voxelGrid -> m_lowerLeft[0] + g_voxelGrid -> m_spacing*i,
                                                 g_voxelGrid -> m_lowerLeft[1] + g_voxelGrid -> m_spacing*j,
                                                 g_voxelGrid -> m_lowerLeft[2] + g_voxelGrid -> m_spacing*k);
                        int intersectionNum = numSurfaceIntersections(voxelPos, direction);
                        if (intersectionNum % 2 == 1) {
                            // odd -> inside
                            g_voxelGrid->setIsInside(i,j,k);
                        }
                    }
                }

                std::cout << "progress: " <<k << " of " << (g_voxelGrid->m_dimZ - 1) << std::endl;
            }
        }
        
        std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

        std::cout << "finished initial loop" << std::endl;

        //Write out voxel data as obj
        saveVoxelsToObj(argv[2], false, false);

        // save voxel data as VOX file
        saveVoxelGrid(voxelSaveFile.c_str(), *g_voxelGrid);
    }
}

void initializePuzzle(unsigned int dim) {
    g_puzzle = new CompFab::Puzzle(dim,dim,dim);
}

void setVoxelColors() {
    for (int k = 0; k < g_voxelGrid -> m_dimZ; k++) {
        for (int j = 0; j < g_voxelGrid -> m_dimY; j++) {
            for (int i = 0; i < g_voxelGrid -> m_dimX; i++) {

                // setting the voxel color here
                unsigned int pieceNum = g_voxelGrid->getPieceNum(i,j,k);
                srand(pieceNum);

                g_voxelGrid->setVoxelColor(i,j,k, getRandomColor(),getRandomColor(),getRandomColor());
            }
        }
    }
}

void doSurfaceAnalysis() {
    bool top, bottom, left, right, fwd, bck;
    
    for (int k = 0; k < g_voxelGrid -> m_dimZ; k++) {
        for (int j = 0; j < g_voxelGrid -> m_dimY; j++) {
            for (int i = 0; i < g_voxelGrid -> m_dimX; i++) {
                // only do this if we have a solid piece, not air
                if (g_voxelGrid->isInside(i, j, k)) {
                    left = g_voxelGrid->isInside(std::max(0, i-1), j, k);
                    right = g_voxelGrid->isInside(std::min((int) g_voxelGrid -> m_dimX -1, i+1), j, k);
                    top = g_voxelGrid->isInside(i, std::max(0, j-1), k);
                    bottom = g_voxelGrid->isInside(i, std::min((int) g_voxelGrid -> m_dimY -1, j+1), k);
                    fwd = g_voxelGrid->isInside(i, j, std::max(0, k-1));
                    bck = g_voxelGrid->isInside(i, j, std::min((int) g_voxelGrid -> m_dimZ -1, k+1));

                    if (!(left && right && top && bottom && fwd && bck)) {
                        g_voxelGrid->setOnSurface(i,j,k);
                    }

                } 
            }
        }
    }

    return;

    // do another layer for structural integrity
    std::vector<int> secondLayer;
    for (int k = 0; k < g_voxelGrid -> m_dimZ; k++) {
        for (int j = 0; j < g_voxelGrid -> m_dimY; j++) {
            for (int i = 0; i < g_voxelGrid -> m_dimX; i++) {
                // only do this if we have a solid piece, not air, but also not on the surface already
                if (g_voxelGrid->isInside(i, j, k) && !g_voxelGrid->isOnSurface(i, j, k)) {
                    left = g_voxelGrid->isOnSurface(std::max(0, i-1), j, k);
                    right = g_voxelGrid->isOnSurface(std::min((int) g_voxelGrid -> m_dimX -1, i+1), j, k);
                    top = g_voxelGrid->isOnSurface(i, std::max(0, j-1), k);
                    bottom = g_voxelGrid->isOnSurface(i, std::min((int) g_voxelGrid -> m_dimY -1, j+1), k);
                    fwd = g_voxelGrid->isOnSurface(i, j, std::max(0, k-1));
                    bck = g_voxelGrid->isOnSurface(i, j, std::min((int) g_voxelGrid -> m_dimZ -1, k+1));

                    if (left || right || top || bottom || fwd || bck) {
                        secondLayer.push_back(k*(g_voxelGrid->m_dimX*g_voxelGrid->m_dimY)+j*g_voxelGrid->m_dimY + i);
                    }

                } 
            }
        }
    }

    std::cout << "second layer has size: " << secondLayer.size() << std::endl;

    for (int i = 0; i < secondLayer.size(); i++ ) {
        g_voxelGrid->setOnSurface(secondLayer[i]);
    }
}

// int main(int argc, char **argv) {
//     std::map<unsigned int, CompFab::PuzzlePiece*> m;

//     m[1] = new CompFab::PuzzlePiece(1);

//     int a = 1/0;
//     return 0;
// }

int main(int argc, char **argv) {
    unsigned int dim = 32; //dimension of voxel grid (e.g. 32x32x32)

    //Load OBJ
    if(argc < 3)
    {
        std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        return 0;
    }
    
    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    loadMesh(argv[1], dim); 
    
    unsigned int PIECE_SIZE = 10;

    std::cout<<"Creating Puzzle : "<<"\n";
    initializePuzzle(dim/PIECE_SIZE);


    std::clock_t start;
    start = std::clock();

    /////////////
    /////////////
    /////////////

    initializeVoxelGrid(argv, dim);
    std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

    std::cout << "starting surface analysis" << std::endl;

    doSurfaceAnalysis();

    std::cout << "finished surface analysis" << std::endl;

    std::string surfacefile = getFileSuffixFromArgs(argv) + "_surface.obj";
    std::cout << surfacefile << std::endl;

    saveVoxelsToObj(surfacefile.c_str(), true, false);

    /////////////
    /////////////
    /////////////

    std::cout << "starting divisions" << std::endl;
    unsigned int newx, newy, newz, pieceNum;
    for (int k = 0; k < g_voxelGrid -> m_dimZ; k++) {
        for (int j = 0; j < g_voxelGrid -> m_dimY; j++) {
            for (int i = 0; i < g_voxelGrid -> m_dimX; i++) {
                // only do this if we have a solid piece, not air
                if (g_voxelGrid->isOnSurface(i, j, k)) {
                    newx = i/PIECE_SIZE;
                    newy = j/PIECE_SIZE;
                    newz = k/PIECE_SIZE;
                    pieceNum = newz*(g_voxelGrid -> m_dimX/PIECE_SIZE*g_voxelGrid ->m_dimY/PIECE_SIZE)+newy*g_voxelGrid ->m_dimY/PIECE_SIZE + newx;
                    g_voxelGrid -> setPieceNum(i, j, k, pieceNum);

                    // create the piece in the puzzle struct
                    if (!g_puzzle->has_piece_at(i,j,k)) {
                        g_puzzle->add_piece(pieceNum);
                    }

                    // add this voxel to that piece
                    // std::vector<int> newVoxels;
                    // newVoxels.push_back(k*(g_voxelGrid -> m_dimX*g_voxelGrid ->m_dimY)+j*g_voxelGrid ->m_dimY + i);
                    // g_puzzle->get_piece_at(pieceNum)->second->add_voxels(newVoxels);

                    // add a single voxel to the piece
                    int new_voxel_id = k*(g_voxelGrid -> m_dimX*g_voxelGrid ->m_dimY)+j*g_voxelGrid ->m_dimY + i;
                    g_puzzle->get_piece_at(pieceNum)->second->add_voxel(new_voxel_id);
                } 
            }
        }
    }

    g_puzzle->check_contiguous(g_voxelGrid -> m_dimX, g_voxelGrid -> m_dimY, g_voxelGrid -> m_dimZ);
    std::cout << "updating the voxel grid" << std::endl;
    g_voxelGrid->updatePiecesFromPuzzle(g_puzzle);

    /////////////
    /////////////
    /////////////

    std::cout << "merging small parts together" << std::endl;
    // iterate over all puzzle pieces
    // if it has too small size, merge with the piece with the largest neighbors

    // std::vector<CompFab::PuzzlePiece> allPieces = g_puzzle->get_pieces();
    // std::vector<CompFab::PuzzlePiece>::iterator itr;
    // for (itr = allPieces.begin(); itr < allPieces.end(); ++itr) {
    //     // CompFab::PuzzlePiece piece = (CompFab::PuzzlePiece) *itr;

    //     // if (piece.size() < 10) {
    //     //     // its too small
    //     //     // std::cout << "found too small piece" << piece.size() << std::endl;
    //     // }
    // }

    // for (CompFab::PuzzlePiece piece : allPieces) {

    // }


    std::cout << "done merging" << std::endl;
    /////////////
    /////////////
    /////////////

    setVoxelColors();

    int surfaces = 0;
    int pieces = 0;

    for (int i = 0; i < g_voxelGrid -> m_size; i++ ){
        if (g_voxelGrid->isOnSurface(i)) {
            surfaces += 1;
        }
        if (g_voxelGrid->getPieceNum(i) != 0) {
            pieces += 1;
        } 
    }

    std::cout<<"found " <<surfaces<<" surface voxels and " << pieces << " piece voxels" << std::endl;

    std::cout << "ending divisions" << std::endl;
    std::string dividedfile = getFileSuffixFromArgs(argv) + "_divided.obj";
    std::cout << dividedfile << std::endl;

    saveVoxelsToObj(dividedfile.c_str(), true, true);
    delete g_voxelGrid;
    delete g_puzzle;
}