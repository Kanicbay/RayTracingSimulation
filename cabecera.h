#include "math.h"
#include <fstream>
#include <iostream>
using namespace std;
const double PI = M_PI;
const int xy = 0;
const int xz = 1;
const int yz = 2;
const int MaxReflections = 51;
const int numReceptors = 1;
int numRoomTriangles = 1;
double maxDistanceRoom = 0.0f;


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// VECTOR CLASS

class vector{
public:
    double x, y, z;         // coordinates
    double m;               // module

    vector(){               // constructor
        x = 0.0;
        y = 0.0;
        z = 0.0;
    };

    void operator=(double num){     // assing same value to all the coordinates
        x = y = z = num;
    };

    double Module(){                // find vector module -> double
        m = sqrt(x*x + y*y + z*z);
        if(m == 0.0f){
            return 0.0f;
        }
        return m;
    };

    vector operator/(double num){   // overload vector / double -> vector
        vector v1;
        v1.x = x / num;
        v1.y = y / num;
        v1.z = z / num;
        return v1;
    };

    vector operator*(double num){   // overload vector * double -> vector
        vector v1;
        v1.x = x * num;
        v1.y = y * num;
        v1.z = z * num;
        return v1;
    };

    vector operator/(vector v){     // overload producto vectorial -> vector
        vector v1;
        v1.x =  y * v.z - z * v.y;
        v1.y = -x * v.z + z * v.x;
        v1.z =  x * v.y - y * v.x;
        return v1;
    };

    vector operator-(vector v){     // overload vector - vector -> vector
        vector v1;
        v1.x = x - v.x;
        v1.y = y - v.y;
        v1.z = z - v.z;
        return v1;
    };

    vector operator+(double num){   // overload vector + double -> vector
        vector v1;
        v1.x = x + num;
        v1.y = y + num;
        v1.z = z + num;
        return v1;
    }

    double operator*(vector v){     // overload producto escalar -> double
        return x*v.x + y*v.y + z*v.z;
    };

    vector unitar(){
        vector v;
        v.x = x;
        v.y = y;
        v.z = z;
        if(v.Module() == 0){
            v = 0;
            return v;
        }
        v = v / v.Module();
        return v;
    };

    void printVector(){
        cout << "Vector:  x: " << x << " , y: " << y << " , z: " << z << endl;
    }

    void clear(){     // Clear coordinates
        x = 0;
        y = 0;
        z = 0;
    };

    bool operator==(vector v1){
        if(v1.x == x && v1.y == y && v1.z == z){
             return true;
        }
        return false;
    };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// POINT CLASS

class point{
public:
    double x, y, z;

    point(){
        x = 0.0;
        y = 0.0;
        z = 0.0;
    };

    void operator=(int num){     // assing same value to all the coordinates
        x = num;
        y = num;
        z = num;
    };

    point operator+(point p){   // overload point + point -> point
        point p1;
        p1.x = p.x + x;
        p1.y = p.y + y;
        p1.z = p.z + z;
        return p1;
    };

    point operator/(double num){    // overload point / double -> point
        point p1;
        p1.x = x / num;
        p1.y = y / num;
        p1.z = z / num;
        return p1;
    };

    vector operator-(point p){      // overload point - vector -> vector
        vector v1;
        v1.x = x - p.x;
        v1.y = y - p.y;
        v1.z = z - p.z;
        return v1;
    };

    point operator+(vector v){      // overload point + vector -> point
        point p1;
        p1.x = x + v.x;
        p1.y = y + v.y;
        p1.z = z + v.z;
        return p1;
    };

    vector operator*(vector v){     // overload point * point -> point
        vector v1;
        v1.x = x * v.x;
        v1.y = y * v.y;
        v1.z = z * v.z;
        return v1;
    };

    point operator+(double num){    // overload point + num -> point
        point p1;
        p1.x = x + num;
        p1.y = y + num;
        p1.z = z + num;
        return p1;
    };

    point operator=(vector v){
        point p1;
        p1.x = v.x;
        p1.y = v.y;
        p1.z = v.z;
        return p1;
    }

    double distance (point p){      // distance between points -> double
        return sqrt((p.x - x)*(p.x - x) + (p.y - y)*(p.y - y) + (p.z - z)*(p.z - z));
    };

    void clear(){     // Clear coordinates
        x = 0;
        y = 0;
        z = 0;
    };

    void printPoint(){
        cout << x << " ," << y << " ," << z << endl;
    }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SOURCE CLASS

class source{
public:
    point sourcePosition;  // Position source also rays origin point
    vector *Rays;          // Rays vector (Direction)
    int NRAYS;             // Num of rays
    double eS;             // Initial Energy for the source, it needs to be distributed for all rays (es/NRays)

    source(){               // Constructor
        sourcePosition = 0;
        Rays = nullptr;
        NRAYS = 0;
        eS = 0.0f;
    }

    void clear() {
        delete[] Rays;
        sourcePosition = 0;
        NRAYS = 0;
        eS = 0.0f;
    };

    void createRays(int NumberOfRays) {       // Create rays vector based on icosahedron polygon (Code provide by the professor)
        NRAYS = NumberOfRays;

        //matriz das Arestas {1o ponto da aresta, 2o ponto da aresta, Posição dos pontos da aresta na matriz Rays}
        int A[30][3]= {{0,1,0}, {0,2,0}, {0,3,0}, {0,4,0}, {0,5,0},
                        {1,6,0}, {2,6,0}, {2,7,0}, {3,7,0}, {3,8,0},
                        {4,8,0}, {4,9,0}, {5,9,0}, {5,10,0},{1,10,0},
                        {6,11,0},{7,11,0},{8,11,0},{9,11,0},{10,11,0},
                        {1,2,0}, {2,3,0}, {3,4,0}, {4,5,0}, {5,1,0},
                        {6,7,0}, {7,8,0}, {8,9,0}, {9,10,0},{10,6,0}
                      };

        //matriz dos triangulos {1a aresta, 2a aresta, [0] V em pé [-1] V de cabeça pra baixo}
        int T[20][3]= {{0,1,0},   {1,2,0},   {2,3,0},   {3,4,0},   {4,0,0},
                        {5,6,-1},  {6,7,0},   {7,8,-1},  {8,9,0},   {9,10,-1},
                        {10,11,0}, {11,12,-1},{12,13,0}, {13,14,-1},{14,5,0},
                        {15,16,-1},{16,17,-1},{17,18,-1},{18,19,-1},{19,15,-1}
                      };

        int i,j,k,n,m,RAY;

        double S,R,xB,yB,zB,xC,yC,zC,c[8];
        //create Rays matrix

        if(NRAYS > 0)
            delete[] Rays;

        n = int(floor(sqrt((NRAYS-2)/10)+0.5));

        NRAYS = int(2+10*pow(n,2));

        Rays = new vector[NRAYS];

        //calculating the icosaedron vertives
        S = 2/sqrt(5);
        R = (5-sqrt(5))/5;
        Rays[0].x = 0;
        Rays[0].y = 0;
        Rays[0].z = 1;
        for(i=1; i<6; i++) {

            Rays[i].x = S*cos((PI*i*72)/180);
            Rays[i].y = S*sin((PI*i*72)/180);
            Rays[i].z = 1-R;
            Rays[i+5].x = S*cos((72*PI*i)/180+(36*PI)/180);
            Rays[i+5].y = S*sin((72*PI*i)/180+(36*PI)/180);
            Rays[i+5].z = R-1;
        }
        Rays[11].x = 0;
        Rays[11].y = 0;
        Rays[11].z = -1;
        RAY = 12;

        //calculating the rays on the icosaedron edges
        for(j=0; j<30; j++) {
            A[j][2] = RAY;
            xB = Rays[A[j][0]].x;
            yB = Rays[A[j][0]].y;
            zB = Rays[A[j][0]].z;
            xC = Rays[A[j][1]].x;
            yC = Rays[A[j][1]].y;
            zC = Rays[A[j][1]].z;
            c[0] = pow(xC, 2) * (pow(yB, 2) + pow(zB, 2)) + pow(yC * zB - yB * zC, 2) - 2 * xB * xC * (yB * yC + zB * zC) + pow(xB, 2) * (pow(yC, 2) + pow(zC, 2));

            c[1] = acos(xB * xC + yB * yC + zB * zC);
            c[2] = -xC * (yB * yC + zB * zC) + xB * (pow(yC, 2) + pow(zC, 2));
            c[3] = xC * (pow(yB, 2) + pow(zB, 2)) - xB * (yB * yC + zB * zC);
            c[4] = pow(xC, 2) * yB - xB * xC * yC + zC * (-yC * zB + yB * zC);
            c[5] = -xB * xC * yB + pow(xB, 2) * yC + zB * (yC * zB - yB * zC);
            c[6] = pow(xC, 2) * zB - xB * xC * zC + yC * (yC * zB - yB * zC);
            c[7] = -xB * xC * zB + pow(xB, 2) * zC + yB * (-yC * zB + yB * zC);

            for(i=1; i<n; i++) {
                Rays[RAY].x = (c[2] * cos(i * c[1] / n) + c[3] * cos((n - i) * c[1] / n)) / c[0];
                Rays[RAY].y = (c[4] * cos(i * c[1] / n) + c[5] * cos((n - i) * c[1] / n)) / c[0];
                Rays[RAY].z = (c[6] * cos(i * c[1] / n) + c[7] * cos((n - i) * c[1] / n)) / c[0];
                RAY++;
            }
        }

        //calculating the rays on the icosaedron faces

        for(k=0; k<20; k++)
            for(j=1; j<n; j++) {
                xB = Rays[A[T[k][0]][2]+j-1].x;
                yB = Rays[A[T[k][0]][2]+j-1].y;
                zB = Rays[A[T[k][0]][2]+j-1].z;
                xC = Rays[A[T[k][1]][2]+j-1].x;
                yC = Rays[A[T[k][1]][2]+j-1].y;
                zC = Rays[A[T[k][1]][2]+j-1].z;
                c[0] = pow(xC, 2) * (pow(yB, 2) + pow(zB, 2)) + pow(yC * zB - yB * zC, 2) - 2 * xB * xC * (yB * yC + zB * zC) + pow(xB, 2) * (pow(yC, 2) + pow(zC, 2));

                c[1] = acos(xB * xC + yB * yC + zB * zC);
                c[2] = -xC * (yB * yC + zB * zC) + xB * (pow(yC, 2) + pow(zC, 2));
                c[3] = xC * (pow(yB, 2) + pow(zB, 2)) - xB * (yB * yC + zB * zC);
                c[4] = pow(xC, 2) * yB - xB * xC * yC + zC * (-yC * zB + yB * zC);
                c[5] = -xB * xC * yB + pow(xB, 2) * yC + zB * (yC * zB - yB * zC);
                c[6] = pow(xC, 2) * zB - xB * xC * zC + yC * (yC * zB - yB * zC);
                c[7] = -xB * xC * zB + pow(xB, 2) * zC + yB * (-yC * zB + yB * zC);

                if(T[k][2]==0)
                    m=j;
                else
                    m=n-j;

                for(i=1; i<m; i++) {
                    Rays[RAY].x = (c[2] * cos(i * c[1] / m) + c[3] * cos((m - i) * c[1] / m)) / c[0];
                    Rays[RAY].y = (c[4] * cos(i * c[1] / m) + c[5] * cos((m - i) * c[1] / m)) / c[0];
                    Rays[RAY].z = (c[6] * cos(i * c[1] / m) + c[7] * cos((m - i) * c[1] / m)) / c[0];
                    RAY++;
                }
            }
    };

};


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// TRIANGLE CLASS       (Provided by the professor)

class triangle {
public:
    point p0,p1,p2,bc;  //triangle Points
    int Projection;     //projection
    double a0;          //a0 constante para cálculos futuros
    int ID;             //identificador único

    triangle(){
        p0=0;
        p1=0;
        p2=0;
        bc=0;
        Projection=0;
        a0=0;
        ID=0;
    };

    void operator=(triangle t){
        p0=t.p0;
        p1=t.p1;
        p2=t.p2;
        bc=t.bc;
        Projection=t.Projection;
        a0=t.a0;
        ID=t.ID;
    };

    void Clear(){
        p0=0;
        p1=0;
        p2=0;
        bc=0;
        Projection=0;
        a0=0;
        ID=0;
    };

    void Centroid(){
        bc=(p0+p1+p2)/3;
    };

    double solidAngle(point b){
        double area=0.0,d=0.2;
        triangle t;
        vector v0,v1,v2;
        v0=p0-b;
        v1=p1-b;
        v2=p2-b;
        v0=v0/v0.Module();
        v1=v1/v1.Module();
        v2=v2/v2.Module();
        t.p0=b+(v0*d);
        t.p1=b+(v1*d);
        t.p2=b+(v2*d);
        area=t.TriangleArea();
        return area;
    };

    double TriangleArea() {
        double a;
        vector v=(p1-p0)/(p2-p0);
        a=0.5*v.Module();
        return a;
    };

    void CalculateProjection() {
        vector n;
        double x0,y0,z0,x1,y1,z1,x2,y2,z2;
        x0=p0.x;
        y0=p0.y;
        z0=p0.z;
        x1=p1.x;
        y1=p1.y;
        z1=p1.z;
        x2=p2.x;
        y2=p2.y;
        z2=p2.z;
        n=(p1-p0)/(p2-p0);
        n.x=n.x*n.x;
        n.y=n.y*n.y;
        n.z=n.z*n.z;
        if((n.x>=n.y)&&(n.x>=n.z)) {                        //projeção yz
            Projection=yz;
            a0=1/(-y1*z0+y2*z0+y0*z1-y2*z1-y0*z2+y1*z2 + 0.000001);
        }
        if((n.y>=n.x)&&(n.y>=n.z)) {                        //projeção xz
            Projection=xz;
            a0=1/(-x1*z0+x2*z0+x0*z1-x2*z1-x0*z2+x1*z2 + 0.000001);
        }
        if((n.z>=n.x)&&(n.z>=n.y)) {                        //projeção xy
            Projection=xy;
            a0=1/(-x1*y0+x2*y0+x0*y1-x2*y1-x0*y2+x1*y2 + 0.000001);
        }
    };

    void baricentro() {
        bc = (p0 + p1 + p2) / 3;
    };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// TRIANGLE NORMAL FUNCTION

vector TriangleNormal(triangle t){          // Calculate triangle normal
    vector n;
    n = (t.p1 - t.p0) / (t.p2 - t.p0);
    if(n.Module() == 0){
        n = 0;
        return n;
    }
    n = n / n.Module();
    return n;
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// POINT GRID CLASS

class pointGrid {
public:
    point **e;   // malla[m][n]
    int m,n;

    pointGrid(){
        e = nullptr;
        m = 0;
        n = 0;
    };

    void initialize(int i, int j){
        m = i;
        n = j;
        e = new point*[m];
        for(int i=0; i<m; i++){
            e[i] = new point[n];
        }
    };
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// PLANE CLASS

class plane {
public:
    int Npoints, Ntriangles, NSubplanes;   // Number of points, number of triangles
    double distance;
    point *p;                  // point vector
    plane *sp;                 // sub plane vector
    triangle *t;               // triangle vector
    vector n;                  // normal vector
    pointGrid pg;

    plane() {                  // Plane constructor
        Npoints = 0;
        Ntriangles = 0;
        distance = 0.0f;
        NSubplanes = 0;

        t = nullptr;
        p = nullptr;
    };

    void clear() {              // Clear the instance
        Npoints = 0;
        Ntriangles = 0;
        delete[] p;
        delete[] t;
    };

    void NewPoints(int n){      // Add new points for the plane
        point *tp;                               // Temporal Point Vector
        tp = new point[Npoints + n];             // Init the vector with (Npoints + n) length
        int i;
        for(i = 0; i < Npoints; i++){            // If there is a previous points
            tp[i] = p[i];                        // copy everything from previous to the
        }                                        // temporal one
        for(i = Npoints; i < Npoints + n; i++){
            tp[i].clear();                       // Clear new points
        }

        if(Npoints > 0){                         // If there is a previos points so
            delete[] p;                          // clear the vector
            p = nullptr;
        }
        p = tp;                                  // Assign new values
        Npoints += n;                            // Update Npoints
    };

    void NewTriangle(int numTriangles){          // Create new triangles
        t = new triangle[numTriangles];          // update triangle vector
        Ntriangles = numTriangles;               // update num of triangles
    };

    void PointGetTriangle() {                    // Find triangle points from cube plane
        NewTriangle(Npoints-2);                  // create num of triangles base of (Number of points into the plane minus 2)
        int i = 1;
        for(int T=0; T<Ntriangles; T++){         // move over the triangles
            i--;                                 // asign i equals to 0 and find the first point
            t[T].p0.x = p[i].x;
            t[T].p0.y = p[i].y;
            t[T].p0.z = p[i].z;
            i++;
            if(i==Npoints) i=0;                  // assing 0 to i when i is equal to Npoints helps to get back to the triangle first point
            t[T].p1.x = p[i].x;
            t[T].p1.y = p[i].y;
            t[T].p1.z = p[i].z;
            i++;
            if(i==Npoints) i=0;
            t[T].p2.x = p[i].x;
            t[T].p2.y = p[i].y;
            t[T].p2.z = p[i].z;
            i++;
        }
    };

    vector *getPerpendicularVectors(){

        vector *tp2;                                        // vector of possible perpendicular vectors
        vector *perpendVector;                              // vector of perpendicular vectors
        int i,j = 0;

        // --------------Verify perpendicular vector--------------
        int numPPV = Npoints-1;                         // Num of possible perpendicular vectors establish to 3
        tp2 = new vector[numPPV];                       // vector of possible perpendicular vectors
        perpendVector = new vector[numPPV-1];           // vector of perpendicular vector establish to 2


        //--------------Save all the posible perpendicular vector--------------
        for(i=1; i<=numPPV; i++){
            tp2[i-1] = (p[i] - p[0]);                   // make vectors from the first point and save it
        }

        //--------------Find the vector two perpendicular vectors and save it as a unitar vectors--------------
        for(i=0; i<(Npoints-1); i++){
            for(j=0; j<(Npoints-1); j++){
                if(tp2[i]*tp2[j] == 0){
                    perpendVector[0] = tp2[j].unitar();             // Iter each possible vector over all the possible vectors
                    perpendVector[1] = tp2[i].unitar();             // And fin if the scalar product gives 0, so it's perpendicular
                    distance = tp2[j].Module();                     // save the two unitar vectors and break the for
                    break;
                }
            }
        }
        return perpendVector;

    }

    void genMoreTriangle(int nDiv){

        if(Npoints == 4){
            pointGrid  np;                                      // create point grid vector
            np.initialize(nDiv + 1, nDiv + 1);                  // init the point grid with nDiv*nDiv dimention
            plane *pt;                                          // create a temporal plane vector
            int npt;                                            // number of planes
            npt = nDiv * nDiv;                                  // number of planes equals to nDiv*nDiv
            pt = new plane[npt];                                // temporal plane of dimension of npt
            NewTriangle(npt*2);                                 // Init the triangles of the plane
            int i,j = 0;                                        // for indexes
            double steps=0;                                     // step variable, decide the number of step for each point
            Ntriangles = 0;                                     // establish triangles to 0
            sp  = new plane[npt];                               // vector of subplanes vectors
            NSubplanes = npt;                                   // update number of subplanes

            //--------------Find perpendicular vectors----------------------------
            vector *perpendVector = getPerpendicularVectors();

            //--------------create the sub-planes and their triangles-------------
            steps = distance/nDiv;                  // Find the steps for each point
            point p0, p1, p2, p3;                   // Point for each sub-plane
            point initPoint = p[0];                 // saves the first point of the first subple
            point initialPoint = p[0];                      // initial point set with the first plane point
            int contador=0;
            int contador2=0;
            int z=0;

            for(j=0; j<nDiv; j++){                  // Iter each column or row, depend of the direction
                for(i=0; i<nDiv; i++){              // Iter each column or row, depend of the direction

                    p0 = initialPoint;                          // calculate the sub-plane first point
                    p1 = p0 + ((perpendVector[0] * steps));     // calculate the sub-plane second point
                    p2 = p0 + ((perpendVector[1] * steps));     // calculate the sub-plane third point
                    p3 = p2 + ((perpendVector[0] * steps));     // calculate the sub-plane foruth point

                    initialPoint = p2;              // change the initialPoint to p2 to continue to the next sub-plane
                    if(i==0) initPoint = p1;        // save p1 to initPoint to continue in the other column or row with the next sub-plane

                    // Save the sub plane
                    pt[contador].NewPoints(4);      // init the sub-plane
                    pt[contador].p[0] = p0;         // save the first point
                    pt[contador].p[1] = p1;         // save the second point
                    pt[contador].p[2] = p3;         // save the third point
                    pt[contador].p[3] = p2;         // save the fourth point

                    // Save points into the grid
                    if(j==0 && i==0){
                        np.e[z][contador2] = p0;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                        np.e[z][contador2] = p2;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                        np.e[z][contador2] = p1;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                        np.e[z][contador2] = p3;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                    }
                    if(j==0 && i!=0){
                        np.e[z][contador2] = p2;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                        np.e[z][contador2] = p3;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                    }
                    if(j!=0 && i==0){
                        np.e[z][contador2] = p1;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                        np.e[z][contador2] = p3;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                    }
                    if(j!=0 && i!=0){
                        np.e[z][contador2] = p3;
                        contador2 = (contador2 >= nDiv) ? 0 : contador2+1; if(contador2==0) z++;
                    }

                    // Create two triangles into the sub plane
                    pt[contador].PointGetTriangle();    // create the two triangles inside the sub-plane


                    // Save the subtriangles into the triangle array
                    t[Ntriangles] = pt[contador].t[0];           //save the first triangle into triangle vector
                    t[Ntriangles+1] = pt[contador].t[1];         //save the second triangle into triangle vector
                    //t[Ntriangles].ID = (numRoomTriangles);         // assing its ID
                   // t[Ntriangles + 1].ID = (numRoomTriangles + 1);     // assing its ID
                    numRoomTriangles+=2;
                    Ntriangles += 2;                             // increments the NTriangles value
                    sp[contador] = pt[contador];
                    contador++;
                }
                initialPoint = initPoint;               // assing initialPoint to the initPoint
            }
            pg = np;
        }
    };

    void printTriangleInformation(){        // print triangle information
        cout << "Number of triangles: " << Ntriangles << endl;
        for(int T=0; T<Ntriangles; T++){
            cout << "Triangle " << (T+1) << ": " << endl;
            cout << "(" << t[T].p0.x << "  ," << t[T].p0.y << "  ," << t[T].p0.z << "), ID: " << t[T].ID << endl;
            cout << "(" << t[T].p1.x << "  ," << t[T].p1.y << "  ," << t[T].p1.z << "), ID: " << t[T].ID << endl;
            cout << "(" << t[T].p2.x << "  ," << t[T].p2.y << "  ," << t[T].p2.z << "), ID: " << t[T].ID << endl;
            cout << "Baricentro: ";
            t[T].bc.printPoint();

        }
        cout << endl;
    };

    void printTriangleInformationTxt(){
        for(int T=0; T<Ntriangles; T++){
            cout << "(" << t[T].p0.x << "," << t[T].p0.y << "," << t[T].p0.z << "), ";
            cout << "(" << t[T].p1.x << "," << t[T].p1.y << "," << t[T].p1.z << "), ";
            cout << "(" << t[T].p2.x << "," << t[T].p2.y << "," << t[T].p2.z << "), ";
            cout<<endl;
        }
    };

    void printSubPlaneInformationTxt(){
        for(int T=0; T<NSubplanes; T++){
            cout << "(" << sp[T].p[0].x << "," << sp[T].p[0].y << "," << sp[T].p[0].z << "), ";
            cout << "(" << sp[T].p[1].x << "," << sp[T].p[1].y << "," << sp[T].p[1].z << "), ";
            cout << "(" << sp[T].p[2].x << "," << sp[T].p[2].y << "," << sp[T].p[2].z << "), ";
            cout << "(" << sp[T].p[3].x << "," << sp[T].p[3].y << "," << sp[T].p[3].z << "), ";
            cout<<endl;
        }
    };

    void printPlaneGridInformationTxt(){
        for(int m=0; m<pg.m ; m++){
            for(int n=0; n<pg.n ; n++){
                cout << "(" << pg.e[m][n].x << "," << pg.e[m][n].y << "," << pg.e[m][n].z << "), ";
            }
        }
        cout<<endl;
    };

    void printNormalInformation(){  // print normal plane information
        cout << "Normal Plane Information:  x: " << n.x << "  y: " << n.y << "  z: " << n.z << endl;
    };

};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// REFLECTION STRUCTURE

struct reflection{
    point incidentPoint[MaxReflections];       // incident point vector
    double d[MaxReflections];                  // distance between points vector
    int idTriangle[MaxReflections];            // id of the collision triangle
    int idPlano[MaxReflections];               // id of the collision plane
    vector vec[MaxReflections];                // indicent vectors
    double mod[MaxReflections];                // indicent vectors distance
    double tim[MaxReflections];                // indicent vectors times

    int N;                                     // number reflection
    bool lost;                                 // lost ray
    int Ray;                                   // id of the ray / number of ray
    double actualEnergy[MaxReflections];       // actual energy vector
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CLASS Double vector
class DoubleVector{
public:
    double **vec;                               // main vector
    int n, m;                                   // main vector dimensions nxm

    DoubleVector(){
        vec = nullptr;
        n = m = 0;
    };

    void operator=(double num){                        // assing a num for all the vectors elements
        if(n != 0){
            for(int i=0; i<n; i++){
                for(int j=0; j<m; j++){
                    vec[i][j] = num;
                }
            }
        }
    };

    void operator=(DoubleVector v){                     // assing vector to another one
        if(n != 0){
            for(int i=0; i<n; i++){
                for(int j=0; j<m; j++){
                    vec[i][j] = v.vec[i][j];
                }
            }
        }
    };

    void initiVector(int row, int col){             // init the vector
        n = row;
        m = col;
        vec = new double*[n];
        for(int i=0; i<n; i++){
            vec[i] = new double[m];
        }
    }

    void printVector(){                             // print the vector
        if(n != 0){
            cout << "[" << endl;
            for(int i=0; i<n; i++){
                cout << "[" << endl;
                for(int j=0; j<m; j++){
                    cout << vec[i][j] << ", ";
                }
                cout << "]," << endl;
            }
            cout << "]";
        }
    }


    void saveData(){                                        // save the data
        std::ofstream ofs("Simulation_Data/Room_Data/MatrixTemporalEspace.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());

        cout << "[" << endl;
        for(int i=0; i<n; i++){
            cout << "[";
            for(int j=0; j<m; j++){
                if(j!=(m-1)){
                    cout << vec[i][j] << ", ";
                } else {
                    cout << vec[i][j];
                }
            }
            cout << "]," << endl;
        }
        cout << "]";
        std::cout.rdbuf(cout_buff);
    }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CLASS vector vector
class VectorV{
public:
    vector **vec;
    int n, m;

    VectorV(){
        vec = nullptr;
        n = m = 0;
    };

    void operator=(VectorV v){
        if(n != 0){
            for(int i=0; i<n; i++){
                for(int j=0; j<m; j++){
                    vec[i][j] = v.vec[i][j];
                }
            }
        }
    };

    void initiVector(int row, int col){
        n = row;
        m = col;
        vec = new vector*[n];
        for(int i=0; i<n; i++){
            vec[i] = new vector[m];
        }
    }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// RECEPTOR CLASS
class receptor {
public:
    point position;                 // receptor position
    double diskRadius;              // radius of the disk
    double *resultEnergy;           // result energy from the Matrix Temporal Espace
    int nEnergy;                    // num of data into the energy result
    double left=0.0f;

    receptor(){
        position = 0;
        diskRadius = 0;
        resultEnergy = nullptr;
        nEnergy=0;
    };

    void initEnergy(int n){
        resultEnergy = new double[n];
        nEnergy = n;
    };

    double solidAngle(point bc){
        vector v = position - bc;
        double r = (diskRadius*0.2) / v.Module();
        return PI*r*r;
    }

    void saveData(){
        string documentName = "Simulation_Data/Receptor_Data/data_" + to_string(position.x) + "_" + to_string(position.y) + "_" + to_string(position.z) + ".txt";
        std::ofstream ofs(documentName);
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());

        cout << "[ ";
        for(int i=0; i<nEnergy; i++){
            cout << resultEnergy[i] << ", ";
        }
        cout << " ]";
        std::cout.rdbuf(cout_buff);
    }

    // Calculate energy for the receptor
    void addEnergy(int index_reflection, int numRay, double V_SON, int tim, double Er, reflection *ry, double distance){
        point rpdis, rpci;
        vector rayUnitary, normalReceptor;

        rpdis = position;
        rpci = ry[numRay].incidentPoint[index_reflection-1];
        rayUnitary = ry[numRay].vec[index_reflection].unitar();
        normalReceptor = rayUnitary * -1;

        double distance_ray_receptor = normalReceptor * rayUnitary;

        if(distance_ray_receptor == 0)
            distance_ray_receptor = -1;
        else
            distance_ray_receptor = (normalReceptor*(rpdis-rpci))/distance_ray_receptor;

        int tim1 = 0;
        if(distance_ray_receptor>0 && distance_ray_receptor<maxDistanceRoom){
            rpdis = rpci + (rayUnitary*distance_ray_receptor);
            double distancia_incidente_receptor = (rpdis-position).Module();
            if(distancia_incidente_receptor<diskRadius){
                distance_ray_receptor = (rayUnitary*distance_ray_receptor).Module();
                tim1 = int(1000*((distance_ray_receptor)/V_SON)) + tim;
                //cout << "Entro con el rayo: " << numRay <<" y con el tiempo: " << tim1 << endl;
                resultEnergy[tim1] += Er;
            }
        }
    }
};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ROOM CLASS

class room {
public:
    plane *p;               // plane vector
    plane *tp;              // plane temporal vector
    int NP;                 // number of planes
    double maxDistance;     // max distance into the room, needs to verify if a ray is lost
    double sourceEnergy;    // source energy, needs to assign energy to each ray
    double alpha, deltha;   // variable to calculate energy
    int **timeFlying;                 // vector to find the flying time
    int **timeFlyingR;
    DoubleVector solidAngle;          // solid angle double vector
    DoubleVector solidAnglesReceptor;
    DoubleVector distanceTriangles;
    DoubleVector distanceReceptors;
    receptor *r;

    room(){
        p = nullptr;
        tp = nullptr;
        NP = 0;
        sourceEnergy = 0.0f;
    };

    void clear() {
        delete[] p;
        delete[] tp;
        NP = 0;
        maxDistance = 0;
        alpha = 0;
        deltha = 0;
    };

    void NewPlanes(int N){      // create new planes based on N such as a numnber
        if(NP > 0){
            tp  = new plane[N+NP];
        } else {
            tp = new plane[N];
            NP = N;
        }
        p = tp;
    };

    void printPlanes(){         // prints the planes information
        for(int i=0; i<NP; i++){
            cout << "Plane " << (i+1) << " --> ";
            for(int j=0; j< p[i].Npoints; j++ ){
                cout << "Point: P" << j <<" x: " << p[i].p[j].x << "  y: " << p[i].p[j].y << "  z: " << p[i].p[j].z << "   ";
            }
            cout << endl;
        }
    };

    void printTriangles(){      // prints triangles information

        for(int i=0; i<NP; i++){
            cout << "Plane " << (i+1) << ": ";
            p[i].printTriangleInformation();
        }
    };

    void printTrianglesTxt(){
        std::ofstream ofs("Simulation_Data/Room_Data/trianglesInformation.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());

        for(int i=0; i<NP; i++){
            p[i].printTriangleInformationTxt();
        }
        std::cout.rdbuf(cout_buff);
    };

    void printSubPlanesTxt(){
        std::ofstream ofs("Simulation_Data/Room_Data/subPlanesInformation.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());

        for(int i=0; i<NP; i++){
            p[i].printSubPlaneInformationTxt();
        }
        std::cout.rdbuf(cout_buff);
    };

    void printPlanesGridTxt(){
        std::ofstream ofs("Simulation_Data/Room_Data/planesGridInformation.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());
        for(int i=0; i<NP; i++){
            p[i].printPlaneGridInformationTxt();
        }
        std::cout.rdbuf(cout_buff);
    };

    vector NormalVector(vector v){  // gets the normal of a vector, it needs into raytracing algorithm
        double module;
        vector normal;
        module = v.Module();        // gets the module of a vector, method from vector
        if(module == 0){
            normal = 0;
        } else {
            normal = v / module;
        }
        return normal;
    };

    void MaxDistance(){             // max distance into the room for each "d" use to calculate the incident point
        maxDistance = 0.0f;         // helps to find a lost ray
        double distance = 0.0f;
        point p1;
        point p2;

        /*
            Distance = point[1].distance[2]
            In order to find the max distance, it needs to iter each possible point into each plane
            with all the points into all the planes. Then get the distance and update the max value

            iter planes:                        (int i)
                iter planes.points:             (int j)
                        planes[i].points[j]                 ---->  point to evaluate        P1
                    iter planes                 (int y)
                        iter planes.points      (int z)
                            planes[y].points[z]             ---->  all available points     P2
                            .
                            .
                            .
                            .
                            update max distance comparing

            if a incident point goes over the cube, so the distance should be greater than the max distance
         */

        for(int i=0; i<NP; i++){                                                // Go over all the planes

            for(int j=0; j < p[i].Npoints; j++){                                // Go over all the points
                    p1 = p[i].p[j];                                             // save the point
                for(int y=0; y<NP; y++){                                        // Go over all the planes again, with each point

                    for(int z=0; z<p[y].Npoints; z++){                          // Go over all the points too
                        p2 = p[y].p[z];                                         // save the point which will be compare to
                        distance = p1.distance(p2);                             // assing to  distance

                        if(maxDistance < distance){                             // evaluate if max changes and update it
                            maxDistance = distance;
                        }
                    }
                }
            }
        }
        maxDistanceRoom = maxDistance;
    }

    bool Inner(point p,triangle t) {            // return true if a point go over a triangle (Provided by the professor)
        double a1,a2,x,y,z,x0,y0,z0,x1,y1,z1,x2,y2,z2;

        x=p.x;
        y=p.y;
        z=p.z;

        x0=t.p0.x;
        y0=t.p0.y;
        z0=t.p0.z;
        x1=t.p1.x;
        y1=t.p1.y;
        z1=t.p1.z;
        x2=t.p2.x;
        y2=t.p2.y;
        z2=t.p2.z;

        if(t.Projection==yz) {                              //Proyeccion yz
            a1=-t.a0*(-y0*z+y2*z+y*z0-y2*z0-y*z2+y0*z2);
            a2=-t.a0*(y0*z-y1*z-y*z0+y1*z0+y*z1-y0*z1);
        }
        if(t.Projection==xz) {                              //Proyeccion xz
            a1=-t.a0*(-x0*z+x2*z+x*z0-x2*z0-x*z2+x0*z2);
            a2=-t.a0*(x0*z-x1*z-x*z0+x1*z0+x*z1-x0*z1);
        }
        if(t.Projection==xy) {                              //Proyeccion xy
            a1=-t.a0*(-x2*y0+x*y0+x0*y2-x*y2-x0*y+x2*y);
            a2=t.a0*(-x1*y0+x*y0+x0*y1-x*y1-x0*y+x1*y);
        }

        if((a1+a2<=1)&&(a1>=0)&&(a2>=0))
            return true;
        else
            return false;
    };

    void printNormalPlanes(){           // print normal planes information
        for(int i=0; i<NP; i++){
            cout << "Plane " << (i+1) << endl;
            p[i].printNormalInformation();
        }
    };

    void getSolidAngleTriangles(){
        int i, j, y, z;                                         // iter variables
        int m = numRoomTriangles;                               // total of triangles in the room
        double *solidAngleArea;                                 // vector total area into solid angle
        solidAngleArea = new double[numRoomTriangles-1];

        triangle t1;                                            // triangle
        triangle t2;                                            // triangle
        solidAngle.initiVector(m-1, m-1);

        for(i=0; i<NP; i++){                                    // iter over all the planes
            for(j=0; j<p[i].Ntriangles; j++){                   // iter over all the triangles
                t1 = p[i].t[j];                                 // save the first triangle
                solidAngleArea[t1.ID-1] = 0;                    // init the area
                for(y=0; y<NP; y++){                            // iter over all the planes
                    for(z=0; z<p[y].Ntriangles; z++){           // iter over all the triangles
                        t2 = p[y].t[z];                         // save the second triangle
                        if(i == y){                             // if the triangles are into the same plane, then the solid angle equals to zero
                            solidAngle.vec[t1.ID-1][t2.ID-1] = 0;

                        } else {
                            solidAngle.vec[t1.ID-1][t2.ID-1] = t2.solidAngle(t1.bc);        // find the slid angle between the first triangle and the center of the second one
                            solidAngleArea[t1.ID-1] += solidAngle.vec[t1.ID-1][t2.ID-1];    // add the found value to the sum of the area
                        }

                    }
                }
            }
        }


        for(i=0; i< p[0].Ntriangles*NP; i++){
            for(j=0; j< p[0].Ntriangles*NP; j++){
                solidAngle.vec[i][j] = solidAngle.vec[i][j] / solidAngleArea[i];       // iter over solid angle vector and normalize the data using solid angle area
            }
        }

    }

    void getSolidAngleReceptors(){
        int i, j;                                               // iter variables
        double *solidAngleArea;                                 // vector total area into solid angle
        solidAngleArea = new double[numRoomTriangles-1];
        int m = numRoomTriangles;

        triangle t1;                                            // triangle
        triangle t2;                                            // triangle
        solidAnglesReceptor.initiVector(numReceptors, m-1);


        for(i=0; i<NP; i++){                                    // iter over all the planes
            for(j=0; j<p[i].Ntriangles; j++){                   // iter over all the triangles
                t1 = p[i].t[j];
                solidAnglesReceptor.vec[numReceptors-1][t1.ID-1] = r[0].solidAngle(t1.bc);        // find the slid angle between the first triangle and the center of the second one
                solidAngleArea[numReceptors-1] += solidAnglesReceptor.vec[numReceptors-1][t1.ID-1];    // add the found value to the sum of the area
            }
        }

        double suma = 0;
        for(i=0; i< p[0].Ntriangles*NP; i++){
            for(j=0; j< numReceptors; j++){
                solidAnglesReceptor.vec[j][i] = solidAnglesReceptor.vec[j][i] / solidAngleArea[j];       // iter over solid angle vector and normalize the data using solid angle area
                suma += solidAngle.vec[j][i];

            }
        }

    }

    void getTrianglesDistances(){
        int i, j, y, z;
        int m = numRoomTriangles;

        distanceTriangles.initiVector(m-1, m-1);
        //vector triangleToTriangle;
        triangle t1;
        triangle t2;

        for(i=0; i<NP; i++){
            for(j=0; j<p[i].Ntriangles; j++){
                t1 = p[i].t[j];
                for(y=0; y<NP; y++){
                    for(z=0; z<p[y].Ntriangles; z++){
                        t2 = p[y].t[z];
                        if(i == y){
                            distanceTriangles.vec[t1.ID-1][t2.ID-1] = 0;
                        } else {
                            distanceTriangles.vec[t1.ID-1][t2.ID-1] = t1.bc.distance(t2.bc);
                        }

                    }
                }
            }
        }
    }

    void getReceptorsDistances(){
        int i, j;
        int m = numRoomTriangles;

        distanceReceptors.initiVector(m-1, numReceptors);
        //vector triangleToTriangle;
        triangle t1;
        triangle t2;

        for(i=0; i<NP; i++){                                    // iter over all the planes
            for(j=0; j<p[i].Ntriangles; j++){                   // iter over all the triangles
                t1 = p[i].t[j];
                for(int R=0; R<numReceptors; R++){
                    distanceReceptors.vec[R][t1.ID-1] = t1.bc.distance(r[R].position);
                }
            }
        }
    }

    void saveDistancesReceptors(){
        std::ofstream ofs("Simulation_Data/Room_Data/distancesReceptorsTriangles.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());
        cout << "[" << endl;
        for(int i=0; i< numReceptors; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << solidAnglesReceptor.vec[i][j] << ",  ";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
        std::cout.rdbuf(cout_buff);
    }

    void saveTimesFlyingR(){
        std::ofstream ofs("Simulation_Data/Room_Data/timeFlyingReceptorsTriangles.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());
        cout << "[" << endl;
        for(int i=0; i< numReceptors; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << timeFlyingR[i][j] << ",  ";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
        std::cout.rdbuf(cout_buff);
    }


    void addReceptor(int n){
        r = new receptor[n];
    }

    void findFlyingTimeReceptors(){
        int m = numRoomTriangles;

        // Init vector
        timeFlyingR = new int*[m];
        for(int n=0; n<numReceptors; n++){
            timeFlyingR[n] = new int[m];
        }

        for(int i=0; i< numRoomTriangles-1; i++){
            for(int j=0; j< numReceptors; j++){
                timeFlyingR[j][i] = 1000 * (distanceReceptors.vec[j][i] / 340);
            }
        }
        cout << endl;
    }

    void findFlyingTime(){
        int m = numRoomTriangles;

        // Init vector
        timeFlying = new int*[m];
        for(int n=0; n<m; n++){
            timeFlying[n] = new int[m];
        }

        for(int i=0; i< numRoomTriangles-1; i++){
            for(int j=0; j< numRoomTriangles-1; j++){
                if(distanceTriangles.vec[i][j] != 0){
                    timeFlying[i][j] = 1000 * (distanceTriangles.vec[i][j] / 340);
                } else {
                    timeFlying[i][j] = 0;
                }

            }
        }
    }


    void printSolidAngles(){
        cout << "[" << endl;
        for(int i=0; i< p[0].Ntriangles*NP; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << solidAngle.vec[i][j] << ",  ";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
    }

    void saveSolidAngles(){
        std::ofstream ofs("Simulation_Data/Room_Data/solidAnglesTriangles.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());
        cout << "[" << endl;
        for(int i=0; i< p[0].Ntriangles*NP; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << solidAngle.vec[i][j] << ",";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
        std::cout.rdbuf(cout_buff);
    }

    void printSolidAnglesR(){
        cout << "[" << endl;
        for(int i=0; i< numReceptors; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << solidAnglesReceptor.vec[i][j] << ",  ";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
    }

    void saveSolidAnglesR(){
        std::ofstream ofs("Simulation_Data/Room_Data/solidAnglesReceptorsTriangles.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());
        cout << "[" << endl;
        for(int i=0; i< numReceptors; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << solidAnglesReceptor.vec[i][j] << ",";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
        std::cout.rdbuf(cout_buff);
    }

    void printTriangleDistances(){
        cout << "[" << endl;
        for(int i=0; i< p[0].Ntriangles*NP; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << distanceTriangles.vec[i][j] << ",  ";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
    }

    void saveDataTriangleDistances(){
        std::ofstream ofs("Simulation_Data/Room_Data/distancesTriangles.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());
        cout << "[" << endl;
        for(int i=0; i< p[0].Ntriangles*NP; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << distanceTriangles.vec[i][j] << ",";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
        std::cout.rdbuf(cout_buff);
    }

    void printTimesFlying(){
        cout << "[" << endl;
        for(int i=0; i< p[0].Ntriangles*NP; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << timeFlying[i][j] << ",  ";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;

    }

    void saveTimesFlying(){
        std::ofstream ofs("Simulation_Data/Room_Data/timeFlyingTriangles.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());
        cout << "[" << endl;
        for(int i=0; i< p[0].Ntriangles*NP; i++){
            cout << "[";
            for(int j=0; j< p[0].Ntriangles*NP; j++){
                cout << timeFlying[i][j] << ",";
            }
            cout << "]," << endl;
        }
        cout << "]" << endl;
        std::cout.rdbuf(cout_buff);
    }

    void saveDataTrianglesBC(){
        std::ofstream ofs("Simulation_Data/Room_Data/TrianglesBC.txt");
        auto cout_buff = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());
        cout << "[";
        for(int i=0; i<NP; i++){
            for(int j=0; j< p[i].Ntriangles; j++){
                p[i].t[j].bc.printPoint();
            }
        }
        cout << "]" << endl;
        std::cout.rdbuf(cout_buff);
    }

    void printTrianglesBC(){
        cout << "[";
        for(int i=0; i<NP; i++){
            for(int j=0; j< p[i].Ntriangles; j++){
                p[i].t[j].bc.printPoint();
            }
        }
        cout << "]" << endl;
    }

    reflection *RayTracing(point sourcePosition, vector *Rays, int NRays, int numReflects){     // ray tracing algorithm
        /*
        In order to find the Point incident, it needs the formula Pi = O (origin) + Vi.d
        Pi = indicent point
        d = distance
        Vi = incident vector / Each Ray
        O = origin (source position)

        d it's found by using d = n.Vd / n.Vi
        Vd = vector resulting from origin point and one of the point of the plano
        */
        double distancia1 = 0;                          // d or distance
        double distancia2 = 0;                          // temporal distance use in the formula
        point punto1;                                   // Pi or incident point
        point punto2;
        point punto3;                                   // point used into the reflection
        vector vi;                                      // vi or incident vector
        point o;                                        // o - origin
        vector vd;                                      // vd or vector from origin to one of the plane points

        bool flagMaxReflection;                         // stop reflections when reachs the max
        int planoIntersecado;                           // intersected plane
        int trianguloIntersecado = -1;                   // intersected triangle
        int indiceReflecion = 1;                        // num of reflections
        int numReflections = numReflects;               // max reflections assigned
        reflection *rayTracing = new reflection[NRays]; // reflection struct vector
        double numerator, denominator;                  // variables to calculate d

        if(NP > 0){
            for(int R=0; R<NRays; R++){                             // go over all the rays

                vi = Rays[R];                                       // incident vector
                o = sourcePosition;                                 // origin -> source position
                trianguloIntersecado = -1;                          // assign value to -1, first value
                planoIntersecado = -1;

                rayTracing[R].N = 0;                                // Empezar con 0 reflexiones
                rayTracing[R].incidentPoint[0] = o;                 // Primer punto el origen de la fuente
                rayTracing[R].d[0] = 0.0;                           // first distance equals to 0
                rayTracing[R].Ray = R;                              // Num Ray
                rayTracing[R].actualEnergy[0] = sourceEnergy/12;    // Actual ray energy before reflection
                rayTracing[R].idTriangle[0] = -1;
                flagMaxReflection = false;

                while(!flagMaxReflection){                          // calculate all reflections, stop in the max assigned
                    MaxDistance();                                  // calculate max distance
                    distancia1 = maxDistance;                       // d = equals to max distance
                    for(int P=0; P<NP; P++){                        // go over the planes

                        if( (vi * p[P].n) < 0){                     // if scalar product < 0, then it could be a collision plane

                            vd = (p[P].p[0] - o);
                            numerator = p[P].n * vd;                // equals to n*vd
                            //cout << "numerator: " << numerator << endl;
                            denominator = p[P].n * vi;              // equals to n*vi
                            //cout << "denominator: " << denominator << endl;

                            if(denominator == 0){
                                distancia2 = -1;                    // if n*vi equals to 0, then d = -1
                            } else {
                                distancia2 = numerator / denominator;   // (n*vd / n*vi)
                            }


                            if((distancia2 > 0.000001) && (distancia2 < distancia1)){   // d can not be negative, and distance2 < distance

                                punto2 = o + (vi*distancia2);        // pi, incident point, calculate according the formula

                                for(int T = 0; T < p[P].Ntriangles; T++){   // go over the triangles

                                    if(Inner(punto2, p[P].t[T])){    // finds if the temporal point go over the triangle
                                        distancia1 = distancia2;     // assing the distance to d
                                        punto1 = punto2;             // assing the incident point to pi
                                        planoIntersecado = P;        // save the intersected plane
                                        trianguloIntersecado = p[P].t[T].ID;    // save the intersected triangle

                                        break;                       // stop the bucle
                                    }
                                }

                            }
                        }
                    }

                    if(distancia1 < maxDistance && planoIntersecado != -1){     // d needs to be < that max distance and plane different from -1
                        /*
                        Use this formula to calculate the result vector
                        Vr = Vi * ( n * ( Vi * n ) * 2 )
                        It's needs to be the normal from the collision plane

                        and this one to calculate the energy
                        E = Ei*(1-alpha)(1-delta)
                        */

                        // Save first information about reflection
                        punto3 = o;                                             // save the original source point
                        rayTracing[R].actualEnergy[indiceReflecion] = rayTracing[R].actualEnergy[indiceReflecion - 1] * (1 - alpha) * (1 - deltha);  // calculate the energy of the reflection
                        rayTracing[R].incidentPoint[indiceReflecion] = punto1;  // save the incident point into the structure
                        vector distance = punto1 - punto3;
                        rayTracing[R].d[indiceReflecion] = distance.Module();   // save the module  into the structure
                        rayTracing[R].N++;                                      // save the number of reflection  into the structure
                        rayTracing[R].idTriangle[indiceReflecion] = trianguloIntersecado;   // save the triangle id  into the structure
                        rayTracing[R].idPlano[indiceReflecion] = planoIntersecado;          // save the plane id  into the structure
                        rayTracing[R].vec[indiceReflecion] = distance;
                        rayTracing[R].mod[indiceReflecion] = distance.Module();
                        rayTracing[R].tim[indiceReflecion] = 1000 * (rayTracing[R].mod[indiceReflecion]/340);

                        // Calculate the next variables for reflection iteration
                        o = punto1;                                             // assing the new source point
                        vi = NormalVector( vi - ( p[planoIntersecado].n * ( vi * p[planoIntersecado].n * 2)));  // calculate the new incident vector or vr
                        indiceReflecion++;                                      // add the reflection number

                        // Verify max reflections
                        if(indiceReflecion > (numReflections-1)){       // stop bucle when reachs the max reflection number
                            flagMaxReflection = true;
                            indiceReflecion = 1;
                        }


                    } else {                // if d > max distance, so the ray is lost
                        indiceReflecion++;
                        rayTracing[R].lost = true;      // save the there are lost rays into the structure
                        punto3 = o + ( vi * maxDistance);
                        rayTracing[R].incidentPoint[indiceReflecion] = punto3;
                        rayTracing[R].N++;
                        flagMaxReflection = true;
                        indiceReflecion = 0;
                    }
                }
            }
        }
        return rayTracing;
    };
};

