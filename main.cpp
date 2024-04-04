#include <iostream>
#include "cabecera.h"
using namespace std;

void createRoom(){

};

int main()
{


    //---------------------------------------------------------------------------------------------------------------------
    // Instances
    room r;                     // Create room (cube in this case)
    source s;                   // Create source

    //double  distance, time;     // Enviorenment variables such as: distance, time, some one add here, etc

    //---------------------------------------------------------------------------------------------------------------------
    // Setting up the source

    point sourcePosition;
    //sourcePosition.x = 1.5f;
    //sourcePosition.y = 1.5f;
    //sourcePosition.z = 1.5f;
    int numRays;                        // Num of Rays, setted to 12
    double alfa;                              // Alpha coeficient
    double delta;                            // Delta coeficient

    s.sourcePosition = sourcePosition;

    cout << "Proyecto Modelos y Simulacion" << endl;
    cout << "Integrantes:" << endl;
    cout << "- Brian Coyago" << endl;
    cout << "- Fabricio Simbania" << endl;
    cout << "- Daniel Lazo" << endl;
    cout << "- Samir Zurita" << endl;
    cout << "- Bruce Soto" << endl;
    cout << endl;

    cout << "------------------Inicio de la simulacion------------------" << endl;

    cout << "Ingrese la posicion de la fuente " << endl;
    cout << "x: ";
    cin >> s.sourcePosition.x;
    cout << "y: ";
    cin >> s.sourcePosition.y;
    cout << "z: ";
    cin >> s.sourcePosition.z;

    cout << endl;
    cout << "Ingrese las variables del entorno" << endl;
    cout << "coeficiente alfa: ";
    cin >> alfa;
    cout << "coeficiente delta: ";
    cin >> delta;
    cout << "Energia de la fuente: ";
    cin >> s.eS;
    cout << "Cantidad de rayos (1-12): ";
    cin >> numRays;

    cout << endl;
    cout << "--------------------Datos del entorno--------------------" << endl;
    cout << "--------------------Rayos creados con sus direcciones--------------------" << endl;
    cout << endl;
    s.NRAYS = numRays;
    s.createRays(numRays);   // Create rays -> 12
    //Displaying results
    cout << "Source created with " << s.eS << " of energy and " << alfa << ", " << delta << " coeficients" << endl;
    cout << "Number of rays: " << s.NRAYS << endl;
    for (int i = 0; i < s.NRAYS; i++) {
        cout << "Ray " << i << ": (" << s.Rays[i].x << ", " << s.Rays[i].y << ", " << s.Rays[i].z << ")" << endl;
    }
    cout << endl;

    //---------------------------------------------------------------------------------------------------------------------
    //Setting up the room. (In this case a cube of 6 planes)
    int numPlanes = 6;
    int numPoints = 4;
    int numSubPlanes;

    cout << "Ingrese el numero de subplanos para la simulacion " << endl;
    cout << "cantidad: ";
    cin >> numSubPlanes;

    r.NewPlanes(numPlanes);                     // Create 6 planes for the room instance
    r.sourceEnergy = s.eS;
    r.alpha = alfa;
    r.deltha = delta;
    //-------------square 1 back
    r.p[0].NewPoints(numPoints);                // Create points for the plane
    r.p[0].p[0].x=-2.0f;
    r.p[0].p[0].y=2.0f;
    r.p[0].p[0].z=2.0f;
    r.p[0].p[1].x=-2.0f;
    r.p[0].p[1].y=-2.0f;
    r.p[0].p[1].z=2.0f;
    r.p[0].p[2].x=-2.0f;
    r.p[0].p[2].y=-2.0f;
    r.p[0].p[2].z=-2.0f;
    r.p[0].p[3].x=-2.0f;
    r.p[0].p[3].y=2.0f;
    r.p[0].p[3].z=-2.0f;
    //r.p[0].PointGetTriangle();
    r.p[0].genMoreTriangle(numSubPlanes);

    //-------------square 2 front
    r.p[1].NewPoints(numPoints);                // Create points for the plane
    r.p[1].p[0].x=2.0f;
    r.p[1].p[0].y=2.0f;
    r.p[1].p[0].z=2.0f;
    r.p[1].p[1].x=2.0f;
    r.p[1].p[1].y=2.0f;
    r.p[1].p[1].z=-2.0f;
    r.p[1].p[2].x=2.0f;
    r.p[1].p[2].y=-2.0f;
    r.p[1].p[2].z=-2.0f;
    r.p[1].p[3].x=2.0f;
    r.p[1].p[3].y=-2.0f;
    r.p[1].p[3].z=2.0f;
    //r.p[1].PointGetTriangle();
    r.p[1].genMoreTriangle(numSubPlanes);

    //-------------square 3 left
    r.p[2].NewPoints(numPoints);                // Create points for the plane
    r.p[2].p[0].x=-2.0f;
    r.p[2].p[0].y=-2.0f;
    r.p[2].p[0].z=2.0f;
    r.p[2].p[1].x=2.0f;
    r.p[2].p[1].y=-2.0f;
    r.p[2].p[1].z=2.0f;
    r.p[2].p[2].x=2.0f;
    r.p[2].p[2].y=-2.0f;
    r.p[2].p[2].z=-2.0f;
    r.p[2].p[3].x=-2.0f;
    r.p[2].p[3].y=-2.0f;
    r.p[2].p[3].z=-2.0f;
    r.p[2].genMoreTriangle(numSubPlanes);
    //r.p[2].PointGetTriangle();

    //-------------square 4 right
    r.p[3].NewPoints(numPoints);                // Create points for the plane
    r.p[3].p[0].x=2.0f;
    r.p[3].p[0].y=2.0f;
    r.p[3].p[0].z=2.0f;
    r.p[3].p[1].x=-2.0f;
    r.p[3].p[1].y=2.0f;
    r.p[3].p[1].z=2.0f;
    r.p[3].p[2].x=-2.0f;
    r.p[3].p[2].y=2.0f;
    r.p[3].p[2].z=-2.0f;
    r.p[3].p[3].x=2.0f;
    r.p[3].p[3].y=2.0f;
    r.p[3].p[3].z=-2.0f;
    r.p[3].genMoreTriangle(numSubPlanes);
    //r.p[3].PointGetTriangle();

    //-------------square 5 top
    r.p[4].NewPoints(numPoints);                // Create points for the plane
    r.p[4].p[0].x=-2.0f;
    r.p[4].p[0].y=-2.0f;
    r.p[4].p[0].z=2.0f;
    r.p[4].p[1].x=-2.0f;
    r.p[4].p[1].y=2.0f;
    r.p[4].p[1].z=2.0f;
    r.p[4].p[2].x=2.0f;
    r.p[4].p[2].y=2.0f;
    r.p[4].p[2].z=2.0f;
    r.p[4].p[3].x=2.0f;
    r.p[4].p[3].y=-2.0f;
    r.p[4].p[3].z=2.0f;
    r.p[4].genMoreTriangle(numSubPlanes);
    //r.p[4].PointGetTriangle();

    //-------------square 1 bottom
    r.p[5].NewPoints(numPoints);                // Create points for the plane
    r.p[5].p[0].x=-2.0f;
    r.p[5].p[0].y=2.0f;
    r.p[5].p[0].z=-2.0f;
    r.p[5].p[1].x=-2.0f;
    r.p[5].p[1].y=-2.0f;
    r.p[5].p[1].z=-2.0f;
    r.p[5].p[2].x=2.0f;
    r.p[5].p[2].y=-2.0f;
    r.p[5].p[2].z=-2.0f;
    r.p[5].p[3].x=2.0f;
    r.p[5].p[3].y=2.0f;
    r.p[5].p[3].z=-2.0f;
    r.p[5].genMoreTriangle(numSubPlanes);
    //r.p[5].PointGetTriangle();

    cout << "--------------------Informacion de los planos--------------------" << endl;
    cout << endl;
    cout << "Number of planes: " << numPlanes << ", with " << numPoints << " points each one" << endl;
    r.printPlanes();                            // print the result of the planes
    r.printSubPlanesTxt();
    r.printPlanesGridTxt();

    // Calcular normales, código entregado por el docente
    //Inicializar normales de los planos
    // se calculan las normales con la normal de su primer triangulo
    // se calcula el baricentro de los triángulos
    int cont_t=1;
    for (int i=0; i<r.NP; i++) {
        r.p[i].n=TriangleNormal(r.p[i].t[0]);
        for (int j=0;j<r.p[i].Ntriangles;j++){
            r.p[i].t[j].CalculateProjection();
            r.p[i].t[j].Centroid();
            r.p[i].t[j].ID = cont_t;
            cont_t++;
        }
    }

    cout << endl;
    cout << endl;
    cout << "Number of planes: " << numPlanes << ", with 2 triangles each one" << endl;
    r.printTrianglesTxt();
    r.printTriangles();

    cout << "--------------------Informacion de las normales--------------------" << endl;
    cout << endl;
    cout << "Number of planes: " << numPlanes << ", with 1 normal each one" << endl;
    r.printNormalPlanes();


    //---------------------------------------------------------------------------------------------------------------------
    // Proceso Raytracing
    cout << endl;
    cout << "--------------------Proceso RayTracing iniciado--------------------" << endl;
    cout << "Datos obtenidos" << endl;
    cout << "Punto origen fuente:  x: " << s.sourcePosition.x << ",  y: " << s.sourcePosition.y << ",  z: " << s.sourcePosition.z << endl;
    int numReflections;
    cout << "Ingrese la cantidad de reflexiones: ";
    cin >> numReflections;
    numReflections++;
    cout << endl;

    reflection *ry;
    ry = r.RayTracing(s.sourcePosition, s.Rays, s.NRAYS, numReflections);

    cout << "--------------------Proceso RayTracing Finalizado--------------------" << endl;
    cout << endl;


    //---------------------------------------Implementaciones deber matrices------------------------------------------------------------------------------
    cout << "----------------Txt file with center of the triangles generated----------------"<<endl;
    r.saveDataTrianglesBC();
    //cout << endl;

    r.getSolidAngleTriangles();
    cout << "----------------Txt file tringles percen generated----------------"<<endl;
    r.saveSolidAngles();
    //cout << endl;

    r.getTrianglesDistances();
    cout << "----------------Txt file with distances between triangles generated----------------"<<endl;
    r.saveDataTriangleDistances();
    //cout << endl;

    r.findFlyingTime();
    cout << "----------------Txt file with flying time generated----------------"<<endl;
    r.saveTimesFlying();
    //cout << endl;

    //---------------------------------------Implementaciones transiciones de energia------------------------------------------------------------------------------

    //cout << "----------------Transicion de energia parte 1----------------"<<endl;


    int timeSimulation = 1000;
    DoubleVector MET;
    MET.initiVector(numRoomTriangles-1, timeSimulation);
    MET = 0;

    // Init vectors
    DoubleVector mod_distance;
    mod_distance.initiVector(numRays, numReflections);
    DoubleVector mod_tiempo;
    mod_tiempo.initiVector(numRays, numReflections);


    receptor r1;
    r1.position.x = 1.5;
    r1.position.y = 1.5;
    r1.position.z = -1.5;
    r1.initEnergy(timeSimulation);
    r1.diskRadius = 0.3;
    r.addReceptor(1);
    r.r[0] = r1;

    // Receptor data
    r.getSolidAngleReceptors();
    r.saveSolidAnglesR();

    r.getReceptorsDistances();
    r.saveDistancesReceptors();

    r.findFlyingTimeReceptors();
    r.saveTimesFlyingR();


    double Er;
    for(int R=0; R<numRays; R++){
        Er = ry[R].actualEnergy[0];
        mod_distance.vec[R] = ry[R].mod;
        mod_tiempo.vec[R] = ry[R].tim;
        mod_distance.vec[R][0] = 0;
        mod_tiempo.vec[R][0] = 0;

        for(int i=1; i<ry[R].N; i++){

            mod_distance.vec[R][i] = mod_distance.vec[R][i-1] + mod_distance.vec[R][i];
            mod_tiempo.vec[R][i] = mod_tiempo.vec[R][i-1] + mod_tiempo.vec[R][i];

            int tri, tim;
            tim = mod_tiempo.vec[R][i-1];
            tri = ry[R].idTriangle[i];

            if(tim > timeSimulation) tim = 999;

            MET.vec[tri-1][tim] += ry[R].actualEnergy[i-1] * (1-alfa)*delta;
            //MET.vec[tri-1][tim] += ry[R].actualEnergy[i-1] * pow((1-alfa)*delta, i-1);

            //calcular para receptor
            r.r[0].addEnergy(i, R, 340, mod_tiempo.vec[R][i-1], ry[R].actualEnergy[i-1], ry, mod_distance.vec[R][i-1]);

        }
    }

    // save MET data into a txt file
    MET.saveData();
    r.r[0].saveData();

    //cout << "----------------Transicion de energia parte 2----------------"<<endl;


    int totalRoomTriangles = numRoomTriangles - 1;
    double absCoefPlane = (1-alfa);

    for(int t=0; t<timeSimulation; t++){
        for(int k=0; k<totalRoomTriangles; k++){
                if(MET.vec[k][t]!=0){
                    int tim;
                    double porcentaje;
                    for(int tri=0; tri<totalRoomTriangles; tri++){
                        if(r.timeFlying[k][tri]!=0){
                            tim = r.timeFlying[k][tri] + t;
                            if(tim>timeSimulation) tim = timeSimulation;
                            porcentaje = r.solidAngle.vec[k][tri];
                            MET.vec[tri][tim] += MET.vec[k][t] * porcentaje *absCoefPlane;
                        }
                    }

                    for(int Recep=0; Recep<numReceptors; Recep++){
                        int tim = r.timeFlyingR[Recep][k] + t;
                        if(tim > timeSimulation) tim = timeSimulation;
                        r.r[0].resultEnergy[tim] += MET.vec[k][t] * r.solidAnglesReceptor.vec[Recep][k];
                    }

                }

        }
    }

    MET.saveData();
    r.r[0].saveData();

    int makeD;
    do{
        cout << endl;
        cout << "------------------------Transicion de energia finalizada------------------"<<endl;
        cout << "Se obtuvo los siguientes datos: " << endl;
        cout << "1. Datos RayTracing" << endl;
        cout << "2. Centroides de los triangulos en la sala" << endl;
        cout << "3. Distancia entre los triangulos en la sala y de los receptores" << endl;
        cout << "4. Porcentaje energia para cada triangulo" << endl;
        cout << "5. Matriz espacio temporal" << endl;
        cout << "6. Energia de los receptores" << endl;
        cout << "7. Salir" << endl;
        cout << "Que datos desea ver: ";
        cin >> makeD;
        cout << endl;

        switch(makeD){
        case 1:
            cout << "--------------------Datos RayTracing--------------------" << endl;
            cout << endl;
            for(int Ray=0; Ray < numRays; Ray++){
                cout << "Ray " << Ray << " , number of reflections: " << ry[Ray].N << endl;
                for(int Reflection=0; Reflection < ry[Ray].N; Reflection++){
                    cout << "Incident point " << (Reflection + 1) << ":  (" << ry[Ray].incidentPoint[Reflection].x << ", " << ry[Ray].incidentPoint[Reflection].y << ", " << ry[Ray].incidentPoint[Reflection].z << ")";
                    cout << ",  Energy: " << ry[Ray].actualEnergy[Reflection] << ", with Plane: " << ry[Ray].idPlano[Reflection] << ", with Triangle: " << ry[Ray].idTriangle[Reflection];
                    cout << ", distance: " << ry[Ray].d[Reflection] << endl;
                }
                cout << endl;
            }
            break;
        case 2:
            cout << "--------------------Baricentros triangulos--------------------" << endl;
            cout << endl;
            r.printTrianglesBC();
            break;
        case 3:
            cout << "--------------------Distancias triangulosXtriangulos--------------------" << endl;
            cout << endl;
            r.printTriangleDistances();
            break;
        case 4:
            cout << "--------------------Angulos solidos de los triangulos--------------------" << endl;
            cout << endl;
            r.printSolidAngles();
            cout << "--------------------Angulos solidos de los receptores--------------------" << endl;
            cout << endl;
            r.printSolidAnglesR();

            break;
        case 5:
            cout << "Dada la cantidad de datos, se creo un documento en la siguiente ruta";
            cout << "./Simulation_Data/Room_Data/matrizEspacioTemporal.txt" << endl;
            break;
        case 6:
            cout << "Dada la cantidad de datos, se creo un documento en la siguiente ruta";
            cout << "./Simulation_Data/Receptor_Data/data_x.y_z.txt" << endl;
            break;
        case 7:
            cout << "Bye bye" << endl;
            break;

        default:
            cout << "Ingrese una opcion valida" << endl;
            break;
        }

    }while(makeD != 7);

    return 0;
}
