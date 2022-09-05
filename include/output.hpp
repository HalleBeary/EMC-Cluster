#ifndef output_HPP
#define output_HPP

#include "treeStructure.hpp"



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


class Input_OutputClass
{

    public:
        void k_reader();
        void output_writer(OctreeClass *tree, int tsim);
    
    
};  

    
    void Input_OutputClass::k_reader()
    {

    ifstream ifile("./100000_ps.csv");

    string line;

    while (getline(ifile, line))
    {
        istringstream linestring{line}; // creates a string from the 'getted' line

        vector <string> k_values;

        string k_val;

        while( getline(linestring, k_val, ' '))
        {
            k_values.push_back(k_val);
        }
        
        k_values.erase(k_values.begin(), k_values.begin() + 9);
        k_values.erase(k_values.end() - 3, k_values.end());
        k_values.erase(k_values.begin() + 1, k_values.begin() + 3);
        k_values.erase(k_values.begin() + 2, k_values.begin() + 4);

        //TODO Convert to Doubles and save into physicalValues :)
 
        for (string elem : k_values)
           cout << "[" << elem << "]" ;
        cout << endl;
    }

    

}

    void Input_OutputClass::output_writer(OctreeClass *tree, int tsim)
    {

    FILE *output1;

    std::string s = std::to_string(tsim / 10);


    char file1[64];
    sprintf(file1, "../../Output/writer_test4/%d_ps.csv", tsim);
    
    output1 = fopen(file1,"w");


    long double* pos_x = new long double ;
    long double* pos_y = new long double ;
    long double* pos_z = new long double ;

    long double* k_x = new long double ;
    long double* k_y = new long double ;
    long double* k_z = new long double ;
    
    long double* potential = new long double;


    tree->forEachLeaf([&](typename OctreeClass::LeafClass* Leaf)
    { 
        for (FSize idxPart = 0 ; idxPart < Leaf->getTargets()->getNbParticles() ; idxPart ++){
        

        *pos_x = Leaf->getTargets()->getPositions()[0][idxPart];
        *pos_y = Leaf->getTargets()->getPositions()[1][idxPart];
        *pos_z = Leaf->getTargets()->getPositions()[2][idxPart];    

        *k_x = Leaf->getTargets()->getPhysicalValues(1,1)[idxPart]; 
        *k_y = Leaf->getTargets()->getPhysicalValues(2,1)[idxPart];
        *k_z = Leaf->getTargets()->getPhysicalValues(3,1)[idxPart];

        *potential = Leaf->getTargets()->getPotentials()[idxPart];
        
        fprintf(output1, "%Lf   %Lf   %Lf   %Lf   %Lf   %Lf   %Lf\n", *pos_x,  *pos_y,  *pos_z, *k_x, *k_y, *k_z, *potential) ;

        }
        
    });


    delete pos_x; 
    delete pos_y;
    delete pos_z;
    delete k_x;
    delete k_y;
    delete k_z;
    delete potential;
    
    }







/*

void writer_output(OctreeClass *tree, int tsim)
{   

    FILE *output1;

    int i = 0; // sizeof particles gives size of particle array per particle which is 8 (pos, forces, pot, charge)

    tree->forEachLeaf([&](typename OctreeClass::LeafClass* Leaf) 
    { 
        for (FSize idxPart = 0 ; idxPart < Leaf->getTargets()->getNbParticles() ; idxPart ++){

            ++i;

        }
    
    });

    std::string s = std::to_string(tsim / 10);


    char file1[64];
    sprintf(file1, "../../Output/writer_test4/%d_ps.csv", tsim);
    
    output1 = fopen(file1,"w");


    long double* pos_x = new long double [i];
    long double* pos_y = new long double [i];
    long double* pos_z = new long double [i];

    long double* mom_x = new long double [i];
    long double* mom_y = new long double [i];
    long double* mom_z = new long double [i];
    
    long double* potential = new long double[i];


    tree->forEachLeaf([&](typename OctreeClass::LeafClass* Leaf) 
    { 
        for (FSize idxPart = 0 ; idxPart < Leaf->getTargets()->getNbParticles() ; idxPart ++){
        

        long double x = Leaf->getTargets()->getPositions()[0][idxPart];
        long double y = Leaf->getTargets()->getPositions()[1][idxPart];
        long double z = Leaf->getTargets()->getPositions()[2][idxPart];    

        long double kx = Leaf->getTargets()->getPhysicalValues(1,1)[idxPart]; // 1 ev => k^2 = 1.3e19 
        long double ky = Leaf->getTargets()->getPhysicalValues(2,1)[idxPart];
        long double kz = Leaf->getTargets()->getPhysicalValues(3,1)[idxPart];

        long double pot = Leaf->getTargets()->getPotentials()[idxPart];

        pos_x[idxPart] = x;
        pos_y[idxPart] = y;
        pos_z[idxPart] = z;

        mom_x[idxPart] = kx;
        mom_y[idxPart] = ky;
        mom_z[idxPart] = kz;

        potential[idxPart] = pot;
        
        fprintf(output1, "%Lf   %Lf   %Lf   %Lf   %Lf   %Lf   %Lf\n", pos_x[idxPart],  pos_y[idxPart],  pos_z[idxPart], mom_x[idxPart], mom_y[idxPart], mom_z[idxPart], potential[idxPart]);



        // fclose(output1);

        }
    
    });
   // std::cout << "output loop counter: " << j << std::endl;
    
    delete[] pos_x; 
    delete[] pos_y;
    delete[] pos_z;
    delete[] mom_x;
    delete[] mom_y;
    delete[] mom_z;
    delete[] potential;
    
}
*/

        

#endif
