#include "Eff.h"
#include <fstream>
#include <iostream>
#include <algorithm>    // std::lower_bound, std::upper_bound, std::sort                                                                                                                                  
#include <vector>       // std::vector

double e(double lambda)
  {
    ifstream MirrorData;
    vector<double *> Holder;
    vector<double *> MirrorLambdas;
    vector<double *> MirrorEffs;
    string str;
    double value;


    inFile.open("MirrorRefl.txt");
    if (!inFile)
      {
        cerr << "Unable to open file MirrorRefl.txt";
        exit(1);   // call system to stop                                                                                                                                                                   
      }
  while (std::getline(MirrorData, str))
    {
      value = strtok(str,", n/");
      while (value != NULL)
        {

	  Holder.push_back (value);
          int i = 0;
          int j = 0;
          while (i < 2*Holder.Size+1)
            {
              MirrorEffs.push_back(Holder[2*i+1]);
              i++;
            }
          while (j < 2*Holder.Size)
            {
              int VecIndex;
              MirrorLambdas.push_back(Holder[2*j]);
              std::vector<int>::iterator up;
              up = std::upper_bound (0,Holder.Size,lambda);
              j++;
            }
          eff = *up;
         }
    inFile.close();
    return eff;
  }
