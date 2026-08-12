#include <iostream>
#include <math.h>
#include <time.h>
#include <vector>
#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <set>
#include <utility>
#include <cstdio>
using namespace std;

const char* filename = "../../b_diff/data_gil_4.425.bin";

void delete_file(const char* filename) {
    if (std::remove(filename) != 0) {
        std::perror("Error deleting file");
    }
}
float find_max(const std::vector<double>& data){
  float max = 0;
  for (int i =0; i < data.size(); i++){
    if (max < data[i]){
      max = data[i];
    }
  }
  return max;
}


std::vector<int> findPeaks(const std::vector<double>& data, int windowSize, int distance) {
    std::vector<double> smoothedData;
    for (int i = 0; i < data.size(); i++) {
        int start = std::max(0, i - windowSize);
        int end = std::min(static_cast<int>(data.size()), i + windowSize + 1);
        double sum = 0.0;
        for (int j = start; j < end; j++) {
            sum += data[j];
        }
        double average = sum / (end - start);
        smoothedData.push_back(average);
    }

    std::vector<int> peaks;
    int prev_idx = 0;
    float threshold = find_max(smoothedData) * 0.3; 
    for (int i = 1; i < smoothedData.size() - 1; i++) {
        if (smoothedData[i] > smoothedData[i - 1] && smoothedData[i] > smoothedData[i + 1] && i - 0 > distance && smoothedData[i] > threshold) {
            peaks.push_back(i);
            prev_idx = i;
        }
    }

    return peaks;
}



std::pair<double, double> f1(double x1,double x2, double a, double b, double c, double d, double v){
  double resp = a*v + c*(1/(v*v))*x1*x1*x2;
  double resn =  (b+d)*x1;
  return std::make_pair(resp, resn);
}

std::pair<double, double> f2(double x1,double x2, double b, double c, double v){
  double resp = b * x1;
  double resn = c/(v*v)*x1*x1*x2;
  return std::make_pair(resp, resn);
}

std::vector<double> fn(double x1, double x2, double y1, double y2, double a, double b, double c, double d, double Vtot){
  double V1 = Vtot*(x1 + y1) / (x1 + x2 + y1 + y2);
  double V2 = Vtot*(x2 + y2) / (x1 + x2 + y1 + y2);
  // std::cout << V1 << "  " << V2 << std::endl;
  //V1 = Vtot;
  //V2 = Vtot;


  vector<double> pval = vector<double>(8,0);
  
  pval[0] = a*V1;
  pval[1] = b*x1;
  pval[2] = d*x1;
  pval[3] = c / (V1*V1) *x1*x1*y1;
  pval[4] = a*V2;
  pval[5] = b*x2;
  pval[6] = d*x2;
  pval[7] = c / (V2*V2) *x2*x2*y2;
  return pval;
}


int single_run(){
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> distrib(0, 1);

  // double a0 = 3;
  // double b0 = 8;
  // double c0 = 3;
  // double d0 = 8;
  // double a0 = 2;
  // double b0 = 8;
  // double c0 = 2;
  // double d0 = 8;
  double a0 = 2.;
  double b0 = 4.425;
  double c0 = 1.0;
  double d0 = 1.0;
  // double a0 = 1;
  // double b0 = 1.8;
  // double c0 = 1;
  // double d0 = 1;
  double Vtot = 5000;

  double RealTime = 0;
  // double vol = .1*pow(10,-14);
  // double NA = 6.03*pow(10,23);
  // double cal = NA*vol*pow(10,-6);

  double x01 = a0/d0*0.3*Vtot ;
  double x02 = b0*d0/(c0*a0)*0.3*Vtot ;

  double x11 = a0/d0*0.8*Vtot ;
  double x12 = b0*d0/(c0*a0)*0.8*Vtot ;

  // x01 = Vtot*.1;
  // x02 = Vtot*.1;
  // x11 = Vtot*.1;
  // x12 = Vtot*.1;


  // double v0 = .1;
  // double v1 = .9;


  int max_step = 100000000;
  vector<double> X01 = vector<double>(max_step,0);
  vector<double> X02 = vector<double>(max_step,0);
  vector<double> X11 = vector<double>(max_step,0);
  vector<double> X12 = vector<double>(max_step,0);
  vector<double> V0 = vector<double>(max_step,0);
  vector<double> V1 = vector<double>(max_step,0);
  vector<double> RT = vector<double>(max_step,0);



  int idx = 0;
  int idx_vec = 0;



  vector<double> dt = vector<double>(8,0);

  // std::cout << v0 << "  "<< v1 << std::endl;
  double tau = 0;

  for (int i =0; i < max_step*100; i++){
     vector<double> RES = fn(x01, x11, x02, x12, a0, b0, c0, d0, Vtot);

     double total_rate = RES[0] + RES[1] + RES[2] + RES[3] + RES[4] + RES[5] + RES[6] + RES[7];
     for (int j = 0; j < 8; j++){
       RES[j] = RES[j] / total_rate;
       // std::cout << RES[j] << std::endl;
     }
     //   dt[j] = 1/RES[j] * log(1/distrib(gen));
     //   std::cout << dt[j] << std::endl;
     // }
     // std::cout << "_______________________" << std::endl;
     tau = -log(distrib(gen)) / total_rate;
     // std::cout << tau << std::endl;
     int max_iter = 0;
     double sum = 0;
     double prob = distrib(gen);
     for (int j = 0; j <8; j++){
       sum = sum + RES[j];
       // std::cout << sum << "   " << prob << std::endl;
       if (sum > prob){
         max_iter = j;
         break;
       }
     }
     if (max_iter == 0){
       x01 += 1;
     }
     else if (max_iter == 1){
       x01 -= 1;
       x02 += 1;
     }
     else if (max_iter == 2){
       x01 -= 1;
     }
     else if (max_iter == 3){
       x01 += 1;
       x02 -= 1;
     }

     else if (max_iter == 4){
       x11 += 1;
     }
     else if (max_iter == 5){
       x11 -= 1;
       x12 += 1;
     }
     else if (max_iter == 6){
       x11 -= 1;
     }
     else if (max_iter == 7){
       x12 -= 1;
       x11 += 1;
     }
     // RealTime += dt[max_iter];
     RealTime += tau;

     if (idx%100 == 0){
       X01[idx_vec] = x01;
       X02[idx_vec] = x02;
       X11[idx_vec] = x11;
       X12[idx_vec] = x12;
       V0[idx_vec] = (x01 + x02) / (x01 + x11 + x02 + x12);
       V1[idx_vec] = (x11 + x12) / (x01 + x11 + x02 + x12);
       RT[idx_vec] = RealTime;
       idx_vec += 1;
     }
    idx += 1;
  
  }
  // delete_file("peaks.bin");
  // std::ofstream outPeak("peaks.bin", std::ios::binary | std::ios::app);
  // std::vector<int> peakInds = findPeaks(X1, 200, 100);
  // outPeak.write(reinterpret_cast<const char*>(peakInds.data()), peakInds.size() * sizeof(int));
  // outPeak.close();

  // delete_file("scatter_gil_500.bin");
  // delete_file("Canard_gil_5000.bin");
  // std::ofstream outFile("scatter_gil_500.bin", std::ios::binary | std::ios::app);
  std::ofstream outFile(filename, std::ios::binary | std::ios::app);
  outFile.write(reinterpret_cast<const char*>(X01.data()), X01.size() * sizeof(double));
  outFile.write(reinterpret_cast<const char*>(X02.data()), X02.size() * sizeof(double));
  outFile.write(reinterpret_cast<const char*>(X11.data()), X11.size() * sizeof(double));
  outFile.write(reinterpret_cast<const char*>(X12.data()), X12.size() * sizeof(double));
  outFile.write(reinterpret_cast<const char*>(V0.data()), V0.size() * sizeof(double));
  outFile.write(reinterpret_cast<const char*>(RT.data()), RT.size() * sizeof(double));
  outFile.close();

  return 1;

}




int main(){

  int res = single_run();

  return 0;

}
