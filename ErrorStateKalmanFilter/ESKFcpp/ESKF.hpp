#ifndef ESKF_HPP // Include guard
#define ESKF_HPP
// #pragma once // Include guard
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include "math_utility.hpp"
#include <vector>



using namespace Eigen;
using namespace std;

// ---------------------------------- Algorithm Parameters ----------------------------------
double latitude = 40.923061376 * DEG_TO_RAD;                          // Current Latitude  [Rad]
double longitude = 29.31453938 * DEG_TO_RAD;                          // Current Longitude [Rad]
double height = 9.9355;                                               // Current Height    [M]
double average_time = 10.0;                                           // Algorithm Initial Avarage Time [Sn]
const int algorithm_frequency = 200;                                  // Algorithm Frequency [Hz]        (Default Value of IFOS~500 -> 200 Hz)
const double sampling_periode = 1.0 / algorithm_frequency;            // Algorithm Sampling periode [Sn] (Default Value of IFOS~500 -> 5 Ms)
int number_of_state = 17;                                             // Algorithm Sate Vector Size      (Default Value of ESKF -> 17 "Attitude(3x1) Pose(3x1) Velocity(3x1) WanderSin(1x1) WanderCos(1x1) Ba(3x1) Bg(3x1)")
int model_order = 3;                                                  // State Transition Taylor Series Order  (Default Value of ESKF -> 3)
double yaw_angle_offset = 360;                                        // Yaw Angle Offset Value [Deg]          (Default Value of ESKF -> 360)
double wander_ang = 0;                                                // Wander Angle Initial Value [Rad]
double wander_check = 1;                                              // Wander Correction Check
double th_min = 1e-1;                                                 // Data Clean Threashold Min Value
double th_max = 2e1;                                                  // Data Clean Threashold Max Value
double roll_ini = 0.0;                                                // Initial Roll Value for Static Alignment
double pitch_ini = 0.0;                                               // Initial Pitch Value for Static Alignment
double yaw_ini = 0.0;                                                 // Initial Yaw Value for Static Alignment



#endif // #ifndef ESKF_HPP