//
// Created for experiment configuration
//

#include "Constants.h"

// Define the global experiment mode variable
// Default to Experiment 3 (current configuration: 6 objectives, MIN reliability, MAX affordability)
namespace Constants {
    int EXPERIMENT_MODE = 3;
    bool INCLUDE_SEVERITY = false; // Toggle severity objective via -V flag
    int TARGET_WSS_ID = 0;        // Target WSS for experiment 6 (default: WSS 0)
}
