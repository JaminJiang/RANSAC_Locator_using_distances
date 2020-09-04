#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "ransac_locator.h"

#define LANDMARK_NUM 5

int main() {
    Position2D landmarks[LANDMARK_NUM];
    // landmarks[0].x = 0;
    // landmarks[0].y = 0;
    // landmarks[1].x = 0;
    // landmarks[1].y = 1;
    // landmarks[2].x = 1;
    // landmarks[2].y = 1;
    // landmarks[3].x = 1;
    // landmarks[3].y = 0;
    // landmarks[4].x = 2;
    // landmarks[4].y = 2;
    landmarks[0].x = 0.1;
    landmarks[0].y = 0;
    landmarks[1].x = 1;
    landmarks[1].y = 1;
    landmarks[2].x = 2;
    landmarks[2].y = 2;
    landmarks[3].x = 3;
    landmarks[3].y = 3;
    landmarks[4].x = 4;
    landmarks[4].y = 4;
    double distances[LANDMARK_NUM];
    srand(0);
    double groundtruth_x = rand() / (double)RAND_MAX;
    double groundtruth_y = rand() / (double)RAND_MAX;

    for (int i = 0; i < LANDMARK_NUM; i++) {
        distances[i] = sqrt(pow(groundtruth_x - landmarks[i].x, 2) + pow(groundtruth_y - landmarks[i].y, 2)) + (double)(rand()) / (double)RAND_MAX / 100;
        // distances[i] = sqrt(pow(groundtruth_x - landmarks[i].x, 2) + pow(groundtruth_y - landmarks[i].y, 2));
    }
    // distances[4] = distances[4] + 0.2;
    double estimated_x, estimated_y;
    bool success = calculate_pos_robust_ransac(landmarks, distances, LANDMARK_NUM, &estimated_x, &estimated_y);
    // bool calculate_pos_for_all_(landmarks, distances, effective_landmark_num, estimated_x, estimated_x);
    
    printf("\n===============\n");
    printf("groundtruth_x:%lf, groundtruth_y:%lf\n", groundtruth_x, groundtruth_y);
    printf("final result: success:%d, estimated_x:%lf, estimated_y:%lf\n", success, estimated_x, estimated_y);

    return 0;
}