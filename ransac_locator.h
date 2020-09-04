#ifndef RANSAC_LOCATOR_H
#define RANSAC_LOCATOR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>
#include <assert.h>

// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>

#define RANSAC_LEAST_CNT 4
#define CALCULABLE_LEAST_CNT 3
#define ENABLE_RANSAC_REFINE 1
#define PI 3.1415926535
#define EPS 0.0000001

typedef struct {
   double x;
   double y;
} Position2D;

double angle_normalize(double angle) {
    return angle >= 0 ? angle : angle + PI;
}

bool fequal(double f1, double f2) {
    return fabs(f1 - f2) < EPS;
}

bool is_landmarks_in_one_line(const Position2D landmarks[], const int effective_landmark_indexes[], int effective_landmark_num) {
    assert(effective_landmark_num > 2);
    double angle_0 = angle_normalize(atan2(landmarks[effective_landmark_indexes[0]].y - landmarks[effective_landmark_indexes[effective_landmark_num-1]].y, 
            landmarks[effective_landmark_indexes[0]].x - landmarks[effective_landmark_indexes[effective_landmark_num-1]].x));
    for (int i = 1; i < effective_landmark_num - 1; i++) {
        double angle_i = angle_normalize(atan2(landmarks[effective_landmark_indexes[i]].y - landmarks[effective_landmark_indexes[effective_landmark_num-1]].y, 
                landmarks[effective_landmark_indexes[i]].x - landmarks[effective_landmark_indexes[effective_landmark_num-1]].x));
        if (!fequal(angle_0, angle_i)) {
            return false;
        }
    }
    return true;
}

bool calculate_pos_with_index_(const Position2D landmarks[], const double distances[], const int effective_landmark_indexes[], int effective_landmark_num, 
        double* x_out, double* y_out) {
    if (is_landmarks_in_one_line(landmarks, effective_landmark_indexes, effective_landmark_num)) {
        printf("all landmarks in one line, cannot locate the position.\n");
        return false;
    }
    int matrix_row_num = effective_landmark_num - 1;
    int coefficient_num = 2;
    // Ax=b
    gsl_matrix* A = gsl_matrix_alloc(matrix_row_num, coefficient_num);
    gsl_matrix* a = gsl_matrix_alloc(matrix_row_num, coefficient_num);
    gsl_vector* b = gsl_vector_alloc(matrix_row_num);
    
    // result
	gsl_vector *sx = gsl_vector_alloc(coefficient_num);

    //Clear data
    gsl_matrix_set_zero(A);
    gsl_vector_set_zero(b);

    int index_last = effective_landmark_indexes[effective_landmark_num - 1];
    for (int i = 0; i < effective_landmark_num - 1; i++) {
        int index_i = effective_landmark_indexes[i];
        gsl_matrix_set(A, i, 0, (landmarks[index_i].x - landmarks[index_last].x) * 2);
        gsl_matrix_set(A, i, 1, (landmarks[index_i].y - landmarks[index_last].y) * 2);
        double bi = landmarks[index_i].x * landmarks[index_i].x - landmarks[index_last].x * landmarks[index_last].x
                        + landmarks[index_i].y * landmarks[index_i].y - landmarks[index_last].y * landmarks[index_last].y
                        + distances[index_last] * distances[index_last] - distances[index_i] * distances[index_i];
        gsl_vector_set(b, i, bi);
    }
    gsl_matrix_memcpy(a, A);
    gsl_vector* tau = gsl_vector_alloc(coefficient_num); //matrix_row_num, coefficient_num
    gsl_vector* residuals = gsl_vector_alloc(matrix_row_num);

    gsl_multifit_linear_workspace *w = gsl_multifit_linear_alloc(matrix_row_num, coefficient_num);
    gsl_multifit_linear_svd(A, w);
    double rcond = gsl_multifit_linear_rcond(w); // reciprocal condition number
    printf("conditional number is:%f\n", 1.0 / rcond);
    if (rcond < 1e-4) {
        printf("conditional number is too large, indicating that the problem is ill-conditioned.\n");
        return false;
    }
    double lambda_gcv = 0.1;
    double chisq, rnorm, snorm;
    gsl_multifit_linear_solve(0.0, A, b, sx, &rnorm, &snorm, w);

    *x_out = gsl_vector_get(sx, 0);
    *y_out = gsl_vector_get(sx, 1);

    gsl_matrix_free(A);
    gsl_matrix_free(a);
    gsl_vector_free(b);
    gsl_vector_free(tau);
    gsl_vector_free(sx);
    gsl_vector_free(residuals);
    return true;
}

bool calculate_pos_for_all_(const Position2D landmarks[], const double distances[], int effective_landmark_num, 
        double* x_out, double* y_out) {
    int* effective_landmark_indexes = (int*)malloc(sizeof(int) * effective_landmark_num);
    for (int i = 0; i < effective_landmark_num; i++) {
        effective_landmark_indexes[i] = i;
    }
    bool res = calculate_pos_with_index_(landmarks, distances, effective_landmark_indexes, effective_landmark_num, x_out, y_out);
    free(effective_landmark_indexes);
    return res;
}

bool getSampleIndexes(int index_size, int sample_indexes[], int each_sample_cnt) {
    // index_size > each_sample_cnt
    int * shuffled_indices = (int*) malloc(sizeof(int) * index_size);
    for (int i = 0; i < index_size; i++) {
        shuffled_indices[i] = i;
    }
    for (int i = 0; i < each_sample_cnt; i++) {
        // swap i and random index afterwards
        int rand_index = i + (rand() % (index_size - i));
        int tmp = shuffled_indices[i];
        shuffled_indices[i] = shuffled_indices[rand_index];
        shuffled_indices[rand_index] = tmp;

        sample_indexes[i] = shuffled_indices[i];
    }
    free(shuffled_indices);
    return true;
}

int get_all_fit_indexes(const Position2D landmarks[], const double distances[], int effective_landmark_num,
        const int sample_indexes[], int sample_size, double x, double y, double threshold, 
        int all_fit_indexes_out[], double* total_residual_out) {
    int* is_fit_vec = (int*)malloc(sizeof(int) * effective_landmark_num);
    memset(is_fit_vec, 0, sizeof(int) * effective_landmark_num);
    double* residuals = (double*)malloc(sizeof(double) * effective_landmark_num);
    memset(residuals, 0, sizeof(double) * effective_landmark_num);
    for (int i = 0; i < effective_landmark_num; i++) {
        if (is_fit_vec[i] == 0) {
            double estimated_distance = sqrt(pow((landmarks[i].x - x), 2) + pow((landmarks[i].y - y), 2));
            double difference_of_distances = fabs(distances[i] - estimated_distance);
            if (difference_of_distances < threshold) {
                is_fit_vec[i] = 1;
            }
            residuals[i] = difference_of_distances * difference_of_distances;
        }
    }
    int j = 0;
    int cnt = 0;
    *total_residual_out = 0;
    for (int i = 0; i < effective_landmark_num; i++) {
        if (is_fit_vec[i] == 1) {
            all_fit_indexes_out[j++] = i;
            cnt++;
            *total_residual_out += residuals[i];
        }
    }
    free(is_fit_vec);
    free(residuals);
    return cnt;
}

bool calculate_pos_robust_ransac(const Position2D landmarks[], const double distances[], int effective_landmark_num, 
        double* x_out, double* y_out) {
    bool success_res = false;
    int max_iterations = 10;
    double threshold = 0.1;
    if (effective_landmark_num < CALCULABLE_LEAST_CNT) {
        printf("effective_landmark_num too small, cant calculate.\n");
        success_res = false;
    } else if (effective_landmark_num <= RANSAC_LEAST_CNT) {
        printf("effective_landmark_num too small, calculate_pos_for_all_.\n");
        success_res = calculate_pos_for_all_(landmarks, distances, effective_landmark_num, x_out, y_out);
    } else { // ransac choose the best, then refine
        double residual_best = DBL_MAX;
        int fit_cnt_best = -1;

        int sample_indexes[RANSAC_LEAST_CNT];
        int* all_fit_indexes = (int*)malloc(sizeof(int) * effective_landmark_num);
        int* all_fit_indexes_best = (int*)malloc(sizeof(int) * effective_landmark_num);
        for (int i = 0; i < max_iterations; i++) {
            // choose sample landmarks
            printf("========iteration:%d===========\nchosen indexes:", i);
            getSampleIndexes(effective_landmark_num, sample_indexes, RANSAC_LEAST_CNT);
            for (int j = 0; j < RANSAC_LEAST_CNT; j++) {
                printf("%d,", sample_indexes[j]);
            }
            printf("\n");

            // printf("sample_indexes:\n");
            // for (int j = 0; j < RANSAC_LEAST_CNT; j++) {
            //     printf("%d\n", sample_indexes[j]);
            // }
            // calculate pos based on the chosen sample
            double x_out_tmp, y_out_tmp;
            bool success_once = calculate_pos_with_index_(landmarks, distances, sample_indexes, RANSAC_LEAST_CNT, &x_out_tmp, &y_out_tmp);
            if (success_once) {
                printf("x_out_tmp:%f, y_out_tmp:%f\n", x_out_tmp, y_out_tmp);

                double total_residual = DBL_MAX;
                double threshold = 0.1;
                int n = get_all_fit_indexes(landmarks, distances, effective_landmark_num, sample_indexes, RANSAC_LEAST_CNT, x_out_tmp, y_out_tmp, threshold, all_fit_indexes, &total_residual);

                printf("fit cnt:%d, residual:%f, all_fit_indexes:", n, total_residual);
                for (int i = 0; i < n; i++) {
                    printf("%d,", all_fit_indexes[i]);
                }
                printf("\n");
                if (n > fit_cnt_best || (n == fit_cnt_best && total_residual < residual_best)) {
                    residual_best = total_residual;
                    fit_cnt_best = n;
                    for (int j = 0; j < n; j++) {
                        all_fit_indexes_best[j] = all_fit_indexes[j];
                    }
                    *x_out = x_out_tmp;
                    *y_out = y_out_tmp;
                    success_res = true;
                }
            }
        }
        // refine model
        if (success_res) {
            printf("==============\nfit_cnt_best:%d, fit_cnt_best smaller than RANSAC_LEAST_CNT:%s\n", fit_cnt_best, (fit_cnt_best < RANSAC_LEAST_CNT) ? "true" : "false");
            printf("before refine[best]. result:%f, %f, residual:%f\n", *x_out, *y_out, residual_best);
            if (ENABLE_RANSAC_REFINE && fit_cnt_best >= RANSAC_LEAST_CNT) {
                double x_out_tmp, y_out_tmp;
                bool success_once = calculate_pos_with_index_(landmarks, distances, all_fit_indexes_best, fit_cnt_best, &x_out_tmp, &y_out_tmp);
                if (success_once) {

                    double residual_refined = 0.0;
                    for (int i = 0; i < fit_cnt_best; i++) {
                        double estimated_distance = sqrt(pow((landmarks[all_fit_indexes_best[i]].x - *x_out), 2) + pow((landmarks[all_fit_indexes_best[i]].y - *y_out), 2));
                        double difference_of_distances = fabs(distances[all_fit_indexes_best[i]] - estimated_distance);
                        residual_refined += difference_of_distances * difference_of_distances;
                    }
                    if (residual_refined < residual_best) { // order sensitive, so residual_refined may be big than residual_best
                        printf("after refine. result:%f, %f, residual:%f\n", *x_out, *y_out, residual_refined);
                        *x_out = x_out_tmp;
                        *y_out = y_out_tmp;
                    }
                }
            }
        }
        free(all_fit_indexes);
        free(all_fit_indexes_best);
    }

    return success_res;
}

#endif // RANSAC_LOCATOR_H