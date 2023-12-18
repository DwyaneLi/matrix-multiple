#include <stdio.h>
#include <malloc.h>
#include <cstdlib>
#include <random>
#include <iostream>
#include <ctime>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <fstream>
#include <unistd.h>
#include <mmintrin.h>
#include <xmmintrin.h>
#include <immintrin.h>
#include <omp.h>

using namespace std;


void perf_on(int ctl_fd) {
    auto s = write(ctl_fd, "enable", 6);
}

void perf_off(int ctl_fd) {
    auto s = write(ctl_fd, "disable", 7);
}

// 普通乘法
template <typename T>
void simple_matrix_multiple(T* m1, T* m2, T* m_res, int n, int ctl_fd) {
    perf_on(ctl_fd);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
                m_res[i * n + j] += m1[i * n + k] * m2[k * n + j]; 
            }
        }
    }
    perf_off(ctl_fd);
}

// 重排序
template <typename T>
void reorder_matrix_multiple(T* m1, T* m2, T* m_res, int n, int ctl_fd) {
    perf_on(ctl_fd);
    for(int i = 0; i < n; i++) {
        for(int k = 0; k < n; k++) {
            for(int j = 0; j < n; j++) {
                m_res[i * n + j] += m1[i * n + k] * m2[k * n + j]; 
            }
        }
    }
    perf_off(ctl_fd);
}

// 分块乘法
template <typename T>
void block_matrix_multiple(T* m1, T* m2, T* m_res, int N, int K, int ctl_fd) {
    int iTile = K, jTile = K, kTile = K;
    int iOuterBound = N / iTile, jOuterBound = N / jTile, kOuterBound = N / kTile;
    perf_on(ctl_fd);
    for (int i_outer = 0; i_outer < iOuterBound; i_outer++) {
        for (int j_outer = 0; j_outer < jOuterBound; j_outer++) {
            for (int k_outer = 0; k_outer < kOuterBound; k_outer++) {
                for (int i_inner = 0; i_inner < iTile; i_inner++) {
                    for (int k_inner = 0; k_inner < kTile; k_inner++) {
                        for (int j_inner = 0; j_inner < jTile; j_inner++) {
                            m_res[(i_outer * iTile + i_inner) * N + (j_outer * jTile + j_inner)] +=
                            m1[(i_outer * iTile + i_inner) * N + (k_outer * kTile + k_inner)] *
                            m2[(k_outer * kTile + k_inner) * N + (j_outer * jTile + j_inner)];
                        }
                    }
                }
            }
        }
    }
    perf_off(ctl_fd);
}

// 向量流（double）
// 一个向量8个元素
void vector_matrix_multiple_double(double* m1, double* m2, double* m_res, int N, int ctl_fd) {
    __m512d c, a, b, d;
    int doubleNum = 512 / (sizeof(double) * 8);
    perf_on(ctl_fd);
    for (int i_inner = 0; i_inner < N; i_inner++) {
        for (int k_inner = 0; k_inner < N; k_inner++) { 
            a = _mm512_set1_pd(*(m1 + i_inner * N + k_inner));
            for (int _j_inner = 0; _j_inner < N / doubleNum; _j_inner++) {
                b = _mm512_load_pd(m2 + k_inner * N +  _j_inner * doubleNum);                
                c = _mm512_mul_pd(b, a);
                d = _mm512_load_pd(m_res + i_inner * N + _j_inner * doubleNum);
                c = _mm512_add_pd(c, d);
                _mm512_store_pd(m_res + i_inner * N + _j_inner * doubleNum, c);
            }
        }
    }
    perf_off(ctl_fd);
}

// 向量流（double）+ omp
// 一个向量8个元素
void vector_matrix_multiple_double_omp(double* m1, double* m2, double* m_res, int N, int ctl_fd) {
    __m512d c, a, b, d;
    int doubleNum = 512 / (sizeof(double) * 8);
    perf_on(ctl_fd);
#pragma omp parallel for num_threads(32)
    for (int i_inner = 0; i_inner < N; i_inner++) {
        for (int k_inner = 0; k_inner < N; k_inner++) { 
            a = _mm512_set1_pd(*(m1 + i_inner * N + k_inner));
            for (int _j_inner = 0; _j_inner < N / doubleNum; _j_inner++) {
                b = _mm512_load_pd(m2 + k_inner * N +  _j_inner * doubleNum);                
                c = _mm512_mul_pd(b, a);
                d = _mm512_load_pd(m_res + i_inner * N + _j_inner * doubleNum);
                c = _mm512_add_pd(c, d);
                _mm512_store_pd(m_res + i_inner * N + _j_inner * doubleNum, c);
            }
        }
    }
    perf_off(ctl_fd);
}

// 向量流（float）
// 一个向量8个元素
void vector_matrix_multiple_float(float* m1, float* m2, float* m_res, int N, int ctl_fd) {
    __m256 c, a, b, d;
    int floatNum = 256 / (sizeof(float) * 8);
    perf_on(ctl_fd);
    for (int i_inner = 0; i_inner < N; i_inner++) {
        for (int k_inner = 0; k_inner < N; k_inner++) { 
            a = _mm256_set1_ps(*(m1 + i_inner * N + k_inner));
            for (int _j_inner = 0; _j_inner < N / floatNum; _j_inner++) {
                b = _mm256_load_ps(m2 + k_inner * N +  _j_inner * floatNum);  
                //std::cout << "lala1" << std::endl;              
                c = _mm256_mul_ps(b, a);
                d = _mm256_load_ps(m_res + i_inner * N + _j_inner * floatNum);
                c = _mm256_add_ps(c, d);
                _mm256_store_ps(m_res + i_inner * N + _j_inner * floatNum, c);
            }
        }
    }
    perf_off(ctl_fd);
}

// 向量流（double）+ 分块
// 一个向量8个元素
void block_vector_matrix_multiple_double(double* m1, double* m2, double* m_res, int N, int K, int ctl_fd) {
    __m512d c, a, b, d;
    int iTile = K, jTile = K, kTile = K;
    int iOuterBound = N / iTile, jOuterBound = N / jTile, kOuterBound = N / kTile;
    perf_on(ctl_fd);
    int doubleNum = 512 / (sizeof(double) * 8);
    for (int i_outer = 0; i_outer < iOuterBound; i_outer++) {
        for (int j_outer = 0; j_outer < jOuterBound; j_outer++) {
            for (int k_outer = 0; k_outer < kOuterBound; k_outer++) {
                for (int i_inner = 0; i_inner < iTile; i_inner++) {
                    for (int k_inner = 0; k_inner < kTile; k_inner++) { 
                        a = _mm512_set1_pd(*(m1 + (i_outer * iTile + i_inner) * N +
                                            (k_outer * kTile + k_inner)));
                        for (int _j_inner = 0; _j_inner < jTile / doubleNum; _j_inner++) {
                            b = _mm512_load_pd(m2 + (k_outer * kTile + k_inner) * N +
                                (j_outer * jTile + _j_inner * doubleNum));                
                            c = _mm512_mul_pd(b, a);
                            d = _mm512_load_pd(m_res + (i_outer * iTile + i_inner) * N +
                                (j_outer * jTile + _j_inner * doubleNum));
                            c = _mm512_add_pd(c, d);
                            _mm512_store_pd(m_res + (i_outer * iTile + i_inner) * N +
                                    (j_outer * jTile + _j_inner * doubleNum), c);
                        }
                    }
                }
            }
        }
    }
    perf_off(ctl_fd);
}

// 向量流（float）+ 分块 + omp
// 一个向量8个元素
void block_vector_matrix_multiple_double_omp(double* m1, double* m2, double* m_res, int N, int K, int ctl_fd) {
    __m512d c, a, b, d;
    int iTile = K, jTile = K, kTile = K;
    int iOuterBound = N / iTile, jOuterBound = N / jTile, kOuterBound = N / kTile;
    perf_on(ctl_fd);
    int doubleNum = 512 / (sizeof(double) * 8);
#pragma omp parallel for num_threads(32)
    for (int i_outer = 0; i_outer < iOuterBound; i_outer++) {
        for (int j_outer = 0; j_outer < jOuterBound; j_outer++) {
            for (int k_outer = 0; k_outer < kOuterBound; k_outer++) {
                for (int i_inner = 0; i_inner < iTile; i_inner++) {
                    for (int k_inner = 0; k_inner < kTile; k_inner++) { 
                        a = _mm512_set1_pd(*(m1 + (i_outer * iTile + i_inner) * N +
                                            (k_outer * kTile + k_inner)));
                        for (int _j_inner = 0; _j_inner < jTile / doubleNum; _j_inner++) {
                            b = _mm512_load_pd(m2 + (k_outer * kTile + k_inner) * N +
                                (j_outer * jTile + _j_inner * doubleNum));                
                            c = _mm512_mul_pd(b, a);
                            d = _mm512_load_pd(m_res + (i_outer * iTile + i_inner) * N +
                                (j_outer * jTile + _j_inner * doubleNum));
                            c = _mm512_add_pd(c, d);
                            _mm512_store_pd(m_res + (i_outer * iTile + i_inner) * N +
                                    (j_outer * jTile + _j_inner * doubleNum), c);
                        }
                    }
                }
            }
        }
    }
    perf_off(ctl_fd);
}

// 向量流（float）+ 分块
// 一个向量8个元素
void block_vector_matrix_multiple_float(float* m1, float* m2, float* m_res, int N, int K, int ctl_fd) {
    __m256 c, a, b, d;
    int iTile = K, jTile = K, kTile = K;
    int iOuterBound = N / iTile, jOuterBound = N / jTile, kOuterBound = N / kTile;
    perf_on(ctl_fd);
    int floatNum = 256 / (sizeof(float) * 8);
    for (int i_outer = 0; i_outer < iOuterBound; i_outer++) {
        for (int j_outer = 0; j_outer < jOuterBound; j_outer++) {
            for (int k_outer = 0; k_outer < kOuterBound; k_outer++) {
                for (int i_inner = 0; i_inner < iTile; i_inner++) {
                    for (int k_inner = 0; k_inner < kTile; k_inner++) { 
                        a = _mm256_set1_ps(*(m1 + (i_outer * iTile + i_inner) * N +
                                            (k_outer * kTile + k_inner)));
                        for (int _j_inner = 0; _j_inner < jTile / floatNum; _j_inner++) {
                            b = _mm256_load_ps(m2 + (k_outer * kTile + k_inner) * N +
                                (j_outer * jTile + _j_inner * floatNum));                
                            c = _mm256_mul_ps(b, a);
                            d = _mm256_load_ps(m_res + (i_outer * iTile + i_inner) * N +
                                (j_outer * jTile + _j_inner * floatNum));
                            c = _mm256_add_ps(c, d);
                            _mm256_store_ps(m_res + (i_outer * iTile + i_inner) * N +
                                    (j_outer * jTile + _j_inner * floatNum), c);
                        }
                    }
                }
            }
        }
    }
    perf_off(ctl_fd);
}      
