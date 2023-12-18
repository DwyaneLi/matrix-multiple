#include <matrix_multiple.hh>

using namespace std;

template <typename T>
void init_matrix(T* matrix, int N) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(1.0, 100.0);
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            matrix[i * N + j] = dis(gen);
        }
    }
}

template <typename T>
void init_matrix_zero(T* matrix, int N) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            matrix[i * N + j] = 0;
        }
    }
}

void matrix_multiple_double(int argc, char** argv, int ctl_fd, int N) {
    double* m1 = (double*)aligned_alloc(64, sizeof(double) * N * N);
    double* m2 = (double*)aligned_alloc(64, sizeof(double) * N * N);
    double* m_res = (double*)aligned_alloc(64, sizeof(double) * N * N);

    init_matrix<double>(m1, N);
    init_matrix<double>(m2, N);
    init_matrix_zero<double>(m_res, N);

    int block_size = 0;

    switch(argv[2][0])
    {
    case 's':
        simple_matrix_multiple<double>(m1, m2, m_res, N, ctl_fd);
        break;

    case 'r':
        reorder_matrix_multiple<double>(m1, m2, m_res, N, ctl_fd);
        break;

    case 'b':
        if(argc != 4) {
            std::cerr << "Wrong number of arguments" << std::endl;
            std::exit(1); 
        }
        block_size = stoi(string(argv[3]));
        block_matrix_multiple<double>(m1, m2, m_res, N, block_size, ctl_fd);
        break;

    case 'v':
        if(argc == 4 && argv[3][0] == 'o') {
            vector_matrix_multiple_double_omp(m1, m2, m_res, N, ctl_fd);
        } else {
            vector_matrix_multiple_double(m1, m2, m_res, N, ctl_fd);
        }
        break;

    case 'w':
        if(argc < 4) {
            std::cerr << "Wrong number of arguments" << std::endl;
            std::exit(1); 
        }

        block_size = stoi(string(argv[3]));
        if(argc == 5 && argv[4][0] == 'o') {
            block_vector_matrix_multiple_double_omp(m1, m2, m_res, N, block_size, ctl_fd);
        } else {
            block_vector_matrix_multiple_double(m1, m2, m_res, N, block_size, ctl_fd);
        }
        break;

    default:
        std::cerr << "Wrong arguments" << std::endl;
        std::exit(1); 
        break;
    }
}

void matrix_multiple_float(int argc, char** argv, int ctl_fd, int N) {
    float* m1 = (float*)aligned_alloc(64, sizeof(float) * N * N);
    float* m2 = (float*)aligned_alloc(64, sizeof(float) * N * N);
    float* m_res = (float*)aligned_alloc(64, sizeof(float) * N * N);

    init_matrix<float>(m1, N);
    init_matrix<float>(m2, N);
    init_matrix_zero<float>(m_res, N);

    int block_size = 0;

    switch(argv[2][0])
    {
    case 's':
        simple_matrix_multiple<float>(m1, m2, m_res, N, ctl_fd);
        break;

    case 'r':
        reorder_matrix_multiple<float>(m1, m2, m_res, N, ctl_fd);
        break;

    case 'b':
        if(argc != 4) {
            std::cerr << "Wrong number of arguments" << std::endl;
            std::exit(1); 
        }
        block_size = stoi(string(argv[3]));
        block_matrix_multiple<float>(m1, m2, m_res, N, block_size, ctl_fd);
        break;

    case 'v':
        vector_matrix_multiple_float(m1, m2, m_res, N, ctl_fd);
        break;

    case 'w':
        if(argc != 4) {
            std::cerr << "Wrong number of arguments" << std::endl;
            std::exit(1); 
        }
        block_size = stoi(string(argv[3]));
        block_vector_matrix_multiple_float(m1, m2, m_res, N, block_size, ctl_fd);
        break;

    default:
        std::cerr << "Wrong arguments" << std::endl;
        std::exit(1); 
        break;
    }
}


int main(int argc, char* argv[]) {
    std::cout << "start" << std::endl;
    if(argc < 3) {
        std::cerr << "Wrong arguments" << std::endl;
        std::exit(1); 
    }
    const char* ctl_fd_env = "CTL_FD";
    auto str = getenv(ctl_fd_env);
    auto ctl_fd = stoi(str);
    string N_str(argv[1]);
    auto N = stoi(N_str);

    if(N == 4096) {
        matrix_multiple_double(argc, argv, ctl_fd, N);
    } else {
        matrix_multiple_float(argc, argv, ctl_fd, N);
    }
    std::cout << "end" << std::endl;
    return 0;
}