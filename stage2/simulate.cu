#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include "constants.hpp"
#include <array>
#include <cuda_runtime.h>

static constexpr int MAX_PER_REGION = (REGIONS_SINGLE * REGIONS_SINGLE / (MIN_MASS * MIN_MASS));

void check_error(cudaError_t err, const char* msg)
{
    if (err != cudaSuccess)
    {
        fprintf(stderr, "%s :%s!\n", msg, cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

struct CudaBalls
{
    float *x;
    float *y;
    float *v_x;
    float *v_y;
    float *mass;
    int *index;

    int size;
};

__device__
bool is_collision(float x1, float y1, float r1, float x2, float y2, float r2)
{
    float sumRadius = r1 + r2;

    float xd = x1 - x2;
    float yd = y1 - y2;

    float sqrRadius = sumRadius * sumRadius;
    float distSqr = (xd * xd) + (yd * yd);

    return distSqr <= sqrRadius; 
}

struct wall_distance
{
    float l_wall, r_wall, t_wall, b_wall;
};

__device__
wall_distance wall_collision_point(float x, float y)
{
    return {x - LEFT_WALL,
            RIGHT_WALL - x,
            TOP_WALL - y,
            y - BOTTOM_WALL};
}

__device__
void collide_wall(float &x, float &y, float m, float &v_x, float &v_y)
{
    auto wall_dists = wall_collision_point(x, y);
    float collision_x = x;
    float collision_y = y;

    if (wall_dists.l_wall < m) {
        x = LEFT_WALL + m;
        collision_x = LEFT_WALL;
    } else if (wall_dists.r_wall < m) {
        x = RIGHT_WALL - m;
        collision_x = RIGHT_WALL;
    /* only consider collision with one wall */
    } else if (wall_dists.t_wall < m) {
        y = TOP_WALL - m;
        collision_y = TOP_WALL;
    } else if (wall_dists.b_wall < m) {
        y = BOTTOM_WALL + m;
        collision_y = BOTTOM_WALL;
    } else {
        return;
    }

    auto dot = v_x * (x - collision_x) + v_y * (y - collision_y);
    auto norm = (x - collision_x) * (x - collision_x) + (y - collision_y) * (y - collision_y);
    auto dot_norm = dot / norm;

    v_x -= 2 * dot_norm * (x - collision_x);
    v_y -= 2 * dot_norm * (y - collision_y);
}

__device__
void collide(float &x1, float &y1, float m1, float &v_x1, float &v_y1, float &x2, float& y2, float m2, float &v_x2, float &v_y2)
{
    if (!is_collision(x1, y1, m1, x2, y2, m2))
        return;

    float overlap_x = x2 - x1;
    float overlap_y = y2 - y1;

    float overlap_dist = sqrt(overlap_x * overlap_x + overlap_y * overlap_y);
    float radius_sum = m1 + m2;
    float coef = radius_sum / overlap_dist;

    /* move x2 away from x1 so that x2 - x1 == r1 + r2 */
    x2 = x1 + coef * overlap_x;
    y2 = y1 + coef * overlap_y;

    auto norm = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    norm = norm != 0 ? norm : 1;

    auto dot1 = (v_x1 - v_x2) * (x1 - x2) + (v_y1 - v_y2) * (y1 - y2);
    auto dot2 = (v_x2 - v_x1) * (x2 - x1) + (v_y2 - v_y1) * (y2 - y1);

    v_x1 -= (2 * m2) / (m1 + m2) * dot1 / norm * (x1 - x2);
    v_y1 -= (2 * m2) / (m1 + m2) * dot1 / norm * (y1 - y2);

    v_x2 -= (2 * m1) / (m1 + m2) * dot2 / norm * (x2 - x1);
    v_y2 -= (2 * m1) / (m1 + m2) * dot2 / norm * (y2 - y1);
}

__global__ void
clear_regions(int *regions_size)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i >= REGIONS_NUM)
        return;

        regions_size[i] = 0;
}

__global__ void
split_to_regions(CudaBalls balls, int *regions, int *regions_size)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i >= balls.size)
        return;

    int region_x = (balls.x[i] + WIDTH / 2) / REGIONS_SINGLE;
    int region_y = (balls.y[i] + HEIGHT / 2) / REGIONS_SINGLE;

    int r = region_x + region_y * REGIONS_X;
    int size = atomicAdd(&regions_size[r], 1);

    regions[r * MAX_PER_REGION + size] = i;
}

__global__ void
do_collide(CudaBalls balls, int *regions, int *regions_size)
{
    int r = blockIdx.x;
    int i = threadIdx.x;

    if (i >= regions_size[r])
        return;

    auto ball_index = regions[r * MAX_PER_REGION + i];

    int neighbour_regions[] = {1, REGIONS_X, REGIONS_X + 1};
    // neighbour regions
    for (int rs = 0; rs < 3; rs++) {
        int neighbour = r + neighbour_regions[rs];

        if (neighbour < 0 || neighbour >= REGIONS_NUM)
            continue;

        for (int j = 0; j < regions_size[neighbour]; j++) {
            auto nball_index = regions[neighbour * MAX_PER_REGION + j];
            collide(balls.x[ball_index],
                    balls.y[ball_index],
                    balls.mass[ball_index],
                    balls.v_x[ball_index],
                    balls.v_y[ball_index],
                    balls.x[nball_index],
                    balls.y[nball_index],
                    balls.mass[nball_index],
                    balls.v_x[nball_index],
                    balls.v_y[nball_index]);
        }
    }

    // same region
    for (int j = i + 1; j < regions_size[r]; j++) {
            auto nball_index = regions[r * MAX_PER_REGION + j];
            collide(balls.x[ball_index],
                    balls.y[ball_index],
                    balls.mass[ball_index],
                    balls.v_x[ball_index],
                    balls.v_y[ball_index],
                    balls.x[nball_index],
                    balls.y[nball_index],
                    balls.mass[nball_index],
                    balls.v_x[nball_index],
                    balls.v_y[nball_index]);
    }
}

__global__ void
collide_wall(CudaBalls _balls)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= _balls.size)
        return;

    collide_wall(_balls.x[i],
        _balls.y[i],
        _balls.mass[i],
        _balls.v_x[i],
        _balls.v_y[i]);
}

__global__ void
advance(CudaBalls balls)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= balls.size)
        return;

    balls.x[i] += balls.v_x[i];
    balls.y[i] += balls.v_y[i];
}

template <typename T>
static T* aligned_alloc(int n) {
    void* tmp = 0;

    auto alignment = sizeof(T) > sizeof(void*) ? sizeof(T) : sizeof(void*);
    if (posix_memalign(&tmp, alignment, sizeof(T) * n)) {
        throw std::bad_alloc();
    }
    return (T*)tmp;
}

struct Balls
{
    Balls(int num = MAX_PER_REGION): size(0), capacity(num) {
        x = aligned_alloc<float>(num);
        y = aligned_alloc<float>(num);
        v_x = aligned_alloc<float>(num);
        v_y = aligned_alloc<float>(num);
        mass = aligned_alloc<float>(num);
        index = aligned_alloc<int>(num);
    }

    ~Balls() {
        std::free(x);
        std::free(y);
        std::free(v_x);
        std::free(v_y);
        std::free(mass);
        std::free(index);
    }

    float *x;
    float *y;
    float *v_x;
    float *v_y;
    float *mass;
    int *index;

    int size;
    int capacity;

    char padding[8];
};

class simulation
{
public:
    simulation(int n_balls, int seed): _balls(n_balls)
    {
        srand(seed);

        x_regions = (WIDTH - 2 * MARGIN - 2 * MAX_MASS) / (2 * MAX_MASS);
        y_regions = (HEIGHT - 2 * MARGIN - 2 * MAX_MASS) / (2 * MAX_MASS);

        int n_regions = x_regions * y_regions;
        std::vector<bool> pos_bitmap(n_regions, false);

        for (int i = 0; i < n_balls; i++) {
            auto position = get_random_pos(pos_bitmap);

            _balls.x[i] = position.first;
            _balls.y[i] = position.second;
            _balls.v_x[i] = get_random(MIN_SPEED, MAX_SPEED);
            _balls.v_y[i] = get_random(MIN_SPEED, MAX_SPEED);
            _balls.mass[i] = get_random(MIN_MASS, MAX_MASS);
            _balls.index[i] = i;
        }

        _balls.size = n_balls;

        threadsPerBlock = 256;
        blocksPerGrid = (n_balls + threadsPerBlock - 1) / threadsPerBlock;

        regions = alloc_cuda_regions();
        regions_size = alloc_cuda_regions_size();
    }

    ~simulation()
    {
        cudaError_t err = cudaSuccess;

        err = cudaFree(regions);
        check_error(err, "free regions");

        err = cudaFree(regions_size);
        check_error(err, "free regions_size");
    }

    void step(int iters = 10)
    {
        CudaBalls c_balls = alloc_cuda_balls(_balls.size);
        memcpy_balls_to_device(c_balls, _balls);

        // #pragma omp for collapse(2)
        // for (int r1 = 0; r1 < REGIONS_X; r1 +=2) {
        //     for (int r2 = 0; r2 < REGIONS_X; r2 +=2) {
        //         do_collide(r1 * REGIONS_X + r2);
        //     }
        // }

        // #pragma omp for collapse(2)
        // for (int r1 = 1; r1 < REGIONS_X; r1 +=2) {
        //     for (int r2 = 0; r2 < REGIONS_X; r2 +=2) {
        //         do_collide(r1 * REGIONS_X + r2);
        //     }
        // }

        // #pragma omp for collapse(2)
        // for (int r1 = 0; r1 < REGIONS_X; r1 +=2) {
        //     for (int r2 = 1; r2 < REGIONS_X; r2 +=2) {
        //         do_collide(r1 * REGIONS_X + r2);
        //     }
        // }

        // #pragma omp for collapse(2)
        // for (int r1 = 1; r1 < REGIONS_X; r1 += 2) {
        //     for (int r2 = 1; r2 < REGIONS_X; r2 +=2) {
        //         do_collide(r1 * REGIONS_X + r2);
        //     }
        // }

        // do_collide<<REGIONS_X, 

        for (int i = 0; i < iters; i++) {
            clear_regions<<<blocksPerGrid, threadsPerBlock>>>(regions_size);
            split_to_regions<<<blocksPerGrid, threadsPerBlock>>>(c_balls, regions, regions_size);

            do_collide<<<REGIONS_NUM, 128>>>(c_balls, regions, regions_size);
            collide_wall<<<blocksPerGrid, threadsPerBlock>>>(c_balls);

            advance<<<blocksPerGrid, threadsPerBlock>>>(c_balls);
            cudaError_t err = cudaGetLastError();
            check_error(err, "count kernel");
        }

        memcpy_balls_to_host(c_balls, _balls);
        free_cuda_balls(c_balls);
    }

    const Balls& balls()
    {
        return _balls;
    }

private:
    int *alloc_cuda_regions()
    {
        int *regions;

        cudaError_t err = cudaSuccess;
        err = cudaMalloc((void **)&regions, REGIONS_NUM * MAX_PER_REGION * sizeof(int));
        check_error(err, "table allocation regions");

        return regions;
    }

    int *alloc_cuda_regions_size()
    {
        int *regions_s;

        cudaError_t err = cudaSuccess;
        err = cudaMalloc((void **)&regions_s, REGIONS_NUM * sizeof(int));
        check_error(err, "table allocation regions_s");

        return regions_s;
    }

    CudaBalls alloc_cuda_balls(int size)
    {
        CudaBalls c_balls;
        c_balls.size = size;

        cudaError_t err = cudaSuccess;
        err = cudaMalloc((void **)&c_balls.x, size * sizeof(float));
        check_error(err, "table allocation x");
        err = cudaMalloc((void **)&c_balls.y, size * sizeof(float));
        check_error(err, "table allocation y");
        err = cudaMalloc((void **)&c_balls.v_x, size * sizeof(float));
        check_error(err, "table allocation vx");
        err = cudaMalloc((void **)&c_balls.v_y, size * sizeof(float));
        check_error(err, "table allocation vy");
        err = cudaMalloc((void **)&c_balls.mass, size * sizeof(float));
        check_error(err, "table allocation mass");
        err = cudaMalloc((void **)&c_balls.index, size * sizeof(int));
        check_error(err, "table allocation index");

        return c_balls;
    }

    void memcpy_balls_to_device(CudaBalls &c_balls, Balls &_balls)
    {
        cudaError_t err = cudaSuccess;
        err = cudaMemcpy(c_balls.x, _balls.x, _balls.size * sizeof(float), cudaMemcpyHostToDevice);
        check_error(err, "table memcpy to device x");
        err = cudaMemcpy(c_balls.y, _balls.y, _balls.size * sizeof(float), cudaMemcpyHostToDevice);
        check_error(err, "table memcpy to device y");
        err = cudaMemcpy(c_balls.v_x, _balls.v_x, _balls.size * sizeof(float), cudaMemcpyHostToDevice);
        check_error(err, "table memcpy to device vx");
        err = cudaMemcpy(c_balls.v_y, _balls.v_y, _balls.size * sizeof(float), cudaMemcpyHostToDevice);
        check_error(err, "table memcpy to device vy");
        err = cudaMemcpy(c_balls.mass, _balls.mass, _balls.size * sizeof(float), cudaMemcpyHostToDevice);
        check_error(err, "table memcpy to device mass");
        err = cudaMemcpy(c_balls.index, _balls.index, _balls.size * sizeof(int), cudaMemcpyHostToDevice);
        check_error(err, "table memcpy to device index");
    }

    void memcpy_balls_to_host(CudaBalls &c_balls, Balls &_balls)
    {
        cudaError_t err = cudaSuccess;
        err = cudaMemcpy(_balls.x, c_balls.x, _balls.size * sizeof(float), cudaMemcpyDeviceToHost);
        check_error(err, "table memcpy to host x");
        err = cudaMemcpy(_balls.y, c_balls.y, _balls.size * sizeof(float), cudaMemcpyDeviceToHost);
        check_error(err, "table memcpy to host y");
        err = cudaMemcpy(_balls.v_x, c_balls.v_x, _balls.size * sizeof(float), cudaMemcpyDeviceToHost);
        check_error(err, "table memcpy to host vx");
        err = cudaMemcpy(_balls.v_y, c_balls.v_y, _balls.size * sizeof(float), cudaMemcpyDeviceToHost);
        check_error(err, "table memcpy to host vy");
        err = cudaMemcpy(_balls.mass, c_balls.mass, _balls.size * sizeof(float), cudaMemcpyDeviceToHost);
        check_error(err, "table memcpy to host mass");
        err = cudaMemcpy(_balls.index, c_balls.index, _balls.size * sizeof(int), cudaMemcpyDeviceToHost);
        check_error(err, "table memcpy to host index");
    }

    void free_cuda_balls(CudaBalls &c_balls)
    {
        cudaError_t err = cudaSuccess;

        err = cudaFree(c_balls.x);
        check_error(err, "free x");
        err = cudaFree(c_balls.y);
        check_error(err, "free y");
        err = cudaFree(c_balls.v_x);
        check_error(err, "free vx");
        err = cudaFree(c_balls.v_y);
        check_error(err, "free vy");
        err = cudaFree(c_balls.mass);
        check_error(err, "free mass");
        err = cudaFree(c_balls.index);
        check_error(err, "free index");
    }

    float get_random(float lo, float hi)
    {
        return lo + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(hi-lo)));
    }

    std::pair<float, float> get_random_pos(std::vector<bool> &pos_bitmap)
    {
            int pos_random;
            do {
                pos_random = rand() % pos_bitmap.size();
            } while (pos_bitmap[pos_random] == true);

            pos_bitmap[pos_random] = true;

            float x = (float) ((pos_random % x_regions) * 2 * MAX_MASS + 2 * MAX_MASS - WIDTH / 2); 
            float y = (float) ((pos_random / x_regions) * 2 * MAX_MASS + 2 * MAX_MASS - HEIGHT / 2);

            return {x, y};
    }

    int x_regions, y_regions;
    Balls _balls;

    int *regions;
    int *regions_size;

    int threadsPerBlock;
    int blocksPerGrid;
};

// int
// main(int argc, char *argv[])
// {
//     // Error code to check return values for CUDA calls
//     cudaError_t err = cudaSuccess;

//     size_t freem, totalm;
//     err = cudaMemGetInfo(&freem, &totalm);
//     check_error(err, "get info");
//     printf("Opt1: Free memory %lu, total memory %lu\n", freem, totalm);

//     // Reset the device and exit
//     err = cudaDeviceReset();
//     check_error(err, "device reset");

//     printf("Done\n");
//     return 0;
// }


