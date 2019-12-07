#pragma once

#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <omp.h>

#include "constants.hpp"
#include <array>

typedef float float8_t __attribute__ ((vector_size (8 * sizeof(float))));
typedef int vint8_t __attribute__ ((vector_size (8 * sizeof(int))));

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
    Balls(int num = REGIONS_SINGLE * REGIONS_SINGLE / (MIN_MASS * MIN_MASS)): size(0), capacity(num) {
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

wall_distance wall_collision_point(float x, float y)
{
    return {x - LEFT_WALL,
            RIGHT_WALL - x,
            TOP_WALL - y,
            y - BOTTOM_WALL};
}

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

void collide2(float &x1, float &y1, float m1, float &v_x1, float &v_y1, float &x2, float& y2, float m2, float &v_x2, float &v_y2)
{
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

class simulation
{
public:
    simulation(int n_balls, int seed): _balls(n_balls)//, _collided_1(n_balls), _collided_2(n_balls)
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

        region_id = aligned_alloc<int>(_balls.size);
        collisions_num = aligned_alloc<int>(_balls.size);

        for (int i = 0; i < MAX_COLLISIONS; i++)
            collisions[i] = aligned_alloc<int>(_balls.size);
    }

    ~simulation()
    {
        std::free(region_id);
        std::free(collisions_num);

        for (int i = 0; i < MAX_COLLISIONS; i++)
            std::free(collisions[i]);
    }

    void step()
    {
        #pragma omp parallel
        {
            #pragma omp for simd
            for (int i = 0; i < REGIONS_NUM; i++) {
                regions[i].size = 0;
            }

            #pragma omp for simd
            for (int i = 0; i <  _balls.size; i++) {
                region_id[i] = 0;
                collisions_num[i] = 0;
            }

            // int threadid = omp_get_thread_num();
            // int threadcount = omp_get_num_threads();
            // int num_per_thread = _balls.size / threadcount;

            // for (int i = 0; i < _balls.size; i++) {
            //     int region_x = (_balls.x[i] + WIDTH / 2) / REGIONS_SINGLE;
            //     int region_y = (_balls.y[i] + HEIGHT / 2) / REGIONS_SINGLE;

            //     int r = region_x + region_y * REGIONS_X;

            //     if (r >= threadid * num_per_thread && r < (threadid + 1) * num_per_thread)
            //         {
            //             auto *balls = &regions[r];
            //             balls->x[balls->size] = _balls.x[i];
            //             balls->y[balls->size] = _balls.y[i];
            //             balls->v_x[balls->size] = _balls.v_x[i];
            //             balls->v_y[balls->size] = _balls.v_y[i];
            //             balls->mass[balls->size] = _balls.mass[i];
            //             balls->index[balls->size] = _balls.index[i];
            //             balls->size++;
            //         }
            // }

            #pragma omp for
            for (int i = 0; i < _balls.size; i++) {
                Balls* balls;

                int region_x = (_balls.x[i] + WIDTH / 2) / REGIONS_SINGLE;
                int region_y = (_balls.y[i] + HEIGHT / 2) / REGIONS_SINGLE;

                int r = region_x + region_y * REGIONS_X;
                balls = &regions[r];

                int in;
                
                #pragma omp atomic capture
                in = balls->size++;

                balls->x[in] = _balls.x[i];
                balls->y[in] = _balls.y[i];
                balls->v_x[in] = _balls.v_x[i];
                balls->v_y[in] = _balls.v_y[i];
                balls->mass[in] = _balls.mass[i];
                balls->index[in] = _balls.index[i];
            }

            // #pragma omp for
            // for (int i = 0; i < _balls.size; i+= 8) {
            //     Balls* balls[8];

            //     for (int j = 0; j < 8; j++) {
            //         int region_x = (_balls.x[i + j] + WIDTH / 2) / REGIONS_SINGLE;
            //         int region_y = (_balls.y[i + j] + HEIGHT / 2) / REGIONS_SINGLE;

            //         int r = region_x + region_y * REGIONS_X;
            //         balls[j] = &regions[r];
            //     }

            //     #pragma omp critical
            //     {
            //         for (int j = 0; j < 8; j++) {
            //             balls[j]->x[balls[j]->size] = _balls.x[i + j];
            //             balls[j]->y[balls[j]->size] = _balls.y[i + j];
            //             balls[j]->v_x[balls[j]->size] = _balls.v_x[i+ j ];
            //             balls[j]->v_y[balls[j]->size] = _balls.v_y[i+ j];
            //             balls[j]->mass[balls[j]->size] = _balls.mass[i+ j];
            //             balls[j]->index[balls[j]->size] = _balls.index[i+ j];
            //             balls[j]->size++;
            //         }
            //     }
            // }

            // #pragma omp for
            // for (int r = 0; r < REGIONS_NUM; r++) {
            //     int neighbour_regions[] = {1, REGIONS_X, REGIONS_X + 1};

            //     // neighbour regions
            //     for (int rs = 0; rs < 3; rs++) {
            //         int neighbour = r + neighbour_regions[rs];

            //         if (neighbour < 0 || neighbour >= REGIONS_NUM)
            //             continue;

            //         for (int i = 0; i < regions[r].size; i++) {
            //             for (int j = 0; j < regions[neighbour].size; j++) {
            //                 if (is_collision(regions[r].x[i],
            //                         regions[r].y[i],
            //                         regions[r].mass[i],
            //                         regions[neighbour].x[j],
            //                         regions[neighbour].y[j],
            //                         regions[neighbour].mass[j])) {
            //                             auto &num = collisions_num[regions[r].index[i]];
            //                             collisions[num][regions[r].index[i]] = regions[neighbour].index[j];
            //                             num++;
            //                         }
            //             }
            //         }
            //     }

            //     for (int i = 0; i < regions[r].size; i++) {
            //         for (int j = i + 1; j < regions[r].size; j++) {
            //                 if (is_collision(regions[r].x[i],
            //                         regions[r].y[i],
            //                         regions[r].mass[i],
            //                         regions[r].x[j],
            //                         regions[r].y[j],
            //                         regions[r].mass[j])) {
            //                             auto &num = collisions_num[regions[r].index[i]];
            //                             collisions[num][regions[r].index[i]] = regions[r].index[j];
            //                             num++;
            //                 }
            //         }
            //     }
            // }

            // for (int i = 0; i < MAX_COLLISIONS; i++)
            // {
            //     Balls collided_1(_balls.size);
            //     Balls collided_2(_balls.size);
            //     int size = 0;

            //     #pragma omp for
            //     for (int j = 0; j < _balls.size; j++) {
            //         if (collisions_num[j] > i) {
            //             collided_1.x[size] = _balls.x[j];
            //             collided_1.y[size] = _balls.y[j];
            //             collided_1.v_x[size] = _balls.v_x[j];
            //             collided_1.v_y[size] = _balls.v_y[j];
            //             collided_1.mass[size] = _balls.mass[j];
            //             collided_1.index[size] = _balls.index[j];

            //             collided_2.x[size] = _balls.x[collisions[i][j]];
            //             collided_2.y[size] = _balls.y[collisions[i][j]];
            //             collided_2.v_x[size] = _balls.v_x[collisions[i][j]];
            //             collided_2.v_y[size] = _balls.v_y[collisions[i][j]];
            //             collided_2.mass[size] = _balls.mass[collisions[i][j]];
            //             collided_2.index[size] = _balls.index[collisions[i][j]];
            //             size++;
            //         }
            //     }

            //     #pragma omp for
            //     for (int j = 0; j < size; j++) {
            //             collide2(collided_1.x[j],
            //                         collided_1.y[j],
            //                         collided_1.mass[j],
            //                         collided_1.v_x[j],
            //                         collided_1.v_y[j],
            //                         collided_2.x[j],
            //                         collided_2.y[j],
            //                         collided_2.mass[j],
            //                         collided_2.v_x[j],
            //                         collided_2.v_y[j]);
            //     }

            //     #pragma omp for
            //     for (int j = 0; j < size; j++) {
            //             auto index = collided_1.index[j];

            //             _balls.x[index] = collided_1.x[j];
            //             _balls.y[index] = collided_1.y[j];
            //             _balls.v_x[index] = collided_1.v_x[j];
            //             _balls.v_y[index] = collided_1.v_y[j];
            //             _balls.mass[index] = collided_1.mass[j];
            //             _balls.index[index] = collided_1.index[j];

            //             index = collided_2.index[j];

            //             _balls.x[index] = collided_2.x[j];
            //             _balls.y[index] = collided_2.y[j];
            //             _balls.v_x[index] = collided_2.v_x[j];
            //             _balls.v_y[index] = collided_2.v_y[j];
            //             _balls.mass[index] = collided_2.mass[j];
            //             _balls.index[index] = collided_2.index[j];
            //     }
            // }

            #pragma omp for collapse(2)
            for (int r1 = 0; r1 < REGIONS_X; r1 +=2) {
                for (int r2 = 0; r2 < REGIONS_X; r2 +=2) {
                    do_collide(r1 * REGIONS_X + r2);
                }
            }

            #pragma omp for collapse(2)
            for (int r1 = 1; r1 < REGIONS_X; r1 +=2) {
                for (int r2 = 0; r2 < REGIONS_X; r2 +=2) {
                    do_collide(r1 * REGIONS_X + r2);
                }
            }

            #pragma omp for collapse(2)
            for (int r1 = 0; r1 < REGIONS_X; r1 +=2) {
                for (int r2 = 1; r2 < REGIONS_X; r2 +=2) {
                    do_collide(r1 * REGIONS_X + r2);
                }
            }

            #pragma omp for collapse(2)
            for (int r1 = 1; r1 < REGIONS_X; r1 += 2) {
                for (int r2 = 1; r2 < REGIONS_X; r2 +=2) {
                    do_collide(r1 * REGIONS_X + r2);
                }
            }

            // assign data back
            #pragma omp for
            for (int r = 0; r < REGIONS_NUM; r++) {
                for (int i = 0; i < regions[r].size; i++) {
                    auto &balls = regions[r];
                    auto index = balls.index[i];

                    _balls.x[index] = balls.x[i];
                    _balls.y[index] = balls.y[i];
                    _balls.v_x[index] = balls.v_x[i];
                    _balls.v_y[index] = balls.v_y[i];
                    _balls.mass[index] = balls.mass[i];
                    _balls.index[index] = balls.index[i];
                }   
            }

            // collide with walls
            #pragma omp for
            for (int i = 0; i < _balls.size; i++) {
                collide_wall(_balls.x[i],
                        _balls.y[i],
                        _balls.mass[i],
                        _balls.v_x[i],
                        _balls.v_y[i]);
            }

            #pragma omp for simd
            for (int i = 0; i < _balls.size; i++) {
                    _balls.x[i] += _balls.v_x[i];
                    _balls.y[i] += _balls.v_y[i];
            }
        }
    }

    const Balls& balls()
    {
        return _balls;
    }

private:
    void do_collide(int r)
    {
        int neighbour_regions[] = {1, REGIONS_X, REGIONS_X + 1};
                // neighbour regions
                for (int rs = 0; rs < 3; rs++) {
                    int neighbour = r + neighbour_regions[rs];

                    if (neighbour < 0 || neighbour >= REGIONS_NUM)
                        continue;

                    for (int i = 0; i < regions[r].size; i++) {
                        for (int j = 0; j < regions[neighbour].size; j++) {
                            collide(regions[r].x[i],
                                    regions[r].y[i],
                                    regions[r].mass[i],
                                    regions[r].v_x[i],
                                    regions[r].v_y[i],
                                    regions[neighbour].x[j],
                                    regions[neighbour].y[j],
                                    regions[neighbour].mass[j],
                                    regions[neighbour].v_x[j],
                                    regions[neighbour].v_y[j]);
                        }
                    }
                }

                // same region
                for (int i = 0; i < regions[r].size; i++) {
                    for (int j = i + 1; j < regions[r].size; j++) {
                            collide(regions[r].x[i],
                                    regions[r].y[i],
                                    regions[r].mass[i],
                                    regions[r].v_x[i],
                                    regions[r].v_y[i],
                                    regions[r].x[j],
                                    regions[r].y[j],
                                    regions[r].mass[j],
                                    regions[r].v_x[j],
                                    regions[r].v_y[j]);
                    }
                }
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

    std::array<Balls, REGIONS_NUM> regions;

    int *region_id;
    int *collisions[MAX_COLLISIONS];
    int *collisions_num;
    int x_regions, y_regions;
    Balls _balls;
   // Balls _collided_1, _collided_2;
};
