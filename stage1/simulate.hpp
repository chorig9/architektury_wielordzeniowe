#pragma once

#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "constants.hpp"
#include <array>

/*
 XXX:
 1. data-oriented
 2. split screen + assign the same balls to same threads
 */

struct vector
{
    float x;
    float y;

    vector& operator-=(const vector &rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    vector& operator+=(const vector &rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    vector operator+(const vector &rhs)
    {
        return {x + rhs.x, y + rhs.y};
    }

    vector operator-(const vector &rhs)
    {
        return {x - rhs.x, y - rhs.y};
    }

    friend vector operator*(float lhs, const vector &rhs)
    {
        return {rhs.x * lhs, rhs.y * lhs};
    }
};

struct Ball
{
    Ball(vector position, vector velocity, float mass, int i): position(position), velocity(velocity), mass(mass), index(i)
    {
    }

    Ball(const Ball &rhs) = default;

    Ball &operator=(const Ball &rhs) = default;

    vector position;
    vector velocity;
    float mass;

    int index;

    const float radius() const
    {
        return mass;
    }
};

float dot(const vector &lhs, const vector &rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

float norm_pow2(const vector &rhs)
{
    return rhs.x * rhs.x + rhs.y * rhs.y;        
}

bool is_collision(Ball &lhs, Ball &rhs)
{
    float xd = lhs.position.x - rhs.position.x;
    float yd = lhs.position.y - rhs.position.y;

    float sumRadius = lhs.radius() + rhs.radius();
    float sqrRadius = sumRadius * sumRadius;

    float distSqr = (xd * xd) + (yd * yd);

    return distSqr <= sqrRadius; 
}

struct wall_distance
{
    float l_wall, r_wall, t_wall, b_wall;
};

wall_distance wall_collision_point(Ball &lhs)
{
    return {lhs.position.x - LEFT_WALL,
            RIGHT_WALL - lhs.position.x,
            TOP_WALL - lhs.position.y,
            lhs.position.y - BOTTOM_WALL};
}

void collide_wall(Ball &lhs)
{
    auto wall_dists = wall_collision_point(lhs);
    vector x2 = {lhs.position.x, lhs.position.y};

    if (wall_dists.l_wall < lhs.radius()) {
        lhs.position.x = LEFT_WALL + lhs.radius();
        x2.x = LEFT_WALL;
    } else if (wall_dists.r_wall < lhs.radius()) {
        lhs.position.x = RIGHT_WALL - lhs.radius();
        x2.x = RIGHT_WALL;
    /* only consider collision with one wall */
    } else if (wall_dists.t_wall < lhs.radius()) {
        lhs.position.y = TOP_WALL - lhs.radius();
        x2.y = TOP_WALL;
    } else if (wall_dists.b_wall < lhs.radius()) {
        lhs.position.y = BOTTOM_WALL + lhs.radius();
        x2.y = BOTTOM_WALL;
    } else {
        return;
    }

    auto &v1 = lhs.velocity;
    auto &m1 = lhs.mass;
    auto &x1 = lhs.position;

    auto norm1 = norm_pow2(x1 - x2);

    v1 -= 2 * dot(v1, x1 - x2) / norm1 * (x1 - x2);
}

void collide(Ball &lhs, Ball &rhs)
{
    if (!is_collision(lhs, rhs))
        return;

    auto &v1 = lhs.velocity;
    auto &v2 = rhs.velocity;
    auto &m1 = lhs.mass;
    auto &m2 = rhs.mass;
    auto &x1 = lhs.position;
    auto &x2 = rhs.position;

    vector overlap = {x2.x - x1.x, x2.y - x1.y};
    float overlap_dist = sqrt(overlap.x * overlap.x + overlap.y * overlap.y);
    float radius_sum = lhs.radius() + rhs.radius();
    float coef = radius_sum / overlap_dist;

    /* move x2 away from x1 so that x2 - x1 == r1 + r2 */
    x2 = x1 + coef * overlap;

    auto norm = norm_pow2(x1 - x2);

    auto new_v1 = v1 - (2 * m2) / (m1 + m2) * dot(v1 - v2, x1 - x2) / norm * (x1 - x2);

    v2 -= (2 * m1) / (m1 + m2) * dot(v2 - v1, x2 - x1) / norm * (x2 - x1);
    v1 = new_v1;
}

void advance(Ball &b)
{
    b.position += b.velocity;
}

class simulation
{
public:
    simulation(int n_balls, int seed): _balls()
    {
        _balls.reserve(n_balls);
        srand(seed);

        x_regions = (WIDTH - 2 * MARGIN - 2 * MAX_MASS) / (2 * MAX_MASS);
        y_regions = (HEIGHT - 2 * MARGIN - 2 * MAX_MASS) / (2 * MAX_MASS);

        int n_regions = x_regions * y_regions;
        std::vector<bool> pos_bitmap(n_regions, false);

        for (int i = 0; i < n_balls; i++) {
            auto position = get_random_pos(pos_bitmap);
            auto velocity = vector{get_random(MIN_SPEED, MAX_SPEED), get_random(MIN_SPEED, MAX_SPEED)};
            auto mass = get_random(MIN_MASS, MAX_MASS);

            _balls.emplace_back(position, velocity, mass, i);
        }
    }

    void step()
    {
        std::array<std::vector<Ball>, REGIONS_NUM> regions;

        // for (int r = 0; r < REGIONS_NUM; r++) {
        //     // XXX: just calculate to which region ball should go? (in parralell)
        //     for (int i = 0; i < _balls.size(); i++) {
        //         int region_x = (_balls[i].position.x + WIDTH / 2) / REGIONS_SINGLE;
        //         int region_y = (_balls[i].position.y + HEIGHT / 2) / REGIONS_SINGLE;

        //         if (region_x == (r % REGIONS_X) && region_y == (r / REGIONS_X)) {
        //             regions[r].push_back(_balls[i]);
        //         }
        //     }
        // }

        for (int i = 0; i < _balls.size(); i++) {
                int region_x = (_balls[i].position.x + WIDTH / 2) / REGIONS_SINGLE;
                int region_y = (_balls[i].position.y + HEIGHT / 2) / REGIONS_SINGLE;

                int r = region_x + region_y * REGIONS_X;
                regions[r].push_back(_balls[i]);
            }

        for (int r = 0; r < REGIONS_NUM; r++) {

            int neighbour_regions[] = {1, REGIONS_X, REGIONS_X + 1};

            for (int rs = 0; rs < 3; rs++) {
                int neighbour = r + neighbour_regions[rs];

                if (neighbour < 0 || neighbour >= REGIONS_NUM)
                    continue;

                for (int i = 0; i < regions[r].size(); i++) {
                    for (int j = 0; j < regions[neighbour].size(); j++) {
                        collide(regions[r][i], regions[neighbour][j]);
                    }
                }
            }

            for (int i = 0; i < regions[r].size(); i++) {
                for (int j = i + 1; j < regions[r].size(); j++) {
                    collide(regions[r][i], regions[r][j]);
                }
            }
        }

        for (int r = 0; r < REGIONS_NUM; r++) {
            for (int i = 0; i < regions[r].size(); i++) {
                _balls[regions[r][i].index] = regions[r][i];
            }   
        } 

        for (int i = 0; i < _balls.size() ; i++) {
            collide_wall(_balls[i]);

            advance(_balls[i]);
        }
    }

    const std::vector<Ball>& balls()
    {
        return _balls;
    }

private:
    float get_random(float lo, float hi)
    {
        return lo + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(hi-lo)));
    }

    vector get_random_pos(std::vector<bool> &pos_bitmap)
    {
            int pos_random;
            do {
                pos_random = rand() % pos_bitmap.size();
            } while (pos_bitmap[pos_random] == true);

            pos_bitmap[pos_random] = true;

            int x = (pos_random % x_regions) * 2 * MAX_MASS + 2 * MAX_MASS - WIDTH / 2; 
            int y = (pos_random / x_regions) * 2 * MAX_MASS + 2 * MAX_MASS - HEIGHT / 2;

            return {x, y};
    }

    int x_regions, y_regions;
    std::vector<Ball> _balls;
};
