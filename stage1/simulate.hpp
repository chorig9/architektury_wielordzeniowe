#pragma once

#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "constants.hpp"

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
    Ball(vector position, vector velocity, float mass): position(position), velocity(velocity), mass(mass)
    {
    }

    vector position;
    vector velocity;
    const float mass;

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

/* XXX-opt: calculate everything and one if later? */
/* return table of distances to walls and call std::min */
vector wall_collision_point(Ball &lhs)
{
    if (lhs.position.x + lhs.radius() >= RIGHT_WALL)
        return {RIGHT_WALL, lhs.position.y};
    if (lhs.position.x - lhs.radius() <= LEFT_WALL)
        return {LEFT_WALL, lhs.position.y};
    if (lhs.position.y + lhs.radius() >= TOP_WALL)
        return {lhs.position.x, TOP_WALL};
    if (lhs.position.y - lhs.radius() <= BOTTOM_WALL)
        return {lhs.position.x, BOTTOM_WALL};
    else
        return {0, 0};
}

void collide_wall(Ball &lhs)
{
    auto &v1 = lhs.velocity;
    auto &m1 = lhs.mass;
    auto &x1 = lhs.position;
    auto x2 = wall_collision_point(lhs);

    /* no collision */
    if (x2.x == 0 && x2.y == 0)
        return;

    

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
    x2.x += OVERLAP_MARGIN;
    x2.y += OVERLAP_MARGIN;

    auto norm1 = norm_pow2(x1 - x2);
    auto norm2 = norm_pow2(x2 - x1);

    auto new_v1 = v1 - (2 * m2) / (m1 + m2) * dot(v1 - v2, x1 - x2) / norm1 * (x1 - x2);

    v2 -= (2 * m1) / (m1 + m2) * dot(v2 - v1, x2 - x1) / norm2 * (x2 - x1);
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

            _balls.emplace_back(position, velocity, mass);
        }
    }

    void step()
    {
        for (int i = 0; i < _balls.size() ; i++) {
            for (int j = i + 1; j < _balls.size(); j++) {
                collide(_balls[i], _balls[j]);
            }

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
