#pragma once

#define WIDTH 640
#define HEIGHT 480
#define MARGIN 10
#define LEFT_WALL (- WIDTH / 2 + MARGIN)
#define RIGHT_WALL (WIDTH / 2 - MARGIN)
#define BOTTOM_WALL (- HEIGHT / 2 + MARGIN)
#define TOP_WALL (HEIGHT / 2 - MARGIN)

static constexpr float OVERLAP_MARGIN = 0.01;

static constexpr float MAX_MASS = 15; // MASS == RADIUS
static constexpr float MIN_MASS  = 5;

static constexpr float MAX_SPEED = 10;
static constexpr float MIN_SPEED = 1;
