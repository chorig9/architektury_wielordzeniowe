#pragma once

#define WIDTH 640
#define HEIGHT 480
#define MARGIN 10
#define LEFT_WALL (- WIDTH / 2 + MARGIN)
#define RIGHT_WALL (WIDTH / 2 - MARGIN)
#define BOTTOM_WALL (- HEIGHT / 2 + MARGIN)
#define TOP_WALL (HEIGHT / 2 - MARGIN)

static constexpr float SIMULATION_STEP = 0.1;

static constexpr float MAX_MASS = 6; // MASS == RADIUS
static constexpr float MIN_MASS  = 3;

static constexpr float MAX_SPEED = 3;
static constexpr float MIN_SPEED = 1;
