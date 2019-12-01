#pragma once

#define WIDTH 1024
#define HEIGHT 720
#define MARGIN 10
#define LEFT_WALL (- WIDTH / 2 + MARGIN)
#define RIGHT_WALL (WIDTH / 2 - MARGIN)
#define BOTTOM_WALL (- HEIGHT / 2 + MARGIN)
#define TOP_WALL (HEIGHT / 2 - MARGIN)

static constexpr int DRAW_STEP = 1;

static constexpr float MAX_MASS = 6; // MASS == RADIUS
static constexpr float MIN_MASS  = 3;

static constexpr float MAX_SPEED = 1;
static constexpr float MIN_SPEED = 0.1;

static constexpr int REGION_WIDTH = 70;
static constexpr int REGION_OVERLAP = MAX_MASS;
static constexpr int REGIONS_SINGLE = REGION_WIDTH - REGION_OVERLAP;
static constexpr int REGIONS_X = (WIDTH - 2 * MARGIN + REGIONS_SINGLE) / REGIONS_SINGLE;
static constexpr int REGIONS_NUM = REGIONS_X * REGIONS_X;
