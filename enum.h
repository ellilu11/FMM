#pragma once

enum class Mode {
    READ,
    GEN
};

enum class Dist {
    UNIFORM,
    GAUSSIAN,
    GRID
};

enum class ChargeDist {
    PLUS,
    MINUS,
    DIP,
    QUAD,
    OCT,
    RAND
};

enum class Precision {
    LOW,
    MEDIUM,
    HIGH
};

enum class Dir {
    W, E, S, N, D, U,
    SW, SE, NW, NE, DW, DE, UW, UE, DS, DN, US, UN,
    DSW, DSE, DNW, DNE, USW, USE, UNW, UNE
};
