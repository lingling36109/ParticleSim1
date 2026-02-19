#include "common.h"
#include <cmath>
#include <vector>

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

using Bin = std::vector<particle_t*>;
std::vector<std::vector<Bin>> bins;
constexpr double BINSIZE = cutoff;
int num_bins;

void helper(particle_t* parts, int num_parts) {
    // Clear bins
    for (int i = 0; i < num_bins; ++i)
        for (int j = 0; j < num_bins; ++j)
            bins[i][j].clear();

    // Insert particles
    for (int i = 0; i < num_parts; ++i) {
        int bx = static_cast<int>(parts[i].x / BINSIZE);
        int by = static_cast<int>(parts[i].y / BINSIZE);

        bx = std::min(bx, num_bins - 1);
        by = std::min(by, num_bins - 1);

        bins[bx][by].push_back(&parts[i]);
    }
}


void init_simulation(particle_t* parts, int num_parts, double size) {
    num_bins = static_cast<int>(std::floor(size / BINSIZE));
    bins.assign(num_bins, std::vector<Bin>(num_bins));

    helper(parts, num_parts);
}

void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Reset accelerations
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;
    }

    // Compute forces using bins
    for (int bx = 0; bx < num_bins; ++bx) {
        for (int by = 0; by < num_bins; ++by) {
            for (particle_t* p : bins[bx][by]) {
                for (int dx = -1; dx <= 1; ++dx) {
                    for (int dy = -1; dy <= 1; ++dy) {

                        int nbx = bx + dx;
                        int nby = by + dy;

                        if (nbx < 0 || nbx >= num_bins ||
                            nby < 0 || nby >= num_bins)
                            continue;

                        for (particle_t* q : bins[nbx][nby]) {
                            apply_force(*p, *q);
                        }
                    }
                }
            }
        }
    }

    // Move particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }

    // Rebuild bins
    helper(parts, num_parts);
}
