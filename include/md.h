#ifndef MD_H
#define MD_H

#include <string>
#include <vector>

#define MD_VERSION "TO.2.3"

constexpr double kPI = 3.14159265358979323846;

// ----------------------
// Detection of one animal at one frame
// ----------------------
struct Detection {
    int ID;
    long custom_frame;
    double cen_x;
    double cen_y;
    double dir_x;
    double dir_y;
    double angle;
    double ahead_density; // density ahead
    double local_density; // density at current position
    double speed;
    double acceleration;
    double angle_change;
    double distance_from_last;
    bool interpolated = false;
};

// ----------------------
// Segment of continuous tracking for one animal
// ----------------------
struct Segment {
    std::vector<Detection> points;
};

// ----------------------
// Track summary over all segments
// ----------------------
struct TrackSummary {
    int ID;
    int segment_index;
    long first_frame;
    long last_frame;
    int n_obs;
    double total_distance;
    double mean_speed;
    long max_gap;
    double mean_angle_change;
    double max_jump;
};

// ----------------------
// Main engine entry point
// ----------------------
int motion_dynamics_run(const std::string& input,
    const std::string& output_prefix,
    long frame_window,
    bool smooth = false); // optional smoothing

#endif